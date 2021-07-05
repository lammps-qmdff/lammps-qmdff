/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "npair_half_bin_newton_tri_intel.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "comm.h"
#include "group.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalfBinNewtonTriIntel::NPairHalfBinNewtonTriIntel(LAMMPS *lmp) : 
  NPairIntel(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction with Newton's 3rd law for triclinic
   each owned atom i checks its own bin and other bins in triclinic stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void NPairHalfBinNewtonTriIntel::build(NeighList *list)
{
  if (nstencil > INTEL_MAX_STENCIL)
    error->all(FLERR, "Too many neighbor bins for USER-INTEL package.");

  #ifdef _LMP_INTEL_OFFLOAD
  if (exclude)
    error->all(FLERR, "Exclusion lists not yet supported for Intel offload");
  #endif

  if (_fix->precision() == FixIntel::PREC_MODE_MIXED)
    hbnti(list, _fix->get_mixed_buffers());
  else if (_fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    hbnti(list, _fix->get_double_buffers());
  else
    hbnti(list, _fix->get_single_buffers());

  _fix->stop_watch(TIME_HOST_NEIGHBOR);
}

template <class flt_t, class acc_t>
void NPairHalfBinNewtonTriIntel::
hbnti(NeighList *list, IntelBuffers<flt_t,acc_t> *buffers) {
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  list->inum = nlocal;

  int host_start = _fix->host_start_neighbor();
  const int off_end = _fix->offload_end_neighbor();

  #ifdef _LMP_INTEL_OFFLOAD
  if (off_end) grow_stencil();
  if (_fix->full_host_list()) host_start = 0;
  int offload_noghost = _fix->offload_noghost();
  #endif

  buffers->grow_list(list, atom->nlocal, comm->nthreads, off_end);

  int need_ic = 0;
  if (atom->molecular)
    dminimum_image_check(need_ic, neighbor->cutneighmax, neighbor->cutneighmax,
			 neighbor->cutneighmax);

  #ifdef _LMP_INTEL_OFFLOAD
  if (need_ic) {
    if (offload_noghost) {
      hbnti<flt_t,acc_t,1,1>(1, list, buffers, 0, off_end);
      hbnti<flt_t,acc_t,1,1>(0, list, buffers, host_start, nlocal, off_end);
    } else {
      hbnti<flt_t,acc_t,0,1>(1, list, buffers, 0, off_end);
      hbnti<flt_t,acc_t,0,1>(0, list, buffers, host_start, nlocal);
    }
  } else {
    if (offload_noghost) {
      hbnti<flt_t,acc_t,1,0>(1, list, buffers, 0, off_end);
      hbnti<flt_t,acc_t,1,0>(0, list, buffers, host_start, nlocal, off_end);
    } else {
      hbnti<flt_t,acc_t,0,0>(1, list, buffers, 0, off_end);
      hbnti<flt_t,acc_t,0,0>(0, list, buffers, host_start, nlocal);
    }
  }
  #else
  if (need_ic)
    hbnti<flt_t,acc_t,0,1>(0, list, buffers, host_start, nlocal);
  else
    hbnti<flt_t,acc_t,0,0>(0, list, buffers, host_start, nlocal);
  #endif
}

template <class flt_t, class acc_t, int offload_noghost, int need_ic>
void NPairHalfBinNewtonTriIntel::
hbnti(const int offload, NeighList *list, IntelBuffers<flt_t,acc_t> *buffers,
      const int astart, const int aend, const int offload_end) {
  if (aend-astart == 0) return;

  const int nall = atom->nlocal + atom->nghost;
  int pad = 1;
  int nall_t = nall;

  #ifdef _LMP_INTEL_OFFLOAD
  if (offload_noghost && offload) nall_t = atom->nlocal;
  if (offload) {
    if (INTEL_MIC_NBOR_PAD > 1)
      pad = INTEL_MIC_NBOR_PAD * sizeof(float) / sizeof(flt_t);
  } else
  #endif
    if (INTEL_NBOR_PAD > 1)
      pad = INTEL_NBOR_PAD * sizeof(float) / sizeof(flt_t);
  const int pad_width = pad;

  const ATOM_T * _noalias const x = buffers->get_x();
  int * _noalias const firstneigh = buffers->firstneigh(list);
  const int e_nall = nall_t;

  const int molecular = atom->molecular;
  int *ns = NULL;
  tagint *s = NULL;
  int tag_size = 0, special_size;
  if (buffers->need_tag()) tag_size = e_nall;
  if (molecular) {
    s = atom->special[0];
    ns = atom->nspecial[0];
    special_size = aend;
  } else {
    s = &buffers->_special_holder;
    ns = &buffers->_nspecial_holder;
    special_size = 0;
  }
  const tagint * _noalias const special = s;
  const int * _noalias const nspecial = ns;
  const int maxspecial = atom->maxspecial;
  const tagint * _noalias const tag = atom->tag;

  int * _noalias const ilist = list->ilist;
  int * _noalias numneigh = list->numneigh;
  int * _noalias const cnumneigh = buffers->cnumneigh(list);
  const int nstencil = this->nstencil;
  const int * _noalias const stencil = this->stencil;
  const flt_t * _noalias const cutneighsq = buffers->get_cutneighsq()[0];
  const int ntypes = atom->ntypes + 1;
  const int nlocal = atom->nlocal;

  #ifndef _LMP_INTEL_OFFLOAD
  int * const mask = atom->mask;
  tagint * const molecule = atom->molecule;
  #endif

  int tnum;
  int *overflow;
  double *timer_compute;
  #ifdef _LMP_INTEL_OFFLOAD
  if (offload) {
    timer_compute = _fix->off_watch_neighbor();
    tnum = buffers->get_off_threads();
    overflow = _fix->get_off_overflow_flag();
    _fix->stop_watch(TIME_HOST_NEIGHBOR);
    _fix->start_watch(TIME_OFFLOAD_LATENCY);
  } else 
  #endif
  {
    tnum = comm->nthreads;
    overflow = _fix->get_overflow_flag();
  }
  const int nthreads = tnum;
  const int maxnbors = buffers->get_max_nbors();
  int * _noalias const atombin = buffers->get_atombin();
  const int * _noalias const binpacked = buffers->get_binpacked();

  const int xperiodic = domain->xperiodic;
  const int yperiodic = domain->yperiodic;
  const int zperiodic = domain->zperiodic;
  const flt_t xprd_half = domain->xprd_half;
  const flt_t yprd_half = domain->yprd_half;
  const flt_t zprd_half = domain->zprd_half;

  #ifdef _LMP_INTEL_OFFLOAD
  const int * _noalias const binhead = this->binhead;
  const int * _noalias const bins = this->bins;
  const int cop = _fix->coprocessor_number();
  const int separate_buffers = _fix->separate_buffers();
  #pragma offload target(mic:cop) if(offload) \
    in(x:length(e_nall+1) alloc_if(0) free_if(0)) \
    in(tag:length(tag_size) alloc_if(0) free_if(0)) \
    in(special:length(special_size*maxspecial) alloc_if(0) free_if(0)) \
    in(nspecial:length(special_size*3) alloc_if(0) free_if(0)) \
    in(bins,binpacked:length(nall) alloc_if(0) free_if(0)) \
    in(binhead:length(mbins+1) alloc_if(0) free_if(0)) \
    in(cutneighsq:length(0) alloc_if(0) free_if(0)) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(cnumneigh:length(0) alloc_if(0) free_if(0)) \
    out(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(ilist:length(0) alloc_if(0) free_if(0)) \
    in(atombin:length(aend) alloc_if(0) free_if(0)) \
    in(stencil:length(nstencil) alloc_if(0) free_if(0)) \
    in(maxnbors,nthreads,maxspecial,nstencil,offload_end,pad_width,e_nall) \
    in(offload,separate_buffers, astart, aend, nlocal, molecular, ntypes) \
    in(xperiodic, yperiodic, zperiodic, xprd_half, yprd_half, zprd_half) \
    out(overflow:length(5) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(tag)
  #endif
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime();
    #endif

    #ifdef _LMP_INTEL_OFFLOAD
    overflow[LMP_LOCAL_MIN] = astart;
    overflow[LMP_LOCAL_MAX] = aend - 1;
    overflow[LMP_GHOST_MIN] = e_nall;
    overflow[LMP_GHOST_MAX] = -1;
    #endif

    int nstencilp = 0;
    int binstart[INTEL_MAX_STENCIL], binend[INTEL_MAX_STENCIL];
    for (int k = 0; k < nstencil; k++) {
      binstart[nstencilp] = stencil[k];
      int end = stencil[k] + 1;
      for (int kk = k + 1; kk < nstencil; kk++) {
        if (stencil[kk-1]+1 == stencil[kk]) {
          end++;
          k++;
        } else break;
      }
      binend[nstencilp] = end;
      nstencilp++;
    }

    #if defined(_OPENMP)
    #pragma omp parallel default(none) \
      shared(numneigh, overflow, nstencilp, binstart, binend)
    #endif
    {
      #ifdef _LMP_INTEL_OFFLOAD
      int lmin = e_nall, lmax = -1, gmin = e_nall, gmax = -1;
      #endif

      const int num = aend - astart;
      int tid, ifrom, ito;
      IP_PRE_omp_range_id(ifrom, ito, tid, num, nthreads);
      ifrom += astart;
      ito += astart;

      int which;

      const int list_size = (ito + tid * 2 + 2) * maxnbors;
      int ct = (ifrom + tid * 2) * maxnbors;
      int *neighptr = firstneigh + ct;
      const int obound = maxnbors * 3;

      for (int i = ifrom; i < ito; i++) {
        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        const int itype = x[i].w;
        const int ioffset = ntypes * itype;

        // loop over all atoms in bins in stencil
        // pairs for atoms j "below" i are excluded
        // below = lower z or (equal z and lower y) or (equal zy and lower x)
        //         (equal zyx and j <= i)
        // latter excludes self-self interaction but allows superposed atoms

        const int ibin = atombin[i];

	int raw_count = maxnbors;
        for (int k = 0; k < nstencilp; k++) {
          const int bstart = binhead[ibin + binstart[k]];
          const int bend = binhead[ibin + binend[k]];
          for (int jj = bstart; jj < bend; jj++) {
            const int j = binpacked[jj];

	    #ifdef _LMP_INTEL_OFFLOAD
	    if (offload_noghost) {
              if (j < nlocal) {
                if (i < offload_end) continue;
              } else if (offload) continue;
            }
	    #endif

            if (x[j].z < ztmp) continue;
            if (x[j].z == ztmp) {
              if (x[j].y < ytmp) continue;
              if (x[j].y == ytmp) {
                if (x[j].x < xtmp) continue;
                if (x[j].x == xtmp && j <= i) continue;
              }
            }

            #ifndef _LMP_INTEL_OFFLOAD
            if (exclude) {
	      const int jtype = x[j].w;
	      if (exclusion(i,j,itype,jtype,mask,molecule)) continue;
	    }
	    #endif

	    neighptr[raw_count++] = j;
	  }
	}
	if (raw_count > obound)
	  *overflow = 1;

        #if defined(LMP_SIMD_COMPILER)
	#ifdef _LMP_INTEL_OFFLOAD
	int vlmin = lmin, vlmax = lmax, vgmin = gmin, vgmax = gmax;
	#if __INTEL_COMPILER+0 > 1499
	#pragma vector aligned
        #pragma simd reduction(max:vlmax,vgmax) reduction(min:vlmin, vgmin)
	#endif
	#else
	#pragma vector aligned
        #pragma simd
	#endif
	#endif
	for (int u = maxnbors; u < raw_count; u++) {
          int j = neighptr[u];
          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
	  const int jtype = x[j].w;
          const flt_t rsq = delx * delx + dely * dely + delz * delz;
          if (rsq > cutneighsq[ioffset + jtype]) 
	    neighptr[u] = e_nall;
	  else {
            if (need_ic) {
	      int no_special;
	      ominimum_image_check(no_special, delx, dely, delz);
	      if (no_special)
	        neighptr[u] = -j - 1;
	    }

            #ifdef _LMP_INTEL_OFFLOAD
            if (j < nlocal) {
              if (j < vlmin) vlmin = j;
              if (j > vlmax) vlmax = j;
            } else {
              if (j < vgmin) vgmin = j;
              if (j > vgmax) vgmax = j;
            }
            #endif
          }
        }

        int n = 0, n2 = maxnbors;
	for (int u = maxnbors; u < raw_count; u++) {
	  const int j = neighptr[u];
	  int pj = j;
	  if (pj < e_nall) {
	    if (need_ic)
	      if (pj < 0) pj = -pj - 1;

	    if (pj < nlocal)
	      neighptr[n++] = j;
	    else
	      neighptr[n2++] = j;
	  }
	}
	int ns = n;
	for (int u = maxnbors; u < n2; u++)
	  neighptr[n++] = neighptr[u];

        ilist[i] = i;
        cnumneigh[i] = ct;
	ns += n2 - maxnbors;

        int edge = (ns % pad_width);
        if (edge) {
          const int pad_end = ns + (pad_width - edge);
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min=1, max=15, avg=8
          #endif
          for ( ; ns < pad_end; ns++)
            neighptr[ns] = e_nall;
        }
        numneigh[i] = ns;

	ct += ns;
	const int alignb = (INTEL_DATA_ALIGN / sizeof(int));
	edge = (ct % alignb);
	if (edge) ct += alignb - edge;
	neighptr = firstneigh + ct;
	if (ct + obound > list_size) {
	  if (i < ito - 1) {
	    *overflow = 1;
	    ct = (ifrom + tid * 2) * maxnbors;
	  }
	}
      }

      if (*overflow == 1)
        for (int i = ifrom; i < ito; i++)
          numneigh[i] = 0;

      #ifdef _LMP_INTEL_OFFLOAD
      if (separate_buffers) {
        #if defined(_OPENMP)
        #pragma omp critical
        #endif
        {
          if (lmin < overflow[LMP_LOCAL_MIN]) overflow[LMP_LOCAL_MIN] = lmin;
          if (lmax > overflow[LMP_LOCAL_MAX]) overflow[LMP_LOCAL_MAX] = lmax;
          if (gmin < overflow[LMP_GHOST_MIN]) overflow[LMP_GHOST_MIN] = gmin;
          if (gmax > overflow[LMP_GHOST_MAX]) overflow[LMP_GHOST_MAX] = gmax;
        }
        #pragma omp barrier
      }

      int ghost_offset = 0, nall_offset = e_nall;
      if (separate_buffers) {
        int nghost = overflow[LMP_GHOST_MAX] + 1 - overflow[LMP_GHOST_MIN];
        if (nghost < 0) nghost = 0;
        if (offload) {
          ghost_offset = overflow[LMP_GHOST_MIN] - overflow[LMP_LOCAL_MAX] - 1;
          nall_offset = overflow[LMP_LOCAL_MAX] + 1 + nghost;
	} else {
          ghost_offset = overflow[LMP_GHOST_MIN] - nlocal;
          nall_offset = nlocal + nghost;
        }
      }
      #endif

      if (molecular) {
        for (int i = ifrom; i < ito; ++i) {
          int * _noalias jlist = firstneigh + cnumneigh[i];
          const int jnum = numneigh[i];
          #if defined(LMP_SIMD_COMPILER)
	  #pragma vector aligned
          #pragma simd
	  #endif
          for (int jj = 0; jj < jnum; jj++) {
            const int j = jlist[jj];
	    if (need_ic && j < 0) {
	      which = 0;
	      jlist[jj] = -j - 1;
	    } else
              ofind_special(which, special, nspecial, i, tag[j]);
            #ifdef _LMP_INTEL_OFFLOAD
	    if (j >= nlocal) {
	      if (j == e_nall)
		jlist[jj] = nall_offset;
	      else if (which) 
		jlist[jj] = (j-ghost_offset) ^ (which << SBBITS);
	      else jlist[jj]-=ghost_offset;
            } else
            #endif
	      if (which) jlist[jj] = j ^ (which << SBBITS);
          }
        }
      }
      #ifdef _LMP_INTEL_OFFLOAD
      else if (separate_buffers) {
	for (int i = ifrom; i < ito; ++i) {
          int * _noalias jlist = firstneigh + cnumneigh[i];
          const int jnum = numneigh[i];
	  int jj = 0;
	  for (jj = 0; jj < jnum; jj++)
	    if (jlist[jj] >= nlocal) break;
	  while (jj < jnum) {
	    if (jlist[jj] == e_nall) jlist[jj] = nall_offset;
	    else jlist[jj] -= ghost_offset;
	    jj++;
	  }
	}
      }
      #endif
    } // end omp
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime() - *timer_compute;
    #endif
  } // end offload

  #ifdef _LMP_INTEL_OFFLOAD
  if (offload) {
    _fix->stop_watch(TIME_OFFLOAD_LATENCY);
    _fix->start_watch(TIME_HOST_NEIGHBOR);
    for (int n = 0; n < aend; n++) {
      ilist[n] = n;
      numneigh[n] = 0;
    }
  } else {
    for (int i = astart; i < aend; i++)
      list->firstneigh[i] = firstneigh + cnumneigh[i];
    if (separate_buffers) {
      _fix->start_watch(TIME_PACK);
      _fix->set_neighbor_host_sizes();
      buffers->pack_sep_from_single(_fix->host_min_local(),
				    _fix->host_used_local(),
				    _fix->host_min_ghost(),
				    _fix->host_used_ghost());
      _fix->stop_watch(TIME_PACK);
    }
  }
  #else
  for (int i = astart; i < aend; i++)
    list->firstneigh[i] = firstneigh + cnumneigh[i];
  #endif
}
