/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(fullExc,AtomVecFullExc)

#else

#ifndef LMP_ATOM_VEC_FULL_EXC_H
#define LMP_ATOM_VEC_FULL_EXC_H

#include "atom_vec_full.h"

namespace LAMMPS_NS {

//TODO: fix MPI data exchange. I.e. (un)pack_exchange(), (un)pack_data

class AtomVecFullExc : public AtomVecFull {
public:
  AtomVecFullExc(LAMMPS * lmp);
  void data_atom(double *coord, imageint imagetmp, char **values) override;
  void grow(int n) override;
  void grow_reset() override;
  void copy(int i, int j, int delflag) override;
  void create_atom(int itype, double *coord) override;
  void unpack_border(int i, int i1, double *pDouble) override;
  int pack_border(int n, int *list, double *buf, int pbc_flag, int *pbc) override;
  int pack_border_vel(int i, int *pInt, double *pDouble, int i1, int *pInt1) override;
  int pack_border_hybrid(int i, int *pInt, double *pDouble) override;
  void unpack_border_vel(int i, int i1, double *pDouble) override;
  int unpack_border_hybrid(int i, int i1, double *pDouble) override;
  int pack_exchange(int i, double *buf) override;
  int unpack_exchange(double *buf) override;
  int size_restart() override;
  int pack_restart(int i, double *pDouble) override;
  int unpack_restart(double *pDouble) override;
  int data_atom_hybrid(int i, char **buf) override;
  int pack_data_one(int, double**) override;
  int pack_data_hybrid(int i, double *pDouble) override;
  void write_data(FILE *file, int i, double **pDouble) override;
  int write_data_hybrid(FILE *file, double *pDouble) override;

protected:
  int* estate;
  double* q_exc;
  double* q_grnd;
};

}

#endif
#endif