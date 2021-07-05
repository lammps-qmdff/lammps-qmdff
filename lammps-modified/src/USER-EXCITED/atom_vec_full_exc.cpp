#include "atom_vec_full_exc.h"

#include <stdlib.h>
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

AtomVecFullExc::AtomVecFullExc(LAMMPS * lmp ) : AtomVecFull(lmp) {
  xcol_data++;
  size_data_atom++;
  size_border += 3;
}

void AtomVecFullExc::data_atom(double *coord, imageint imagetmp, char **values) {
  int nlocal = atom->nlocal;
  AtomVecFull::data_atom(coord, imagetmp, values);
  q_exc[nlocal] = atof(values[4]);
  q_grnd[nlocal] = q[nlocal];
  estate[nlocal] = 0;
}

void AtomVecFullExc::grow(int n) {
  AtomVecFull::grow(n);
  //AtomVecFullExc changed
  q_grnd = memory->grow(atom->q_grnd, nmax, "atom:q_grnd");
  q_exc = memory->grow(atom->q_exc, nmax, "atom:q_exc");
  estate = memory->grow(atom->estate, nmax, "atom:estate");
}
void AtomVecFullExc::grow_reset() {
  AtomVecFull::grow_reset();
  q_exc = atom->q_exc;
  q_grnd = atom->q_grnd;
  estate = atom->estate;
}

void AtomVecFullExc::copy(int i, int j, int delflag) {
  q_exc[j] = q_exc[i];
  q_grnd[j] = q_grnd[i];
  estate[j] = estate[i];
  AtomVecFull::copy(i, j, delflag);
}

void AtomVecFullExc::unpack_border(int n, int first, double *buf) {
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax)
      grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    q[i] = buf[m++];
    q_grnd[i] = buf[m++];
    q_exc[i] = buf[m++];
    estate[i] = (int) ubuf(buf[m++]).i;
    molecule[i] = (tagint) ubuf(buf[m++]).i;
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n, first, &buf[m]);
}

void AtomVecFullExc::create_atom(int itype, double *coord) {
  int nlocal = atom->nlocal;
  AtomVecFull::create_atom(itype, coord);
  estate[nlocal] = 0;
  q_exc[nlocal] = 0;
  q_grnd[nlocal] = 0;
}

int AtomVecFullExc::pack_border(int n, int *list, double *buf,
                                int pbc_flag, int *pbc) {
  int i, j, m;
  double dx, dy, dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = q[j];
      buf[m++] = q_grnd[j];
      buf[m++] = q_exc[j];
      buf[m++] = ubuf(estate[j]).d;
      buf[m++] = ubuf(molecule[j]).d;
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = q[j];
      buf[m++] = q_grnd[j];
      buf[m++] = q_exc[j];
      buf[m++] = ubuf(estate[j]).d;
      buf[m++] = ubuf(molecule[j]).d;
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n, list, &buf[m]);

  return m;
}

int AtomVecFullExc::pack_border_vel(int i, int *pInt, double *pDouble, int i1, int *pInt1) {
  error->all(FLERR, __PRETTY_FUNCTION__);
  return AtomVecFull::pack_border_vel(i, pInt, pDouble, i1, pInt1);
}

int AtomVecFullExc::pack_border_hybrid(int i, int *pInt, double *pDouble) {
  error->all(FLERR, __PRETTY_FUNCTION__);
  return AtomVecFull::pack_border_hybrid(i, pInt, pDouble);
}

void AtomVecFullExc::unpack_border_vel(int i, int i1, double *pDouble) {
  error->all(FLERR, __PRETTY_FUNCTION__);
  AtomVecFull::unpack_border_vel(i, i1, pDouble);
}

int AtomVecFullExc::unpack_border_hybrid(int i, int i1, double *pDouble) {
  error->all(FLERR, __PRETTY_FUNCTION__);
  return AtomVecFull::unpack_border_hybrid(i, i1, pDouble);
}

int AtomVecFullExc::pack_exchange(int i, double *buf) {
  auto m = AtomVecFull::pack_exchange(i, buf);
  buf[m++] = ubuf(estate[i]).d;
  buf[m++] = q_exc[i];
  buf[m++] = q_grnd[i];
  buf[0] = m;
  return m;
}

int AtomVecFullExc::unpack_exchange(double *buf) {
  int nlocal = atom->nlocal;
  auto m = AtomVecFull::unpack_exchange(buf);
  estate[nlocal] = (int) ubuf(buf[m++]).i;
  q_exc[nlocal] = buf[m++];
  q_grnd[nlocal] = buf[m++];
  return m;
}

int AtomVecFullExc::size_restart() {
  int nlocal = atom->nlocal;
  auto ret = AtomVecFull::size_restart();
  ret += 3*nlocal;
  return ret;
}

int AtomVecFullExc::pack_restart(int i, double *buf) {
  int m = AtomVecFull::pack_restart(i, buf);
  buf[m++] = ubuf(estate[i]).d;
  buf[m++] = q_exc[i];
  buf[m++] = q_grnd[i];
  buf[0] = m;
  return m;
}

int AtomVecFullExc::unpack_restart(double *buf) {
  int nlocal = atom->nlocal;
  int m = AtomVecFull::unpack_restart(buf);

  estate[nlocal] = (int) ubuf(buf[m++]).i;
  q_exc[nlocal] = buf[m++];
  q_grnd[nlocal] = buf[m++];
  return m;
}

int AtomVecFullExc::data_atom_hybrid(int i, char **buf) {
  error->all(FLERR, __PRETTY_FUNCTION__);
  return AtomVecFull::data_atom_hybrid(i, buf);
}

int AtomVecFullExc::pack_data_one(int i, double **buf) {
  //Save current atom charge and restore it to ground state value to preserve consistency
  double qtmp = q[i];
  q[i] = q_grnd[i];
  int m = AtomVecFull::pack_data_one(i, buf);
  q[i] = qtmp;

  buf[i][m++] = q_exc[i];
  return m;
}

int AtomVecFullExc::pack_data_hybrid(int i, double *pDouble) {
  error->all(FLERR, __PRETTY_FUNCTION__);
  return AtomVecFull::pack_data_hybrid(i, pDouble);
}

void AtomVecFullExc::write_data(FILE *fp, int n, double **buf) {
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT "\t" TAGINT_FORMAT
              "\t%d\t%-1.6f\t%-1.6f\t%-1.3f\t%-1.3f\t%-1.3f\t%d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(tagint) ubuf(buf[i][1]).i,
            (int) ubuf(buf[i][2]).i,
            buf[i][3], buf[i][11],buf[i][4],buf[i][5],buf[i][6],
            (int) ubuf(buf[i][7]).i,(int) ubuf(buf[i][8]).i,
            (int) ubuf(buf[i][9]).i);
}

int AtomVecFullExc::write_data_hybrid(FILE *file, double *pDouble) {
  error->all(FLERR, __PRETTY_FUNCTION__);
  return AtomVecFull::write_data_hybrid(file, pDouble);
}


