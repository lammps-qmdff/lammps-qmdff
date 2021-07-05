#include "improper_fourier_E.h"
#include "force.h"
#include "error.h"
#include "atom.h"
#include "memory.h"

using namespace LAMMPS_NS;

void ImproperFourierE::coeff(int narg, char **arg) {
  if ( narg != 9 && narg != 10 ) error->all(FLERR,"Incorrect args for improper coefficients");

  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nimpropertypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double C0_one = force->numeric(FLERR,arg[2]);
  double C1_one = force->numeric(FLERR,arg[3]);
  double C2_one = force->numeric(FLERR,arg[4]);
  double  k_oneE = force->numeric(FLERR,arg[5]);
  double C0_oneE = force->numeric(FLERR,arg[6]);
  double C1_oneE = force->numeric(FLERR,arg[7]);
  double C2_oneE = force->numeric(FLERR,arg[8]);
  int all_one = 1;
  if ( narg == 10 ){
    all_one = force->inumeric(FLERR,arg[9]);
  }

  // convert w0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    C0[i] = C0_one;
    C1[i] = C1_one;
    C2[i] = C2_one;
    all[i] = all_one;
    setflag[i] = 1;
    k[i+exc_offset] = k_oneE;
    C0[i+exc_offset] = C0_oneE;
    C1[i+exc_offset] = C1_oneE;
    C2[i+exc_offset] = C2_oneE;
    setflag[i+exc_offset] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for improper coefficients");
}

ImproperFourierE::ImproperFourierE(LAMMPS_NS::LAMMPS * lmp) : ImproperFourier(lmp) {
}

void ImproperFourierE::allocate() {
  allocated = 1;
  int n = atom->nimpropertypes*2;
  exc_offset = atom->nimpropertypes;

  memory->create(k,n+1,"improper:k");
  memory->create(C0,n+1,"improper:C0");
  memory->create(C1,n+1,"improper:C1");
  memory->create(C2,n+1,"improper:C2");
  memory->create(all,n+1,"improper:C2");

  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

double ImproperFourierE::getK(int i1, int i2, int i3, int i4, int type) const {
  int* estate = atom->estate;

  if(!(estate[i1]==estate[i2] && estate[i2] == estate[i3] && estate[i3] == estate[i4])) {
    error->all(FLERR, "Sanity check Failed! Requested Improper between excited and non-excited atoms ");
  }

  return k[(estate[i1] == 1 ? type+exc_offset : type)];
}

double ImproperFourierE::getC0(int i1, int i2, int i3, int i4, int type) const {
  int* estate = atom->estate;

  if(!(estate[i1]==estate[i2] && estate[i2] == estate[i3] && estate[i3] == estate[i4])) {
    error->all(FLERR, "Sanity check Failed! Requested Improper between excited and non-excited atoms ");
  }

  return C0[(estate[i1] == 1 ? type+exc_offset : type)];
}

double ImproperFourierE::getC1(int i1, int i2, int i3, int i4, int type) const {
  int* estate = atom->estate;

  if(!(estate[i1]==estate[i2] && estate[i2] == estate[i3] && estate[i3] == estate[i4])) {
    error->all(FLERR, "Sanity check Failed! Requested Improper between excited and non-excited atoms ");
  }

  return C1[(estate[i1] == 1 ? type+exc_offset : type)];;
}

double ImproperFourierE::getC2(int i1, int i2, int i3, int i4, int type) const {
  int* estate = atom->estate;

  if(!(estate[i1]==estate[i2] && estate[i2] == estate[i3] && estate[i3] == estate[i4])) {
    error->all(FLERR, "Sanity check Failed! Requested Improper between excited and non-excited atoms ");
  }

  return C2[(estate[i1] == 1 ? type+exc_offset : type)];
}
