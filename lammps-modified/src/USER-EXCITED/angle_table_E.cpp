#include "pointers.h"
#include "error.h"
#include "angle_table_E.h"
#include "atom.h"

using namespace LAMMPS_NS;

AngleTableE::AngleTableE(LAMMPS_NS::LAMMPS * lmp) : AngleTable(lmp){auto style = atom->style_match("fullExc");
  if(style == nullptr){
    error->one(FLERR, "PairTableE can be used only with fullExc atom style");
  }

}

void AngleTableE::coeff(int narg, char **arg) {
  if( narg != 4 ) error->all(FLERR, "Illegal argument count to AngleTableE");
  char *arg1[3]{
    arg[0],
    arg[1],
    arg[2],
  };
  AngleTable::coeff(3, arg1);
  arg1[2] = arg[3];
  AngleTable::coeff(3, arg1);
}
AngleTable::Table *AngleTableE::get_table(int type, int i, int j, int k) {
  if(i == -1 || j == -1 || k == -1) {
    return AngleTable::get_table(type, i, j, k);
  }

  if( setflag[type] != 2 ) {
    char err[255];
    sprintf(err, "Angle type %i doesn't have an excited state table", type);
    error->one(FLERR, err);
  }

  if( atom->estate[i] != atom->estate[j] ||  atom->estate[j] != atom->estate[k] ) {
    error->all(FLERR, "Sanity check Failed! Requested angle between excited and non-excited atoms ");
  }

  return AngleTable::get_table(type, i, j, k) + (atom->estate[i]);
}

