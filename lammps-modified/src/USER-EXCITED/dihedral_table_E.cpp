#include <bitset>
#include "dihedral_table_E.h"
#include "error.h"
#include "atom.h"

using namespace LAMMPS_NS;

DihedralTableE::DihedralTableE(LAMMPS * lmp) : DihedralTable(lmp) {
  auto style = atom->style_match("fullExc");
  if(style == nullptr){
    error->one(FLERR, "PairTableE can be used only with fullExc atom style");
  }
}

void DihedralTableE::coeff(int narg, char **arg) {
  if( narg != 4 ) error->all(FLERR, "Illegal argument count to BondTableE");
  char *arg1[3]{
    arg[0],
    arg[1],
    arg[2],
  };
  DihedralTable::coeff(3, arg1);
  arg1[2] = arg[3];
  DihedralTable::coeff(3, arg1);
}

DihedralTable::Table *DihedralTableE::get_table(int type, int i1, int i2, int i3, int i4) {
  if(i1 == -1 || i2 == -1 || i3 == -1 || i4 == -1) {
    return DihedralTable::get_table(type, i1, i2, i3, i4);
  }

  if( setflag[type] != 2 ) {
    char err[255];
    sprintf(err, "Dihedral type %i doesn't have an excited state table", type);
    error->one(FLERR, err);
  }

  if( atom->estate[i1] ^ atom->estate[i2] ^ atom->estate[i3] ^ atom->estate[i4] ) {
    error->all(FLERR, "Sanity check Failed! Requested dihedral between excited and non-excited atoms ");
  }

  return DihedralTable::get_table(type, i1, i2, i3, i4) + ( atom->estate[i1] | atom->estate[i2] | atom->estate[i3] | atom->estate[i4] );
}

