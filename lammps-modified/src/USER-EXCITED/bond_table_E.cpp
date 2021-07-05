
#include "bond_table_E.h"
#include <bitset>
#include "pointers.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

BondTableE::BondTableE(LAMMPS * lmp) : BondTable(lmp) {
  auto style = atom->style_match("fullExc");
  if(style == nullptr){
    error->one(FLERR, "PairTableE can be used only with fullExc atom style");
  }
}

BondTableE::~BondTableE() = default;

void BondTableE::coeff(int narg, char **arg) {
  if( narg != 4 ) error->all(FLERR, "Illegal argument count to BondTableE");
  char *arg1[3]{
    arg[0],
    arg[1],
    arg[2],
  };
  BondTable::coeff(3, arg1);
  arg1[2] = arg[3];
  BondTable::coeff(3, arg1);
}

BondTable::Table *BondTableE::get_table(int type, int i, int j) {
  if(i == -1 || j == -1) {
    return BondTable::get_table(type, i, j);
  }

  if( setflag[type] != 2 ) {
    char err[255];
    sprintf(err, "Bond type %i doesn't have an excited state table", type);
    error->one(FLERR, err);
  }

  if( atom->estate[i] != atom->estate[j] ) {
    error->all(FLERR, "Sanity check Failed! Requested bond between excited and non-excited atoms ");
  }

  return BondTable::get_table(type, i, j) + (atom->estate[i]);
}

