#include <cstring>
#include <vector>
#include "fix.h"
#include "fix_alter.h"
#include "force.h"
#include "atom.h"
#include "pair_hybrid.h"
#include "style_pair.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

int LAMMPS_NS::FixAlter::setmask() {
  int mask = 0;
  return mask;
}

FixAlter::FixAlter(class LAMMPS *lammps, int narg, char **args) :
  Fix(lammps, narg, args)
{
  auto style = atom->style_match("fullExc");
  if(style == nullptr){
    error->one(FLERR, "PairTableE can be used only with fullExc atom style");
  }

  //process atom list
  build_mol_list(narg - 3, args + 3);
}

void FixAlter::build_mol_list(int narg, char **arg)
{
  fprintf(screen, "FixAlter got %i args\n", narg);
  for( int i = 0; i < narg; ++i ) {
    fprintf(screen, "%s ", arg[i]);
  }
  fprintf(screen, "\n");

  for(int iarg = 0; iarg < narg; iarg++) {
    int atom_tag = force->inumeric(FLERR, arg[iarg]);
    int iatom = atom->map(atom_tag);
    if( iatom != -1 ) {
      auto mol = atom->molecule[iatom];
      fprintf(screen, "Switch mol %i to T1\n", mol);
      if( mol == 0 ) {
        error->one(FLERR, "Atom is not a part of a molecule");
      }
      excite_mols.insert(mol);
    } else {
      char err[256];
      sprintf(err, "Unknown atom specified");
      error->one(FLERR, err);
    }
  }
}

void FixAlter::setup(int vflag)
{
  for (int i = 0; i < atom->nlocal + atom->nghost; i++) {
    auto mol = atom->molecule[i];
    atom->q[i] = atom->q_grnd[i];
    atom->estate[i] = 0;
    if( excite_mols.find(mol) != excite_mols.end() ) {
        fprintf(screen, "Switch atom %i to T1\n", i);
        atom->estate[i] = 1;
        atom->q[i] = atom->q_exc[i];
    }
  }
}

int FixAlter::modify_param(int narg, char **arg)
{
  build_mol_list(narg, arg);
  return narg;
}

FixAlter::~FixAlter() {
  for (int i = 0; i < atom->nlocal + atom->nghost; i++) {
    auto mol = atom->molecule[i];
    if(excite_mols.find(mol) != excite_mols.end()) {
      atom->q[i] = atom->q_grnd[i];
      atom->estate[i] = 0;
    }
  }
}

void FixAlter::min_setup(int vflag) {
  this->setup(vflag);
}
