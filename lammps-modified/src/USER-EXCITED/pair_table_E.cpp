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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "pair_table_E.h"
#include "string.h"
#include "memory.h"

using namespace LAMMPS_NS;

enum{NONE,RLINEAR,RSQ,BMP};

#define MAXLINE 1024
#define EPSILONR 1.0e-6

/* ---------------------------------------------------------------------- */

PairTableE::PairTableE(class LAMMPS * lmp) : PairTable(lmp) {
  auto style = atom->style_match("fullExc");
  if(style == nullptr){
    error->one(FLERR, "PairTableE can be used only with fullExc atom style");
  }
}

PairTableE::~PairTableE() {
  fprintf(screen, "TableE dtor\n");
}
void PairTableE::coeff(int narg, char **arg) {
  if( narg != 7 ) error->all(FLERR, "Illegal argument count to PairTableE");
  char * arg1[6]{
           arg[0],
           arg[1],
           arg[2],
           arg[3],
           arg[5],
           arg[6]
          };
  PairTable::coeff(6, arg1);
  arg1[3] = arg[4];
  PairTable::coeff(6, arg1);
}

PairTable::Table *PairTableE::get_table(int i, int j, int itype, int jtype) {
  if(i == -1 || j == -1) {
    //Special case for Pair::init()
    return PairTable::get_table(i, j, itype, jtype);
  }

  if( setflag[itype][jtype] != 2 && setflag[jtype][itype] != 2 ) {
    char err[255];
    sprintf(err, "Pair type (%i:%i) doesn't have an excited state table", itype, jtype);
    error->one(FLERR, err);
  }

  return PairTable::get_table(i, j, itype, jtype) + (atom->estate[i] | atom->estate[j]);
}