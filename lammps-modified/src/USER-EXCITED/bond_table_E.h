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

#ifdef BOND_CLASS

BondStyle(tableE,BondTableE)

#else

#ifndef LMP_BOND_TABLE_E_H
#define LMP_BOND_TABLE_E_H

#include <stdio.h>
#include "bond_table.h"

namespace LAMMPS_NS {

class BondTableE : public BondTable {
public:
  BondTableE(class LAMMPS *);
  virtual ~BondTableE();
  void coeff(int i, char **pString) override;
  using Table = BondTable::Table;
protected:
  Table *get_table(int type, int i, int j) override;
};

}

#endif
#endif