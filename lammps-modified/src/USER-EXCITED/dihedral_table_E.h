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

/* ----------------------------------------------------------------------
   Contributing author: Andrew Jewett (jewett.aij at g mail)
------------------------------------------------------------------------- */

#ifdef DIHEDRAL_CLASS

DihedralStyle(tableE,DihedralTableE)

#else

#ifndef LMP_DIHEDRAL_TABLE_E_H
#define LMP_DIHEDRAL_TABLE_E_H

#include "dihedral_table.h"

namespace LAMMPS_NS {

class DihedralTableE : public DihedralTable {
public:
  DihedralTableE(LAMMPS * lmp);
  virtual ~DihedralTableE() = default;

  void coeff(int i, char **pString) override;
protected:
  Table *get_table(int type, int i1, int i2, int i3, int i4) override;
};

}

#endif
#endif