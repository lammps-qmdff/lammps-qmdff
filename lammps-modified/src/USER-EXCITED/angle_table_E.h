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

#ifdef ANGLE_CLASS

//AngleStyle(tableE,AngleTableE)

#else

#ifndef LMP_ANGLE_TABLE_E_H
#define LMP_ANGLE_TABLE_E_H

#include "angle_table.h"

namespace LAMMPS_NS {

class AngleTableE : public AngleTable {
public:
  AngleTableE(LAMMPS *);
  virtual void coeff(int narg, char **arg) override;
  using Table = AngleTable::Table;
protected:
  Table *get_table(int type, int i, int j, int k) override;

};
}

#endif

#endif
