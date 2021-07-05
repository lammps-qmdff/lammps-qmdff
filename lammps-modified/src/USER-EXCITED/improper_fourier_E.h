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

#ifdef IMPROPER_CLASS

ImproperStyle(fourierE,ImproperFourierE)

#else

#ifndef LMP_IMPROPER_FOURIER_E_H
#define LMP_IMPROPER_FOURIER_E_H

#include "improper_fourier.h"

namespace LAMMPS_NS {

class ImproperFourierE : public ImproperFourier {
protected:
  void allocate() override;
public:
  ImproperFourierE(LAMMPS *);
  virtual void coeff(int narg, char **arg) override;
  double getK(int i1, int i2, int i3, int i4, int type) const override;
  double getC0(int i1, int i2, int i3, int i4, int type) const override;
  double getC1(int i1, int i2, int i3, int i4, int type) const override;
  double getC2(int i1, int i2, int i3, int i4, int type) const override;
  int exc_offset;
};
}

#endif
#endif