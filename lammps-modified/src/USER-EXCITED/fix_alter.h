#ifdef FIX_CLASS

FixStyle(alter,FixAlter)

#else

#ifndef LMP_FIX_ALTER_H
#define LMP_FIX_ALTER_H

#include <set>
#include "fix.h"



namespace LAMMPS_NS {

class FixAlter : public Fix {
public:
  FixAlter(class LAMMPS *lammps, int narg, char **args);
  ~FixAlter();

  int modify_param(int narg, char **arg) override;

  virtual int setmask() override;
  void setup(int i) override;
  void min_setup(int i) override;
  void build_mol_list(int narg, char **arg);
  std::set<tagint> excite_mols;
};

}

#endif
#endif
