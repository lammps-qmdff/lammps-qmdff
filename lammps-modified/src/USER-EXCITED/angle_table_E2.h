#ifdef ANGLE_CLASS

AngleStyle(tableE2,AngleTableE2)

#else

#ifndef LMP_ANGLE_TABLE_E2_H
#define LMP_ANGLE_TABLE_E2_H

#include "angle_table_E.h"
#include "bond_table_E.h"

namespace LAMMPS_NS {

class BondTableExternals;

class AngleTableE2 : public AngleTableE {
private:
  struct DualTable {
    int ninput;
    AngleTableE::Table ang_table;
    BondTableE::Table bond_table;
  };

  DualTable * d_tables = nullptr;
  void null_table(DualTable* table);
  const double MAXDIST = 10.5; //Maximum distance to compute pair force for

  BondTableExternals* bond;

public:
  void settings(int i, char **pString) override;
  double equilibrium_angle(int i) override;
  void write_restart(FILE *file) override;
  void read_restart(FILE *file) override;
  double single(int type, int i1, int i2, int i3) override;
  DualTable* allocate_table();

  explicit AngleTableE2(LAMMPS * lmp);
  ~AngleTableE2() override;

  void coeff(int narg, char **arg) override;

  virtual void read_dual_table(DualTable* tb, char* file, char* keyword);
  enum class TableType : int {LINEAR,SPLINE};

  std::map<std::string, long> fill_seek_cache(FILE *fp) override;
  void param_extract(DualTable *table, char *string);
  void bcast_table(DualTable *tb);

  void angle_table_fixup(AngleTableE::Table* tb);

  void coeff_one(int narg, char **arg);
  Table *get_table(int type, int i, int j, int k) override;
  void free_table(DualTable *table);

  BondTableE::Table* get_dist_table(int i, int j, int type);
  void compute(int i, int i1) override;
protected:
  void allocate() override;

  void bond_ev_tally(int i, int j, int nlocal, int newton_bond,
                                   double ebond, double fbond);

  void bond_ev_setup(int eflag, int vflag, int alloc = 1);

  friend class BondTableExternals;
};

}
#endif

#endif