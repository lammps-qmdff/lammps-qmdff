
#include <cstring>
#include "angle_table_E2.h"
#include "math_const.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "atom.h"
#include <cmath>
#include <cerrno>
#include "neighbor.h"
#include "comm.h"
#include "bond_table_E.h"

using namespace LAMMPS_NS;
using namespace MathConst;

//TODO null_table
//TODO uf_loockup/u_loockup for both angle and pair
//TODO change get_table to generic pointer
//TODO Inherit AngleTableE2 from Angle?
//TODO see param_extract

//TODO see into settings() method. Theres something wrong with it causing BondTable to break

#define EPSILONR 1.0e-6
#define BIGNUM 1.0e300

namespace LAMMPS_NS {

class BondTableExternals : Pointers {
public:
  AngleTableE2 &owner;
  using Table = BondTableE::Table;

  explicit BondTableExternals(AngleTableE2 &own, LAMMPS *lmp) : owner(own), Pointers(lmp) {};
  double *r0;

  int r0idx = -1;
  double emin = BIGNUM;

  enum class TableType : int {
    NONE, LINEAR, SPLINE
  };

  TableType tabstyle;
  int tablength = -1;

  void spline(double *x, double *y, int n, double yp1, double ypn, double *y2) {
    int i,k;
    double p,qn,sig,un;
    double *u = new double[n];

    if (yp1 > 0.99e30) y2[0] = u[0] = 0.0;
    else {
      y2[0] = -0.5;
      u[0] = (3.0/(x[1]-x[0])) * ((y[1]-y[0]) / (x[1]-x[0]) - yp1);
    }
    for (i = 1; i < n-1; i++) {
      sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
      p = sig*y2[i-1] + 2.0;
      y2[i] = (sig-1.0) / p;
      u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
      u[i] = (6.0*u[i] / (x[i+1]-x[i-1]) - sig*u[i-1]) / p;
    }
    if (ypn > 0.99e30) qn = un = 0.0;
    else {
      qn = 0.5;
      un = (3.0/(x[n-1]-x[n-2])) * (ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]));
    }
    y2[n-1] = (un-qn*u[n-2]) / (qn*y2[n-2] + 1.0);
    for (k = n-2; k >= 0; k--) y2[k] = y2[k]*y2[k+1] + u[k];

    delete [] u;
  }

  double splint(double *xa, double *ya, double *y2a, int n, double x) {
    int klo, khi, k;
    double h, b, a, y;

    klo = 0;
    khi = n - 1;
    while (khi - klo > 1) {
      k = (khi + klo) >> 1;
      if (xa[k] > x)
        khi = k;
      else
        klo = k;
    }
    h = xa[khi] - xa[klo];
    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    y = a * ya[klo] + b * ya[khi] +
      ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
    return y;
  }

  void begin_table_read(Table *tb, char* name) {
    memory->create(tb->rfile, tb->ninput, "bond:rfile");
    memory->create(tb->efile, tb->ninput, "bond:efile");
    memory->create(tb->ffile, tb->ninput, "bond:ffile");

//    strcpy(tb->name, name);
  }

  void end_table_read(Table *tb) {
    // infer r0 from minimum of potential, if not given explicitly

    if ((tb->r0 == 0.0) && (r0idx >= 0))
      tb->r0 = tb->rfile[r0idx];

    // warn if force != dE/dr at any point that is not an inflection point
    // check via secant approximation to dE/dr
    // skip two end points since do not have surrounding secants
    // inflection point is where curvature changes sign

    double r, e, f, rprev, rnext, eprev, enext, fleft, fright;

    int ferror = 0;
    for (int i = 1; i < tb->ninput - 1; i++) {
      r = tb->rfile[i];
      rprev = tb->rfile[i - 1];
      rnext = tb->rfile[i + 1];
      e = tb->efile[i];
      eprev = tb->efile[i - 1];
      enext = tb->efile[i + 1];
      f = tb->ffile[i];
      fleft = -(e - eprev) / (r - rprev);
      fright = -(enext - e) / (rnext - r);
      if (f < fleft && f < fright)
        ferror++;
      if (f > fleft && f > fright)
        ferror++;
      //printf("Values %d: %g %g %g\n",i,r,e,f);
      //printf("  secant %d %d %g: %g %g %g\n",i,ferror,r,fleft,fright,f);
    }
  }

  void spline_table(Table *tb) {
    memory->create(tb->e2file, tb->ninput, "bond:e2file");
    memory->create(tb->f2file, tb->ninput, "bond:f2file");

    double ep0 = -tb->ffile[0];
    double epn = -tb->ffile[tb->ninput - 1];
    spline(tb->rfile, tb->efile, tb->ninput, ep0, epn, tb->e2file);

    if (tb->fpflag == 0) {
      tb->fplo = (tb->ffile[1] - tb->ffile[0]) / (tb->rfile[1] - tb->rfile[0]);
      tb->fphi = (tb->ffile[tb->ninput - 1] - tb->ffile[tb->ninput - 2]) /
        (tb->rfile[tb->ninput - 1] - tb->rfile[tb->ninput - 2]);
    }

    double fp0 = tb->fplo;
    double fpn = tb->fphi;
    spline(tb->rfile, tb->ffile, tb->ninput, fp0, fpn, tb->f2file);
  }

  void compute_table(Table *tb) {
    // delta = table spacing for N-1 bins
    int tlm1 = tablength - 1;

    tb->delta = (tb->hi - tb->lo) / tlm1;
    tb->invdelta = 1.0 / tb->delta;
    tb->deltasq6 = tb->delta * tb->delta / 6.0;

    // N-1 evenly spaced bins in r from min to max
    // r,e,f = value at lower edge of bin
    // de,df values = delta values of e,f
    // r,e,f are N in length so de,df arrays can compute difference

    memory->create(tb->r, tablength, "bond:r");
    memory->create(tb->e, tablength, "bond:e");
    memory->create(tb->de, tlm1, "bond:de");
    memory->create(tb->f, tablength, "bond:f");
    memory->create(tb->df, tlm1, "bond:df");
    memory->create(tb->e2, tablength, "bond:e2");
    memory->create(tb->f2, tablength, "bond:f2");

    double a;
    for (int i = 0; i < tablength; i++) {
      a = tb->lo + i * tb->delta;
      tb->r[i] = a;
      tb->e[i] = splint(tb->rfile, tb->efile, tb->e2file, tb->ninput, a);
      tb->f[i] = splint(tb->rfile, tb->ffile, tb->f2file, tb->ninput, a);
    }

    for (int i = 0; i < tlm1; i++) {
      tb->de[i] = tb->e[i + 1] - tb->e[i];
      tb->df[i] = tb->f[i + 1] - tb->f[i];
    }

    double ep0 = -tb->f[0];
    double epn = -tb->f[tlm1];
    spline(tb->r, tb->e, tablength, ep0, epn, tb->e2);
    spline(tb->r, tb->f, tablength, tb->fplo, tb->fphi, tb->f2);
  }

  void table_fixup(Table *tb) {
    tb->lo = tb->rfile[0];
    tb->hi = tb->rfile[tb->ninput - 1];
    if (tb->lo >= tb->hi)
      error->all(FLERR, "Bond table values are not increasing");
  }

  void table_bind(Table *tb) {

  }

  void bcast_table(Table *tb) {
    MPI_Bcast(&tb->ninput, 1, MPI_INT, 0, world);
    //MPI_Bcast(&tb->r0,1,MPI_INT,0,world);
    MPI_Bcast(&tb->r0, 1, MPI_DOUBLE, 0, world); // CHANGED from previous line by sasha

    int me;
    MPI_Comm_rank(world, &me);
    if (me > 0) {
      memory->create(tb->rfile, tb->ninput, "angle:rfile");
      memory->create(tb->efile, tb->ninput, "angle:efile");
      memory->create(tb->ffile, tb->ninput, "angle:ffile");
    }

    MPI_Bcast(tb->rfile, tb->ninput, MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->efile, tb->ninput, MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->ffile, tb->ninput, MPI_DOUBLE, 0, world);

    MPI_Bcast(&tb->fpflag, 1, MPI_INT, 0, world);
    if (tb->fpflag) {
      MPI_Bcast(&tb->fplo, 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&tb->fphi, 1, MPI_DOUBLE, 0, world);
    }
    //MPI_Bcast(&tb->r0,1,MPI_INT,0,world);
    MPI_Bcast(&tb->r0, 1, MPI_DOUBLE, 0, world); // CHANGED from previoud lin eby SASHA
  }

  void uf_lookup(int type, double x, double &u, double &f, int i, int j) {
    int itable;
    double fraction, a, b;
    char estr[128];

    Table *tb = get_table(type, i, j);
    if (x < tb->lo) {
      sprintf(estr, "Bond length < table inner cutoff: "
        "type %d length %g", type, x);
      error->one(FLERR, estr);
    }
    if (x > tb->hi) {
      sprintf(estr, "Bond length > table outer cutoff: "
        "type %d length %g", type, x);
      error->one(FLERR, estr);
    }

    if (tabstyle == TableType::LINEAR) {
      itable = static_cast<int> ((x - tb->lo) * tb->invdelta);
      fraction = (x - tb->r[itable]) * tb->invdelta;
      u = tb->e[itable] + fraction * tb->de[itable];
      f = tb->f[itable] + fraction * tb->df[itable];
    } else if (tabstyle == TableType::SPLINE) {
      itable = static_cast<int> ((x - tb->lo) * tb->invdelta);
      fraction = (x - tb->r[itable]) * tb->invdelta;

      b = (x - tb->r[itable]) * tb->invdelta;
      a = 1.0 - b;
      u = a * tb->e[itable] + b * tb->e[itable + 1] +
        ((a * a * a - a) * tb->e2[itable] + (b * b * b - b) * tb->e2[itable + 1]) *
          tb->deltasq6;
      f = a * tb->f[itable] + b * tb->f[itable + 1] +
        ((a * a * a - a) * tb->f2[itable] + (b * b * b - b) * tb->f2[itable + 1]) *
          tb->deltasq6;
    }
  }

  Table *get_table(int type, int i, int j) {
    return owner.get_dist_table(i, j, type);
  }

  void u_lookup(int type, double x, double &u) {
    int itable;
    double fraction, a, b;

    Table *tb = get_table(type, -1, -1);
    x = MAX(x, tb->lo);
    x = MIN(x, tb->hi);

    if (tabstyle == TableType::LINEAR) {
      itable = static_cast<int> ((x - tb->lo) * tb->invdelta);
      fraction = (x - tb->r[itable]) * tb->invdelta;
      u = tb->e[itable] + fraction * tb->de[itable];
    } else if (tabstyle == TableType::SPLINE) {
      itable = static_cast<int> ((x - tb->lo) * tb->invdelta);
      fraction = (x - tb->r[itable]) * tb->invdelta;

      b = (x - tb->r[itable]) * tb->invdelta;
      a = 1.0 - b;
      u = a * tb->e[itable] + b * tb->e[itable + 1] +
        ((a * a * a - a) * tb->e2[itable] + (b * b * b - b) * tb->e2[itable + 1]) *
          tb->deltasq6;
    }
  }

  void fill_line(Table *tb, char *tline, int i) {
    int e;
    double dtmp = strtod(tline, &tline);
    tb->rfile[i] = dtmp;
    e = errno;
    tb->efile[i] = strtod(tline, &tline);
    e = errno;
    tb->ffile[i] = strtod(tline, &tline);
    e = errno;

    if (tb->efile[i] < emin) {
      emin = tb->efile[i];
      r0idx = i;
    }
  }

  void compute(int i1, int i2, int type, double &ebond, double& fbond, int eflag, int vflag) {
    int nlocal = atom->nlocal;
    double **x = atom->x;
    double u, mdu;
    int newton_bond = force->newton_bond;

    double *f1n = atom->f[i1];
    double *f3n = atom->f[i2];

    double delx = x[i1][0] - x[i2][0];
    double dely = x[i1][1] - x[i2][1];
    double delz = x[i1][2] - x[i2][2];

    double rsq = delx * delx + dely * dely + delz * delz;
    double r = sqrt(rsq);

    // force & energy

    uf_lookup(type, r, u, mdu, i1, i2);
    fbond = mdu / r;
    ebond = u;


    if (newton_bond || i1 < nlocal) {
      f1n[0] += delx * fbond;
      f1n[1] += dely * fbond;
      f1n[2] += delz * fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f3n[0] -= delx * fbond;
      f3n[1] -= dely * fbond;
      f3n[2] -= delz * fbond;
    }
  }

  double single(int type, int i, int j, double &fforce) {
    double **x = atom->x;
    double delx = x[i][0] - x[j][0];
    double dely = x[i][1] - x[j][1];
    double delz = x[i][2] - x[j][2];

    double rsq = delx * delx + dely * dely + delz * delz;
    double r = sqrt(rsq);

    double u;
    double mdu;
    uf_lookup(type, r, u, mdu, -1, -1);
    fforce = mdu / r;
    return u;
  }

  void free_table(Table *tb) {
    memory->destroy(tb->rfile);
    memory->destroy(tb->efile);
    memory->destroy(tb->ffile);
    memory->destroy(tb->e2file);
    memory->destroy(tb->f2file);

    memory->destroy(tb->r);
    memory->destroy(tb->e);
    memory->destroy(tb->de);
    memory->destroy(tb->f);
    memory->destroy(tb->df);
    memory->destroy(tb->e2);
    memory->destroy(tb->f2);
  }

  void allocate() {
    int n = atom->nangletypes;
    memory->create(r0, n + 1, "bond:r0");
  }

  void settings(int tablength, TableType tabstyle) {
    this->tablength = tablength;
    this->tabstyle = tabstyle;
  }


  void dump_table(Table* tb, int idx) {
    int tlm1 = tablength - 1;
    char fname[25];
    sprintf(fname, "dump_tbl_E2_%i", idx );
    FILE* dump = fopen(fname, "w+");
    fprintf(dump, "\n\n");

    fprintf(dump, "ninput=%i, fpflag=%i, fplo=%f, fphi=%f, r0=%f\n", tb->ninput, tb->fpflag, tb->fplo, tb->fphi, tb->r0);
    fprintf(dump, "delta=%.5f, invdelta=%.5f, deltasq6=%.5f\n",tb->delta, tb->invdelta, tb->deltasq6 );

    for( int i = 0; i < tlm1; i++ ) {
      double rfile = tb->rfile? tb->rfile[i] : FP_NAN;
      double efile = tb->efile? tb->efile[i] : FP_NAN;
      double ffile = tb->ffile? tb->ffile[i] : FP_NAN;
      double e2file = tb->e2file? tb->e2file[i] : FP_NAN;
      double f2file = tb->f2file? tb->f2file[i] : FP_NAN;
      double r = tb->r? tb->r[i] : FP_NAN;
      double e = tb->e? tb->e[i] : FP_NAN;
      double de = tb->de? tb->de[i] : FP_NAN;
      double f = tb->f? tb->f[i] : FP_NAN;
      double df = tb->df? tb->df[i] : FP_NAN;
      double e2 = tb->e2? tb->e2[i] : FP_NAN;
      double f2 = tb->f2? tb->f2[i] : FP_NAN;

      fprintf(dump, "rfile=%.5f, efile=%.5f, ffile=%.5f, e2file=%.5f, f2file=%.5f, r=%.5f, e=%.5f, de=%.5f, f=%.5f, df=%.5f, e2=%.5f, f2=%.5f\n",
              rfile,
              efile,
              ffile,
              e2file,
              f2file,
              r,
              e,
              de,
              f,
              df,
              e2,
              f2
      );
    }

    fclose(dump);
  }

};

}



AngleTableE2::AngleTableE2(LAMMPS *lmp) : AngleTableE(lmp), bond(new BondTableExternals(*this, lmp))
{
  fprintf(screen, "AngleTableE2\n");
}

const int MAXLINE = 1024;

std::map<std::string, long> AngleTableE2::fill_seek_cache(FILE* fp) {
  std::map<std::string, long> seekpos;
  char line[MAXLINE];
  DualTable* tb = allocate_table();

  while (!feof(fp)) {
    fgets(line, MAXLINE, fp);
    if (strspn(line," \t\n\r") == strlen(line)) continue;  // blank line
    if (line[0] == '#') continue;                          // comment
    char *word = strtok(line," \t\n\r");
    long pos = ftell(fp);
    std::string key(word);
    seekpos[key] = pos;

    fgets(line,MAXLINE,fp);                         // no match, skip section
    param_extract(tb,line);
    fgets(line,MAXLINE,fp);
    for (int i = 0; i < tb->ninput+1; i++) fgets(line, MAXLINE, fp);
  }
  memory->sfree(tb);
  return seekpos;
}


static const double SMALL = 0.001;
static const double TINY = 1.E-10;


void AngleTableE2::angle_table_fixup(AngleTableE::Table* tb) {
  if (tb->ninput <= 1) error->one(FLERR,"Invalid angle table length");

  double alo,ahi;
  alo = tb->afile[0];
  ahi = tb->afile[tb->ninput-1];
  if (fabs(alo-0.0) > TINY || fabs(ahi-180.0) > TINY)
    error->all(FLERR,"Angle table must range from 0 to 180 degrees");

  // convert theta from degrees to radians

  for (int i = 0; i < tb->ninput; i++){
    tb->afile[i] *= MY_PI/180.0;
    tb->ffile[i] *= 180.0/MY_PI;
  }
}

void AngleTableE2::coeff(int narg, char **arg) {
  if (narg != 4 && narg != 3)
    error->all(FLERR, "Illegal argument count to AngleTableE");
  char *arg1[3]{
    arg[0],
    arg[1],
    arg[2],
  };
  coeff_one(3, arg1);
  if ( narg == 4 ) {
    arg1[2] = arg[3];
    coeff_one(3, arg1);
  }
}

void AngleTableE2::coeff_one(int narg, char **arg) {
  if (narg != 3) error->all(FLERR,"Illegal angle_coeff command");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nangletypes,ilo,ihi);

  int me;
  MPI_Comm_rank(world,&me);
  d_tables = (DualTable *)
    memory->srealloc(d_tables,(ntables+1)*sizeof(DualTable),"angle:d_tables");

  DualTable *tb = &d_tables[ntables];

  null_table(tb);

  if ( me == 0 ) read_dual_table(tb, arg[1], arg[2]);
  bcast_table(tb);

  // error check on table parameters
  angle_table_fixup(&tb->ang_table);
  bond->table_fixup(&tb->bond_table);
  // spline read-in and compute a,e,f vectors within table

  tb->ninput = tb->ang_table.ninput;



  spline_table(&tb->ang_table);
  compute_table(&tb->ang_table);

  if(ntables == 0 ) {
    bond->dump_table(&tb->bond_table, 1);
  }

  bond->spline_table(&tb->bond_table);

  if(ntables == 0 ) {
    bond->dump_table(&tb->bond_table, 2);
  }

  bond->compute_table(&tb->bond_table);

  if(ntables == 0 ) {
    bond->dump_table(&tb->bond_table, 3);
  }

  // store ptr to table in tabindex

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    if( setflag[i] == 0 ) tabindex[i] = ntables;
    setflag[i]++;
    theta0[i] = tb->ang_table.theta0;
    bond->table_bind(&tb->bond_table);
    count++;
  }
  ntables++;



  if (count == 0) error->all(FLERR,"Illegal angle_coeff command");
}

void AngleTableE2::bcast_table(DualTable *tb) {
  AngleTable::Table* ang_tb = &tb->ang_table;
  BondTableE::Table* bond_tb = &tb->bond_table;

  AngleTable::bcast_table(ang_tb);
  bond->bcast_table(bond_tb);
}

void AngleTableE2::param_extract(DualTable *dtb, char *line)
{
  //TODO additional params for force/energy table
//  fprintf(stderr, "%s(%p, \"%s\")\n", __func__, dtb, line);
  Table* tb = &dtb->ang_table;

  tb->ninput = 0;
  tb->fpflag = 0;
  tb->theta0 = 180.0;

  char *word = strtok(line," \t\n\r\f");
  while (word) {
    if (strcmp(word,"N") == 0) {
      word = strtok(NULL," \t\n\r\f");
      tb->ninput = atoi(word);
    } else if (strcmp(word,"FP") == 0) {
      tb->fpflag = 1;
      word = strtok(NULL," \t\n\r\f");
      tb->fplo = atof(word);
      word = strtok(NULL," \t\n\r\f");
      tb->fphi = atof(word);
      tb->fplo *= (180.0/MY_PI)*(180.0/MY_PI);
      tb->fphi *= (180.0/MY_PI)*(180.0/MY_PI);
    } else if (strcmp(word,"EQ") == 0) {
      word = strtok(NULL," \t\n\r\f");
      tb->theta0 = atof(word);
    } else {
      error->one(FLERR,"Invalid keyword in angle table parameters");
    }
    word = strtok(NULL," \t\n\r\f");
  }

  if (tb->ninput == 0) error->one(FLERR,"Angle table parameters did not set N");
  dtb->ninput = dtb->bond_table.ninput = dtb->ang_table.ninput;
}

void AngleTableE2::read_dual_table(DualTable *tb, char *file, char *keyword) {
  char line[MAXLINE];

  // open file

  FILE *fp = force->open_potential(file);
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(FLERR,str);
  }

  std::string filename(file);
  auto it = seek_pos_cache.find(filename);
  if( it == seek_pos_cache.end() ) {
    // reading file for first time, fill keyword_pos
    it = seek_pos_cache.insert(std::make_pair(filename, fill_seek_cache(fp))).first;
  }

  std::map<std::string, long> &keyword_pos = it->second;

  std::string current_key(keyword);

  auto pos_it = keyword_pos.find(current_key);
  if( pos_it != keyword_pos.end() ) {
    fseek(fp, pos_it->second, SEEK_SET);
  } else
  {
    char err[255];
    sprintf(err, "Did not find keyword \"%s\" in table file", keyword);
    error->one(FLERR, err);
    abort();
  }

  // read args on 2nd line of section
  // allocate table arrays for file values

  BondTableE::Table * bond_table = &tb->bond_table;
  AngleTableE::Table * angle_table = &tb->ang_table;


  fgets_unlocked(line,MAXLINE,fp);
  param_extract(tb,line);
  bond->begin_table_read(bond_table, keyword);
  memory->create(angle_table->afile,angle_table->ninput,"angle:afile");
  memory->create(angle_table->efile,angle_table->ninput,"angle:efile");
  memory->create(angle_table->ffile,angle_table->ninput,"angle:ffile");

  double rfile,rnew;

  int rerror = 0;
  int cerror = 0;

  fgets(line,MAXLINE,fp);

  for (int i = 0; i < tb->ninput; i++) {
    if (NULL == fgets(line, MAXLINE, fp))
      error->one(FLERR, "Premature end of file in angle table");

//    if (4 != sscanf(line,"%d %lg %lg %lg",
//                    &itmp,&rfile,&tb->efile[i],&tb->ffile[i]))  ++cerror;

    {
      char *tline = line;
      strtol(tline, &tline, 10);
      angle_table->afile[i] = strtod(tline, &tline);
      angle_table->efile[i] = strtod(tline, &tline);
      angle_table->ffile[i] = strtod(tline, &tline);
      bond->fill_line(bond_table, tline, i);
    }
  }

  bond->end_table_read(bond_table);

  fclose(fp);

}

AngleTableE2::DualTable *AngleTableE2::allocate_table() {
  DualTable * tb = (DualTable*)memory->smalloc(sizeof(DualTable), "temp_pair");
  tb->ninput = tb->ang_table.ninput;
  return tb;
}


void AngleTableE2::null_table(DualTable * tbl) {
  tbl->ninput = 0;
  memset(&tbl->ang_table, 0, sizeof(tbl->ang_table));
  memset(&tbl->bond_table, 0, sizeof(tbl->bond_table));
}

void AngleTableE2::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double eangle,f1[3],f3[3];
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;
  double theta,u,mdu; //mdu: minus du, -du/dx=f

  double ebond;
  double fbond;
  double bond_f1[3], bond_f2[3];

  eangle = 0.0;
  if (eflag || vflag){
    ev_setup(eflag,vflag);
//Bond::ev_setup() and Angle::ev_setup() do the same thing so I assume it's safe to call only one
//    bond_ev_setup(eflag, vflag);
  }
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // angle (cos and sin)

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;

    // tabulated force & energy

    theta = acos(c);
    uf_lookup(type,theta,u,mdu, i1, i2, i3);

    if (eflag) eangle = u;

    a = mdu * s;
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;
    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    //add force and energy values from pair interaction

    bond->compute(i1, i3, type, ebond, fbond, eflag, vflag);


    // apply force to each of 3 atoms
    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);

    if (evflag) bond_ev_tally(i1, i3, nlocal, newton_bond, ebond, fbond);

  }
}

AngleTableE2::~AngleTableE2() {
  for (int m = 0; m < ntables; m++) free_table(&d_tables[m]);
  memory->destroy(d_tables);
  delete bond;
  ntables = 0;
}

void AngleTableE2::free_table(DualTable *table) {
  AngleTableE::Table* ang_table = &table->ang_table;
  memory->destroy(ang_table->afile);
  memory->destroy(ang_table->efile);
  memory->destroy(ang_table->ffile);
  memory->destroy(ang_table->e2file);
  memory->destroy(ang_table->f2file);

  memory->destroy(ang_table->ang);
  memory->destroy(ang_table->e);
  memory->destroy(ang_table->de);
  memory->destroy(ang_table->f);
  memory->destroy(ang_table->df);
  memory->destroy(ang_table->e2);
  memory->destroy(ang_table->f2);

  BondTableE::Table* bond_table = &table->bond_table;
  bond->free_table(bond_table);
}

void AngleTableE2::settings(int narg, char **arg) {
//  AngleTable::settings(narg,arg);
//  return;

  if (narg != 2) error->all(FLERR,"Illegal angle_style command");

  if (strcmp(arg[0],"linear") == 0) tabstyle = int(TableType::LINEAR);
  else if (strcmp(arg[0],"spline") == 0) tabstyle = int(TableType::SPLINE);
  else error->all(FLERR,"Unknown table style in angle style table");

  tablength = force->inumeric(FLERR,arg[1]);
  if (tablength < 2) error->all(FLERR,"Illegal number of angle table entries");


  BondTableExternals::TableType bond_tabstyle;

  switch(TableType(tabstyle)) {
    case TableType::LINEAR:
      bond_tabstyle = BondTableExternals::TableType::LINEAR;
      break;
    case TableType::SPLINE:
      bond_tabstyle = BondTableExternals::TableType::SPLINE;
      break;
    default:
      bond_tabstyle = BondTableExternals::TableType::NONE;
      break;
  }

  bond->settings(tablength, bond_tabstyle);

  // delete old tables, since cannot just change settings

  for (int m = 0; m < ntables; m++) AngleTable::free_table(&d_tables[m].ang_table);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(tabindex);
  }
  allocated = 0;

  ntables = 0;
  d_tables = NULL;
  tables = NULL;
}

double AngleTableE2::equilibrium_angle(int i) {
  error->one(FLERR, (std::string(__PRETTY_FUNCTION__) + " unimplemented").c_str());
  return AngleTable::equilibrium_angle(i);
}

void AngleTableE2::write_restart(FILE *file) {
  error->one(FLERR, (std::string(__PRETTY_FUNCTION__) + " unimplemented").c_str());
  AngleTable::write_restart(file);
}

void AngleTableE2::read_restart(FILE *file) {
  error->one(FLERR, (std::string(__PRETTY_FUNCTION__) + " unimplemented").c_str());
  AngleTable::read_restart(file);
}


double AngleTableE2::single(int type, int i1, int i2, int i3) {
  double e = AngleTable::single(type, i1, i2, i3);
  double bond_fforce;
  e += bond->single(type, i1, i3, bond_fforce);
  return e;
}


AngleTableE::Table* AngleTableE2::get_table(int type, int i, int j, int k) {
  if(i == -1 || j == -1 || k == -1 || atom->estate == NULL) {
    return &(d_tables + tabindex[type])->ang_table;
  }

  if( atom->estate[i] != atom->estate[j] ||  atom->estate[j] != atom->estate[k] ) {
    error->all(FLERR, "Sanity check Failed! Requested angle between excited and non-excited atoms ");
  }

  if( (atom->estate[i] != 0 ) && (setflag[type] != 2) ) {
    char err[255];
    sprintf(err, "Angle type %i doesn't have an excited state table", type);
    error->one(FLERR, err);
  }

  return &(d_tables + tabindex[type] + (atom->estate[i]))->ang_table;
}

BondTableE::Table* AngleTableE2::get_dist_table(int i, int j, int type) {
  if(i == -1 || j == -1 || atom->estate == NULL) {
    return &(d_tables + tabindex[type])->bond_table;
  }

  if( atom->estate[i] != atom->estate[j] ) {
    error->all(FLERR, "Sanity check Failed! Requested angle between excited and non-excited atoms ");
  }

  if(( atom->estate[i] != 0 ) && ( setflag[type] != 2 ) ) {
    char err[255];
    snprintf(err, 255, "Angle type %i doesn't have an excited state table", type);
    error->one(FLERR, err);
  }

  return &(d_tables + tabindex[type] + (atom->estate[i]))->bond_table;
}

void AngleTableE2::allocate() {
  AngleTable::allocate();
  bond->allocate();
}

void AngleTableE2::bond_ev_tally(int i, int j, int nlocal, int newton_bond,
                                 double ebond, double fbond) {

  double** x = atom->x;

  double delx = x[i][0] - x[j][0];
  double dely = x[i][1] - x[j][1];
  double delz = x[i][2] - x[j][2];

  double ebondhalf,v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) energy += ebond;
      else {
        ebondhalf = 0.5*ebond;
        if (i < nlocal) energy += ebondhalf;
        if (j < nlocal) energy += ebondhalf;
      }
    }
    if (eflag_atom) {
      ebondhalf = 0.5*ebond;
      if (newton_bond || i < nlocal) eatom[i] += ebondhalf;
      if (newton_bond || j < nlocal) eatom[j] += ebondhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx*delx*fbond;
    v[1] = dely*dely*fbond;
    v[2] = delz*delz*fbond;
    v[3] = delx*dely*fbond;
    v[4] = delx*delz*fbond;
    v[5] = dely*delz*fbond;

    if (vflag_global) {
      if (newton_bond) {
        virial[0] += v[0];
        virial[1] += v[1];
        virial[2] += v[2];
        virial[3] += v[3];
        virial[4] += v[4];
        virial[5] += v[5];
      } else {
        if (i < nlocal) {
          virial[0] += 0.5*v[0];
          virial[1] += 0.5*v[1];
          virial[2] += 0.5*v[2];
          virial[3] += 0.5*v[3];
          virial[4] += 0.5*v[4];
          virial[5] += 0.5*v[5];
        }
        if (j < nlocal) {
          virial[0] += 0.5*v[0];
          virial[1] += 0.5*v[1];
          virial[2] += 0.5*v[2];
          virial[3] += 0.5*v[3];
          virial[4] += 0.5*v[4];
          virial[5] += 0.5*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        vatom[i][0] += 0.5*v[0];
        vatom[i][1] += 0.5*v[1];
        vatom[i][2] += 0.5*v[2];
        vatom[i][3] += 0.5*v[3];
        vatom[i][4] += 0.5*v[4];
        vatom[i][5] += 0.5*v[5];
      }
      if (newton_bond || j < nlocal) {
        vatom[j][0] += 0.5*v[0];
        vatom[j][1] += 0.5*v[1];
        vatom[j][2] += 0.5*v[2];
        vatom[j][3] += 0.5*v[3];
        vatom[j][4] += 0.5*v[4];
        vatom[j][5] += 0.5*v[5];
      }
    }
  }
}

void AngleTableE2::bond_ev_setup(int eflag, int vflag, int alloc) {
  int i,n;

  evflag = 1;

  eflag_either = eflag;
  eflag_global = eflag % 2;
  eflag_atom = eflag / 2;

  vflag_either = vflag;
  vflag_global = vflag % 4;
  vflag_atom = vflag / 4;

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nmax > maxeatom) {
    maxeatom = atom->nmax;
    if (alloc) {
      memory->destroy(eatom);
      memory->create(eatom,comm->nthreads*maxeatom,"angle:eatom");
    }
  }
  if (vflag_atom && atom->nmax > maxvatom) {
    maxvatom = atom->nmax;
    if (alloc) {
      memory->destroy(vatom);
      memory->create(vatom,comm->nthreads*maxvatom,6,"angle:vatom");
    }
  }

  // zero accumulators

  if (eflag_global) energy = 0.0;
  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton_bond) n += atom->nghost;
    for (i = 0; i < n; i++) eatom[i] = 0.0;
  }
  if (vflag_atom && alloc) {
    n = atom->nlocal;
    if (force->newton_bond) n += atom->nghost;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
  }
}