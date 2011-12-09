/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 Giovanni Bussi (bussi@sissa.it)

   This file is part of plumed, version 2.0.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#define PLUMED_GENERIC_CXX 1
#include "metadyn.h"

#include "tools/Communicator.h"
#include "core/ActionRegister.h"
#include "bias/Bias.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <sys/stat.h>
#include <sstream>
#include <cassert>


using namespace std;
namespace PLMD{
namespace plumed1{

//+PLUMEDOC BIAS BIASPLUMED1
/**
Use PLUMED 1.3 bias inside PLUMED2. See also \ref COLVARPLUMED1.

\par Syntax
\verbatim
BIASPLUMED1 WHAT=what [SIGMA=sigma]
BIASPLUMED1 DO=UMBRELLA KAPPA=kappa AT=at [SLOPE=slope]
BIASPLUMED1 DO=STEER KAPPA=kappa VEL=vel TO=to [FROM=from]
BIASPLUMED1 DO=ABMD KAPPA=kappa TO=to [RESTART]
BIASPLUMED1 DO=EXTERNAL FILE=file
BIASPLUMED1 ...
  DO=HILLS HEIGHT=height W_STRIDE=w_stride SIGMA=sigma
  [RESTART]
  [WELLTEMPERED=bf SIMTEMP=simtemp]
  [GRID0=nbin,min,max[,[PBC|NOPBC]]]
  [GRID1=nbin,min,max[,[PBC|NOPBC]]]
  [WRITE_GRID_FILE=file WRITE_GRID_STRIDE=stride]
  [READ_GRID_FILE=file]
... BIASPLUMED1
\endverbatim

Check PLUMED-1.3 manual and guess the meaning of each keyword.

Some PLUMED-1.3 directives are still not mapped (e.g. STEERPLAN,DAFED,...).
They can be emulated using the generic "WHAT=" syntax (see \ref COLVARPLUMED1).

\attention Still experimental


*/
//+ENDPLUMEDOC

class BiasPlumed1 :
  public bias::Bias,
  public PlumedGenericCxx
{
  string         dir;
  vector<double> fakeAtoms;
  vector<double> fakeCell;
  vector<double> fakeMass;
  vector<double> fakeCharge;
  vector<double> fakeForces;
  vector<bool> torsion;
public:
  BiasPlumed1(const ActionOptions&);
  void calculate();
  void pbc(const double*,const double*,double*)const;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(BiasPlumed1,"BIASPLUMED1")

void BiasPlumed1::registerKeywords(Keywords& keys){
   Bias::registerKeywords(keys);
   keys.use("ARG");
   keys.add("optional","WHAT","");
   keys.add("optional","SIGMA","");
   keys.add("optional","DO","");
   keys.add("optional","KAPPA","");
   keys.add("optional","AT","");
   keys.add("optional","SLOPE","");
   keys.add("optional","VEL","");
   keys.add("optional","TO","");
   keys.add("optional","FROM","");
   keys.add("optional","FILE","");
   keys.add("optional","WELLTEMPERED","");
   keys.add("optional","SIMTEMP","");
   keys.add("optional","WRITE_GRID_FILE","");
   keys.add("optional","READ_GRID_FILE","");
   keys.add("optional","WRITE_GRID_STRIDE","");
   keys.add("numbered","GRID","");
   keys.addFlag("RESTART",false,"");
}


BiasPlumed1::BiasPlumed1(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao)
{
  const unsigned n=getNumberOfArguments();

  string what;
  string todo;
  vector<double> sigma;
  parse("DO",todo);
  if(todo.length()==0){
    parse("WHAT",what);
    for(int i=0;i<what.length();i++)
      if(what[i]=='+')     what[i]=' ';
      else if(what[i]=='^') what[i]='\n';
    parseVector("SIGMA",sigma);
    assert(sigma.size()==n || sigma.size()==0);
  } else if(todo=="UMBRELLA"){
    ostringstream line;
    line.precision(14);
    vector<double> kappa,at,slope;
    parseVector("KAPPA",kappa); assert(kappa.size()==n);
    parseVector("AT",at);       assert(at.size()==n);
    parseVector("SLOPE",slope); assert(slope.size()==0 || slope.size()==n);
    for(unsigned i=0;i<n;i++){
      line<<"UMBRELLA CV "<<(i+1)<<" KAPPA "<<kappa[0]<<" AT "<<at[0];
      if(slope.size()==n) line<<" SLOPE "<<slope[i];
      line<<"\n";
    }
    what=line.str();
  } else if(todo=="STEER"){
    assert(n==1);
    ostringstream line;
    line.precision(14);
    vector<double> kappa,to,vel,from;
    parseVector("KAPPA",kappa); assert(kappa.size()==1);
    parseVector("VEL",vel);     assert(vel.size()==1);
    parseVector("TO",to);       assert(to.size()==1);
    parseVector("FROM",from);       assert(from.size()==0 || from.size()==1);
    line<<"STEER CV 1 TO "<<to[0]<<" VEL "<<vel[0]<<" KAPPA "<<kappa[0];
    if(from.size()==1) line<<" FROM "<<from[0];
    line<<"\n";
    what=line.str();
  } else if(todo=="ABMD"){
    bool restart=false;
    ostringstream line;
    line.precision(14);
    assert(n==1);
    vector<double> to,kappa;
    parseVector("KAPPA",kappa); assert(kappa.size()==n);
    parseVector("TO",to);       assert(to.size()==n);
    parseFlag("RESTART",restart);
    line<<"ABMD CV 1 TO "<<to[0]<<" KAPPA "<<kappa[0];
    if(restart) line<<" RESTART";
    line<<"\n";
    what=line.str();
  } else if(todo=="EXTERNAL"){
    ostringstream line;
    line.precision(14);
    string file("");
    parse("FILE",file);
    assert(file.length()>0);
    line<<"EXTERNAL FILENAME "<<file;
    line<<" NCV "<<n<< " CV";
    for(unsigned i=0;i<n;i++) line<<" "<<(i+1);
    line<<"\n";
    what=line.str();
  } else if(todo=="HILLS"){
    double height=-1.0; parse("HEIGHT",height);     assert(height>0);
    int w_stride=0;     parse("W_STRIDE",w_stride); assert(w_stride>0);
    ostringstream line;
    line.precision(14);
    line<<"HILLS HEIGHT "<<height<<" W_STRIDE "<<w_stride;
    bool restart=false;
    parseFlag("RESTART",restart);
    if(restart)line<<" RESTART";
    line<<"\n";
    double biasfactor=-1;
    parse("WELLTEMPERED",biasfactor);
    if(biasfactor>0.0){
      line<<"WELLTEMPERED BIASFACTOR "<<biasfactor;
      double simtemp=-1;
      parse("SIMTEMP",simtemp);
      if(simtemp>0.0) line<<" SIMTEMP "<<simtemp;
      line<<"\n";
    }
    bool mwalwers=false;
    parseFlag("WALKERS",mwalwers);
    if(mwalwers){
       line<<"MULTIPLE_WALKERS";
       std::string dir;
       parse("WALKERS_DIR",dir);
       assert(dir.length()>0);
       int stride=0;
       parse("WALKERS_STRIDE",stride);
       assert(stride>0);
       int nwalkers=0;
       parse("WALKERS_NUM",nwalkers);
       assert(nwalkers>0);
       int id=-1;
       parse("WALKERS_ID",id);
       assert(id>=0);
       line<<" HILLS_DIR "<<dir<<" R_STRIDE "<<stride<<" NWALKERS "<<nwalkers<<" ID "<<id<<"\n";
    }
    for(int i=0;i<n;i++){
      bool pbc;
      double min,max;
      int nbin;
      string s; Tools::convert(i,s); 
      vector<string> grid;
      parseVector("GRID"+s,grid);
      if(grid.size()==0)continue;
      assert(grid.size()>=1 && grid.size()<=4);
      Tools::convert(grid[0],nbin);
      Tools::convert(grid[1],min);
      Tools::convert(grid[2],max);
      pbc=false;
      if(grid.size()==4){
        if(grid[3]=="PBC") pbc=true;
        else if(grid[3]=="NOPBC") pbc=false;
        else assert(0);
      }
      line<<"GRID CV "<<(i+1)<<" NBIN "<<nbin<<" MIN "<<min<<" MAX "<<max;
      if(pbc) line<<" PBC";
      line<<"\n";
    }

    parseVector("SIGMA",sigma);
    assert(sigma.size()==n);

    string gfile;
    parse("WRITE_GRID_FILE",gfile);
    if(gfile.length()>0){
      int str=-1;
      parse("WRITE_GRID_STRIDE",str);
      assert(str>0);
      line<<"WRITE_GRID FILENAME "<<gfile<<" W_STRIDE "<<str<<"\n";
    }

    string rfile;
    parse("READ_GRID_FILE",rfile);
    if(rfile.length()>0){
      line<<"READ_GRID FILENAME "<<rfile<<"\n";
    }
    what=line.str();

// STILL MISSING IN HILLS: INVERT, ...

  } else assert(0);

// OTHER MISSING DIRECTIVES: STEERPLAN, UWALL, LWALL, COMMITMENT

  checkRead();

  if(sigma.size()>0){
    log.printf("  with sigmas");
    for(int i=0;i<n;i++) log.printf(" %f",sigma[i]);
    log.printf("\n");
  }

  log.printf("  doing:\n");
  log.printf("%s\n",what.c_str());

// We use 4 fake atoms per CV (see below)
  fakeAtoms.assign(12*n,0.0);
  fakeForces.assign(12*n,0.0);
  fakeMass.assign(4*n,1.0);
  fakeCharge.assign(4*n,1.0);
  fakeCell.assign(9,0.0);
  torsion.assign(n,false);

  for(unsigned i=0;i<n;i++){
    Value*v=getArguments()[i];
    if(v->isPeriodic()){
      double min,max;
      v->getDomain(min,max);
      if(pow(min+3.14159265358979323844,2.0)<1000*epsilon && pow(max-3.14159265358979323844,2.0)<1000*epsilon) torsion[i]=true;
      else assert(0);
    }
  }

  char* t="plumed1.dat";

  dir="plumed1."+getLabel()+plumed.getSuffix();
  log.printf("  in directory %s\n",dir.c_str());

  mkdir(dir.c_str(),0777);


  char pwd[10000];
  getcwd(pwd,10000);
  chdir(dir.c_str());

// I use the std version to avoid suffix appending
  FILE* tmp=std::fopen(t,"w");
  fprintf(tmp,"# FAKE VARIABLES:\n");
// This variables are added to trick plumed1.
// The rule is:
// - for non-periodic CV, we define a fake atom an use its X coordinate
// - for periodic CVs with period [-pi,pi] we define four atoms forming
//   a dihedral angle equal to the wanted CV and use TORSION as a CV
  for(unsigned i=0;i<n;i++){
    if(torsion[i]) fprintf(tmp,"TORSION LIST %d %d %d %d",4*i+1,4*i+2,4*i+3,4*i+4);
    else           fprintf(tmp,"POSITION LIST %d DIR X",4*i+1);
    if(sigma.size()>0) fprintf(tmp," SIGMA %f",sigma[i]);
    fprintf(tmp,"\n");
  }
  fprintf(tmp,"# FREE ENERGY METHODS:\n");
  fprintf(tmp,"%s\n",what.c_str());
  std::fclose(tmp);

  double dt=getTimeStep();
  int nn=n;
#ifdef __PLUMED_MPI
  set_comm(comm.Get_comm());
#endif
  init_metadyn(4*nn,dt,&fakeMass[0], &fakeCharge[0],plumed.getAtoms().getKBoltzmann(),1.0,t);
  chdir(pwd);
  addComponent("bias");
  log<<"  Bibliography "<<plumed.cite("Bonomi, Branduardi, Bussi, Camilloni, Provasi, Raiteri, "
                                      "Donadio, Marinelli, Pietrucci, Broglia, and Parrinello, "
                                      "Comp. Phys. Comm. 180, 1961 (2009)")<<"\n";
  log<<"  Have a look to the "<<dir<<"/plumed.out file for other relevant references\n";
}

void BiasPlumed1::calculate(){
  const unsigned n=getNumberOfArguments();
  for(unsigned i=0;i<n;i++){
    if(torsion[i]){
      fakeAtoms[i*12+0]=1.0;
      fakeAtoms[i*12+2]=-1.0;
      fakeAtoms[i*12+5]=-1.0;
      fakeAtoms[i*12+8]=1.0;
      fakeAtoms[i*12+9]=cos(getArgument(i));
      fakeAtoms[i*12+10]=sin(getArgument(i));
      fakeAtoms[i*12+11]=1.0;
    } else {
      fakeAtoms[i*12]=getArgument(i);
    }
  }
  for(unsigned i=0;i<12*n;i++) fakeForces[i]=0.0;
  int step=getStep();
  double energy=0.0;
  char pwd[10000];
  getcwd(pwd,10000);
  chdir(dir.c_str());
  meta_force_calculation(step,&fakeAtoms[0],&fakeForces[0],&energy);

  chdir(pwd);
  for(unsigned i=0;i<n;i++) {
    if(torsion[i]){
      double f;
      f=-(fakeForces[i*12+9])*sin(getArgument(i))+(fakeForces[i*12+10])*cos(getArgument(i));
      setOutputForce(i,f);
    } else {
      setOutputForce(i,fakeForces[i*12]);
    }
  }
  getPntrToComponent("bias")->set(energy);
}

void BiasPlumed1::pbc(const double*p1,const double*p2,double*d)const{
  d[0]=p1[0]-p2[0];
  d[1]=p1[1]-p2[1];
  d[2]=p1[2]-p2[2];
}


}
}


