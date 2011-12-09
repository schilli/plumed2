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
#include "colvar/Colvar.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <sys/stat.h>
#include <cassert>


using namespace std;
namespace PLMD{
namespace plumed1{

//+PLUMEDOC COLVAR COLVARPLUMED1
/**
Use PLUMED 1.3 collective variables inside PLUMED2.
The PLUMED_GENERIC_CXX interface of PLUMED 1.3 is ran in
a separate directory using some fake plumed.dat input.

\par Syntax
\verbatim
COLVARPLUMED1 ATOMS=atoms WHAT=what [NOPBC] [PERIODIC=periodic] [KEEPNUMBERS]
\endverbatim
- atoms is the ordered list of atoms passed to plumed 1.
- what is a string which is used in the plumed 1 input file. It is parsed as follows:
  - '+' are replaced with spaces
  - '^' are replaced with newlines
- if NOPBC is present, current cell is not passed to plumed 1, which just sees a
  gas phase simulation
- PERIODIC can be used to instruct plumed if this variable is periodic or not. It should be
  provided if the variable has to be biased, and can either be
  - PERIODIC=NO for a non periodic variable
  - PERIODIC=a,b for a variable with periodicity domain [a,b]
  Notice that this usually makes only sense for TORSION or PUCKERING
- if KEEPNUMBERS is present, atoms indexing inside plumed1 is consistent with
  the entire system. This helps for RMSD based variables, or, in general, to
  reuse plumed groups. Note that also when using this option it is 
  compulsory to select at least the required atoms with the ATOMS keyword,
  and that if many more atoms are selected there will be a considerable
  overhead in parallel simulations (since all the ATOMS are passed to plumed1).
- if ATOMS is absent, that numbering is kept consistent (as in KEEPNUMBERS)
  and list of atoms is automatically created. This is the preferred approach,
  albeit it does not allow accessing virtual atoms.

\par Example
\verbatim
COLVARPLUMED1 ...
  NOPBC
  PERIODIC=NO
  LABEL=rmsd
  ATOMS=339,349,369,379,403,413,436,446,469,479,2358,2368,2391,2401,2421,2431,2451,2461,2482,2492
# filename is indicated with a leading ../ so as to access to in the parent directory
# this is necessary as plumed-1 runs in a separate directory
# also notice that atom numbers in the pdb file are relative to the list of passed atoms
# unless the KEEPNUMBERS flag is used
  WHAT=TARGETED+TYPE+RMSD+FRAMESET+../refc.pdb
... COLVARPLUMED1
COLVARPLUMED1 ...
# it is compulsory to state periodicity to add a bias on this variable, similarly to functions
  PERIODIC=-PI,+PI
  LABEL=torsion
  ATOMS=60,65,10,20
# atom numbers passed here are *relative* numbers among the passed list
# unless the KEEPNUMBERS flag is used
  WHAT=TORSION+LIST+1+2+3+4
... COLVARPLUMED1
\endverbatim

\note It may be useful to do a full mapping of the keyword (at least
those documented in the PLUMED 1.3 manual) so as to allow the user to.
avoid the nasty 'WHAT=' syntax.

\attention Still experimental


*/
//+ENDPLUMEDOC

class ColvarPlumed1 :
  public colvar::Colvar,
  public PlumedGenericCxx
{
  string         dir;
  vector<AtomNumber> atoms;
  vector<double> fakeAtoms;
  vector<double> fakeCell;
  vector<double> fakeMass;
  vector<double> fakeCharge;
  vector<double> fakeForces;
  bool nopbc;
  bool first;
  bool keepnumbers;
  bool automatic;
public:
  ColvarPlumed1(const ActionOptions&);
  void calculate();
  void pbc(const double*,const double*,double*)const;
  virtual void prepare();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ColvarPlumed1,"COLVARPLUMED1")

void ColvarPlumed1::registerKeywords(Keywords& keys){
   Colvar::registerKeywords(keys);
   keys.add("compulsory","WHAT","");
   keys.add("atoms","ATOMS","the list of atoms involved in this collective variable");
   keys.add("optional","PERIODIC","");
   keys.addFlag("KEEPNUMBERS",false,"");
}


ColvarPlumed1::ColvarPlumed1(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
nopbc(false),
first(true),
keepnumbers(false),
automatic(false)
{
  addValueWithDerivatives();

  string what;
  parse("WHAT",what);
  vector<double> sigma;
  parseAtomList("ATOMS",atoms);
  vector<string> period;
  double min(0),max(0);
  parseVector("PERIODIC",period);
  if(period.size()==0){
  }else if(period.size()==1 && period[0]=="NO"){
    setNotPeriodic();
  } else if(period.size()==2){
    setPeriodic(period[0],period[1]);
  } else assert(0);
  parseFlag("NOPBC",nopbc);
  if(nopbc){
    log.printf("  atoms PBC are switched off\n");
    log.printf("  this allows a proper virial calculation\n");
  } else {
    log.printf("  atoms PBC are switched on\n");
    log.printf("  this does not allow proper virial calculation\n");
  }
  parseFlag("KEEPNUMBERS",keepnumbers);
  if(atoms.size()==0) automatic=true;
  if(automatic) keepnumbers=true;
  checkRead();

  if(automatic){
    log.printf("  atom numbers are automatically got from plumed 1\n");
  } else {
    log.printf("  with atoms");
    for(int i=0;i<atoms.size();i++) log.printf(" %d",atoms[i].serial());
    log.printf("\n");
  }
  if(keepnumbers) log.printf("  atom numbers are kept consistent in plumed 1\n");

  log.printf("  doing %s\n",what.c_str());
  log.printf("  which is better read as:\n");
  for(int i=0;i<what.length();i++)
    if(what[i]=='+')     what[i]=' ';
    else if(what[i]=='^') what[i]='\n';
  log.printf("%s\n",what.c_str());

  unsigned n=atoms.size();
  unsigned nn=atoms.size();
  if(automatic){
    nn=plumed.getAtoms().getNatoms();
    atoms.resize(nn);
    for(unsigned i=0;i<nn;i++) atoms[i]=AtomNumber::index(i);
    log.printf("  virtual plumed is receiving %d atoms\n",nn);
    log.printf("  some of which are automatically set at each step\n");
  } else if(keepnumbers){
// n is set to the largest atom index required
    nn=0;
    for(int i=0;i<atoms.size();i++) if(atoms[i].index()>nn) nn=atoms[i].index();
    nn++; // +1 (for array boundary)
    log.printf("  virtual plumed is receiving %d atoms\n",nn);
    log.printf("  %d of which are actually set\n",n);
  }

  fakeAtoms.assign(3*nn,0.0);
  fakeForces.assign(3*nn,0.0);
  fakeMass.assign(nn,1.0);
  fakeCharge.assign(nn,1.0);
  fakeCell.assign(9,0.0);
  char* t="plumed1.dat";

  dir="plumed1."+getLabel()+plumed.getSuffix();
  log.printf("  in directory %s\n",dir.c_str());


  mkdir(dir.c_str(),0777);

  char pwd[10000];
  getcwd(pwd,10000);
  chdir(dir.c_str());

// I use the std version to avoid suffix appending
  FILE* tmp=std::fopen(t,"w");

  fprintf(tmp,"# COLLECTIVE VARIABLE:\n");
  fprintf(tmp,"%s\n",what.c_str());
  fprintf(tmp,"# FAKE UMBRELLA SAMPLING:\n");
  fprintf(tmp,"UMBRELLA CV 1 KAPPA 0.0 SLOPE 1.0 AT 0.0\n");
  std::fclose(tmp);

  chdir(pwd);
// CALL TO init_metadyn IS DEFERRED TO FIRST calculate(), AS HERE WE DO NOT 
// KNOW MASSES/CHARGES YET
  requestAtoms(atoms);

  log<<"  Bibliography "<<plumed.cite("Bonomi, Branduardi, Bussi, Camilloni, Provasi, Raiteri, "
                                      "Donadio, Marinelli, Pietrucci, Broglia, and Parrinello, "
                                      "Comp. Phys. Comm. 180, 1961 (2009)")<<"\n";
  log<<"  Have a look to the "<<dir<<"/plumed.out file for other relevant references\n";
}

void ColvarPlumed1::calculate(){
  unsigned nn=fakeMass.size();
  unsigned n=atoms.size();
  for(unsigned i=0;i<n;i++){
    unsigned j=i; if(keepnumbers) j=atoms[i].index();
    fakeAtoms[j*3+0]=getPosition(i)[0];
    fakeAtoms[j*3+1]=getPosition(i)[1];
    fakeAtoms[j*3+2]=getPosition(i)[2];
  }
  for(unsigned j=0;j<3*nn;j++) fakeForces[j]=0.0;
  int step=getStep();
  double energy=0.0;
  char pwd[10000];
  getcwd(pwd,10000);
  chdir(dir.c_str());
  if(first){
    char* t="plumed1.dat";
    double dt=getTimeStep();
    for(unsigned i=0;i<n;i++) {
      unsigned j=i; if(keepnumbers) j=atoms[i].index();
      fakeMass[j]=getMass(i);
      fakeCharge[j]=getCharge(i);
    }
#ifdef __PLUMED_MPI
    set_comm(comm.Get_comm());
#endif
    double factor_angstrom=0.1/plumed.getAtoms().getUnits().getLength();
    init_metadyn(nn,dt,&fakeMass[0], &fakeCharge[0],plumed.getAtoms().getKBoltzmann(),factor_angstrom,t);
    first=false;
  }
  meta_force_calculation(step,&fakeAtoms[0],&fakeForces[0],&energy);

  chdir(pwd);
  Tensor virial;
  for(unsigned i=0;i<n;i++){
    unsigned j=i; if(keepnumbers) j=atoms[i].index();
    Vector p(-fakeAtoms[3*j],-fakeAtoms[3*j+1],-fakeAtoms[3*j+2]);
    Vector d(-fakeForces[3*j],-fakeForces[3*j+1],-fakeForces[3*j+2]);
    setAtomsDerivatives(i,d);
    if(nopbc) virial+=Tensor(p,d);
  }
  setBoxDerivatives(virial);
  setValue(energy);
}

void ColvarPlumed1::pbc(const double*p1,const double*p2,double*d)const{
  Vector p1v(p1[0],p1[1],p1[2]);
  Vector p2v(p2[0],p2[1],p2[2]);
  Vector dv;
  if(!nopbc) dv=pbcDistance(p2v,p1v);
  else      dv=delta(p2v,p1v);
  d[0]=dv[0];
  d[1]=dv[1];
  d[2]=dv[2];
}

void ColvarPlumed1::prepare(){
 if(first) return;
 if(!automatic) return;
 const std::set<unsigned> & at(PlumedGenericCxx::getAtoms());
 atoms.clear();
 for(std::set<unsigned>::const_iterator p=at.begin();p!=at.end();++p){
    atoms.push_back(AtomNumber::index(*p));
 }
 requestAtoms(atoms);
}



}
}


