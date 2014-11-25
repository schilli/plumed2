/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

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
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/File.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Units.h"
#include <cstdio>
#include "core/SetupMolInfo.h"
#include "core/ActionSet.h"
#include "MultiColvarBase.h"
#include "tools/Grid.h"
#include "tools/KernelFunctions.h"
#include "vesselbase/ActionWithInputVessel.h"

using namespace std;

namespace PLMD
{
namespace multicolvar {

//+PLUMEDOC ANALYSIS MULTICOLVARDENS
/*
Dump atom positions and multicolvar on a file.

\par Examples


*/
//+ENDPLUMEDOC

class MultiColvarDensity :
  public ActionPilot,
  public ActionAtomistic,
  public vesselbase::ActionWithInputVessel
{
  std::string kerneltype;
  bool nomemory;
  double norm;
  unsigned rstride;
  std::string filename;
  Grid* gg;
  MultiColvarBase* mycolv; 
  std::vector<double> bw;
  std::vector<unsigned> directions;
public:
  MultiColvarDensity(const ActionOptions&);
  ~MultiColvarDensity();
  static void registerKeywords( Keywords& keys );
  void calculate(){}
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ){ plumed_error(); }
  void apply(){}
  void update();
};

PLUMED_REGISTER_ACTION(MultiColvarDensity,"MULTICOLVARDENS")

void MultiColvarDensity::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithInputVessel::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the grid");
  keys.add("compulsory","RUN","the frequency with which the density profile is written out");
  keys.add("atoms","ORIGIN","we will use the position of this atom as the origin");
  keys.add("compulsory","DIR","the direction in which to calculate the density profile");
  keys.add("compulsory","NBINS","the number of bins to use to represent the density profile");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using.  More details on  the kernels available "
                                            "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("compulsory","OFILE","density","the file on which to write the profile");
  keys.addFlag("NOMEMORY",false,"do a block averaging rather than a cumulative average");
}

MultiColvarDensity::MultiColvarDensity(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  ActionWithInputVessel(ao),
  norm(0)
{

  std::vector<AtomNumber> atom;
  parseAtomList("ORIGIN",atom);
  if( atom.size()!=1 ) error("should only be one atom specified");
  log.printf("  origin is at position of atom : %d\n",atom[0].serial() );

  readArgument("store");
  mycolv = dynamic_cast<MultiColvarBase*>( getDependencies()[0] );
  plumed_assert( getDependencies().size()==1 ); 
  if(!mycolv) error("action labeled " + mycolv->getLabel() + " is not a multicolvar");

  parse("OFILE",filename); parse("RUN",rstride);
  if(filename.length()==0) error("name out output file was not specified");

  std::vector<unsigned> nbins; parseVector("NBINS",nbins); parseFlag("NOMEMORY",nomemory);
  parse("KERNEL",kerneltype); parseVector("BANDWIDTH",bw); 
  std::string direction; parse("DIR",direction);
  log.printf("  calculating density profile along ");
  if( direction=="x" ){
      if( bw.size()!=1 && nbins.size()!=1 ) error("BANDWIDTH or NBINS has wrong dimensionality");
      log.printf("x axis");
      directions.resize(1); directions[0]=0;
  } else if( direction=="y" ){
      if( bw.size()!=1 && nbins.size()!=1 ) error("BANDWIDTH or NBINS has wrong dimensionality");
      log.printf("y axis");
      directions.resize(1); directions[0]=1;
  } else if( direction=="z" ){
      if( bw.size()!=1 && nbins.size()!=1 ) error("BANDWIDTH or NBINS has wrong dimensionality");
      log.printf("z axis");
      directions.resize(1); directions[0]=2;
  } else if( direction=="xy" ){
      if( bw.size()!=2 && nbins.size()!=2 ) error("BANDWIDTH or NBINS has wrong dimensionality");
      log.printf("x and y axes");
      directions.resize(2); directions[0]=0; directions[1]=1;
  } else if( direction=="xz" ){
      if( bw.size()!=2 && nbins.size()!=2 ) error("BANDWIDTH or NBINS has wrong dimensionality");
      log.printf("x and z axes");
      directions.resize(2); directions[0]=0; directions[1]=2;
  } else if( direction=="yz" ){
      if( bw.size()!=2 && nbins.size()!=2 ) error("BANDWIDTH or NBINS has wrong dimensionality");
      log.printf("y and z axis");
      directions.resize(2); directions[0]=1; directions[1]=2;
  } else if( direction=="xyz" ){
      if( bw.size()!=3 && nbins.size()!=3 ) error("BANDWIDTH or NBINS has wrong dimensionality");
      log.printf("x, y and z axes");
      directions.resize(3); directions[0]=0; directions[1]=1; directions[2]=2;
  } else {
     error( direction + " is not valid gradient direction");
  } 
  log.printf(" for colvars calculated by action %s \n",mycolv->getLabel().c_str() );

  checkRead(); requestAtoms(atom); 
  // Stupid dependencies cleared by requestAtoms - why GBussi why? That's got me so many times
  addDependency( mycolv );

  std::vector<bool> pbc(nbins.size());
  std::vector<std::string> args(nbins.size()), gmin(nbins.size()), gmax(nbins.size());;
  for(unsigned i=0;i<directions.size();++i){
      gmin[i]="-0.5"; gmax[i]="0.5"; pbc[i]=true;
      if( directions[i]==0 ) args[i]="x";
      else if( directions[i]==1 ) args[i]="y";
      else if( directions[i]==2 ) args[i]="z";
      else plumed_error();
  }   
  
  if( plumed.getRestart() ){   
     error("restarting of MultiColvarDensity is not yet implemented");
  } else {
     // Setup the grid
     std::string funcl=mycolv->getLabel() + ".dens";
     gg = new Grid(funcl,args,gmin,gmax,nbins,true,true,true,pbc,gmin,gmax);
  }
}

MultiColvarDensity::~MultiColvarDensity(){
  delete gg;
}

void MultiColvarDensity::update(){

  Vector origin = getPosition(0); std::vector<double> pp( directions.size() );
  std::vector<double> cvals( mycolv->getNumberOfQuantities()-4 );
  for(unsigned i=0;i<mycolv->getFullNumberOfTasks();++i){
      mycolv->getValueForTask( i, cvals );
      Vector apos = pbcDistance( mycolv->getCentralAtomPosition(i), origin );
      Vector fpos = getPbc().realToScaled( apos );
      for(unsigned j=0;j<directions.size();++j) pp[j]=fpos[ directions[j] ];
      KernelFunctions kernel( pp, bw, kerneltype, false, cvals[0], true );
      gg->addKernel( kernel );
      norm += 1.0;    // This should be replaced by the proper weight
  }

  // Output and reset the counter if required
  if( getStep()%rstride==0 ){  // && getStep()>0 ){
      // Normalise prior to output
      gg->scaleAllValuesAndDerivatives( 1.0 / norm );

      OFile gridfile; gridfile.link(*this); gridfile.setBackupString("analysis");
      gridfile.open( filename ); gg->writeToFile( gridfile ); gridfile.close();

      if( nomemory ){ 
        gg->clear(); norm=0.0; 
      } else {
        // Unormalise after output
        gg->scaleAllValuesAndDerivatives( norm );
      }
  }

}

}
}
