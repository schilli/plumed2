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
#include "tools/Vector.h"
#include "tools/AtomNumber.h"
#include "tools/Tools.h"
#include "core/Atoms.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/SetupMolInfo.h"

#include <vector>
#include <string>
#include <limits>

using namespace std;

namespace PLMD {
namespace generic{

//+PLUMEDOC GENERIC WRAPAROUND
/*
Rebuild periodic boundary conditions around chosen atoms.

Modify position of atoms indicated by ATOMS by shifting them by lattice vectors so that they are
as close as possible to the atoms indicated by AROUND. More precisely, for every atom i
in the ATOMS list the following procedure is performed:
- The atom j among those in the AROUND list is searched that is closest to atom i.
- The atom i is replaced with its periodic image that is closest to atom j.

This action works similarly to \ref WHOLEMOLECULES in that it replaces atoms coordinate. Notice that only
atoms specified with ATOMS are replaced, and that, at variance with \ref WHOLEMOLECULES,
the order in which atoms are specified is irrelevant.

\attention This directive modifies the stored position at the precise moment it is executed. This means that only collective variables which are below it in the input script will see the corrected positions. As a general rule, put it at the top of the input file. Also, unless you know exactly what you are doing, leave the default stride (1), so that this action is performed at every MD step.

Consider that the computational cost grows with the product
of the size of the two lists (ATOMS and AROUND), so that this action can become very expensive.
If you are using it to analyse a trajectory this is usually not a big problem. If you use it to
analyze a simulation on the fly, e.g. with \ref DUMPATOMS to store a properly wrapped trajectory,
consider the possibility of using the STRIDE keyword here (with great care).
\par Examples

This command instructs plumed to move all the ions to their periodic image that is as close as possible to
the rna group.

\verbatim
rna: GROUP ATOMS=1-100
ions: GROUP ATOMS=101-110
# first make the rna molecule whole
WHOLEMOLECULES ENTITY0=rna
WRAPAROUND ATOMS=ions AROUND=rna
DUMPATOMS FILE=dump.xyz ATOMS=rna,ions
\endverbatim
(see also \ref GROUP and \ref DUMPATOMS)

This command instructs plumed to center all atoms around the center of mass of a solute molecule.

\verbatim
solute: GROUP ATOMS=1-100
all: GROUP ATOMS=1-1000
# first make the solute whole
# this is necessary to compute its center of mass properly
WHOLEMOLECULES ENTITY0=solute
# center of the solute:
com: COM ATOMS=solute
# notice that we wrap around a single atom. this should be fast
WRAPAROUND ATOMS=all AROUND=com
DUMPATOMS FILE=dump.xyz ATOMS=all
\endverbatim
(see also \ref GROUP \ref COM \ref DUMPATOMS)

Notice that whereas \ref WHOLEMOLECULES is designed to make molecules whole,
\ref WRAPAROUND can easily break molecules. In the last example,
if solvent (atoms 101-1000) is made e.g. of water, then water
molecules could be broken by \ref WRAPAROUND (hydrogen could end up
in an image and oxygen in another one).

*/
//+ENDPLUMEDOC


class WrapAround:
  public ActionPilot,
  public ActionAtomistic
{
  vector<AtomNumber> atoms;
  vector<AtomNumber> reference;
public:
  WrapAround(const ActionOptions&ao);
  static void registerKeywords( Keywords& keys );
  void calculate();
  void apply(){}
};

PLUMED_REGISTER_ACTION(WrapAround,"WRAPAROUND")

void WrapAround::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which molecules are reassembled.  Unless you are completely certain about what you are doing leave this set equal to 1!");
  keys.add("atoms","AROUND","reference ions");
  keys.add("atoms","ATOMS","atoms numbers");

}

WrapAround::WrapAround(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionAtomistic(ao)
{
  parseAtomList("ATOMS",atoms);
  parseAtomList("AROUND",reference);

  log.printf("  atoms in reference :");
  for(unsigned j=0;j<reference.size();++j) log.printf(" %d",reference[j].serial() );
  log.printf("\n");
  log.printf("  atoms to be wrapped :");
  for(unsigned j=0;j<atoms.size();++j) log.printf(" %d",atoms[j].serial() );
  log.printf("\n");
  
  checkRead();

  Tools::removeDuplicates(atoms);
  Tools::removeDuplicates(reference);

  vector<AtomNumber> merged(atoms.size()+reference.size());
  merge(atoms.begin(),atoms.end(),reference.begin(),reference.end(),merged.begin());
  Tools::removeDuplicates(merged);
  requestAtoms(merged);
  doNotRetrieve();
  doNotForce();
}

void WrapAround::calculate(){
  for(unsigned i=0;i<atoms.size();++i){
    Vector & first (modifyPosition(atoms[i]));
    double mindist2=std::numeric_limits<double>::max();
    int closest=-1;
    for(unsigned j=0;j<reference.size();++j){
      Vector & second (modifyPosition(reference[j]));
      Vector distance=pbcDistance(first,second);
      double distance2=modulo2(distance);
      if(distance2<mindist2){
        mindist2=distance2;
        closest=j;
      }
    }
    plumed_massert(closest>=0,"closest not found");
    Vector & second (modifyPosition(reference[closest]));
    first=second+pbcDistance(second,first);
  }
}



}

}
