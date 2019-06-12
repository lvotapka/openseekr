/*
   Copyright 2019 by Andy Stokely
   All rights reserved
   
 */

#include "SeekrForce.h"
#include "internal/SeekrForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include <iostream>

using namespace SeekrPlugin;
using namespace OpenMM;
using namespace std;

SeekrForce::SeekrForce() {
}

int SeekrForce::getNumPlanarZMilestones() const {
  return planarZMilestones.size();
}

int SeekrForce::getNumPlanarZAtomIndices() const {
  int total = 0;
  for(std::vector<PlanarZMilestoneInfo>::const_iterator iter=planarZMilestones.begin();	
      iter != planarZMilestones.end(); ++iter) {
  total += iter->numIndices1;
  total += iter->numIndices2;
  }
  return total;
}

int SeekrForce::getPlanarZNumIndices(int forceIndex, int molecule) const {
  if (molecule == 1) {
    return planarZMilestones[forceIndex].numIndices1;
  } else if (molecule == 2) {
    return planarZMilestones[forceIndex].numIndices2;
  } else {
    throw OpenMMException("Molecule index must be 1 or 2");
  }
}

float SeekrForce::getPlanarZOffset(int forceIndex, int milestone_id) const {
  if (milestone_id == 1) {
    return planarZMilestones[forceIndex].offset1;
  } else if (milestone_id == 2) {
    return planarZMilestones[forceIndex].offset2;
  } else if (milestone_id == 3) {
    return planarZMilestones[forceIndex].offset3;
  } else {
    throw OpenMMException("Molecule index must be 1, 2, or 3");
  }
}

void SeekrForce::getPlanarZMilestoneAtoms(int forceIndex, int atomIndex, int& atom_id, int molecule) const {
cout << "ForceIndex: "<< forceIndex << " Atom Index: " << atomIndex << endl;
cout << "Get PlanarZ Num Indicies: " << getPlanarZNumIndices(forceIndex, molecule) << endl; 
  if (molecule == 1) {
    atom_id = planarZMilestones[forceIndex].atomIndices1[atomIndex];
  } else if (molecule == 2) {
    atom_id = planarZMilestones[forceIndex].atomIndices2[atomIndex];
  } else {
    throw OpenMMException("Molecule index must be 1 or 2");
  }
}

bool SeekrForce::getEndOnMiddleCrossing() const {
  return endOnMiddleCrossing;
}

std::string SeekrForce::getDataFileName(int forceIndex) const {
  return planarZMilestones[forceIndex].dataFileName;
}

void SeekrForce::addPlanarZMilestone(int numIndices1, int numIndices2, float offset1, 
                              float offset2, float offset3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2, bool endOnMiddleCrossingArg,
                              std::string dataFileName) 
{
  /*cout << "Adding planarZ milestone. numIndices1: " << numIndices1 << " numIncices: ";
  cout << numIndices2 << " radii: (" << offset1 << ", " << offset2 << ", " << offset3;
  cout << "). atomIndices1: [";
  for (int i=0; i < numIndices1; i++) {
    cout << atomIndices1[i] << " ";
  }
  cout << "] atomIndices2: [";
  for (int i=0; i < numIndices2; i++) {
    cout << atomIndices2[i] << " ";
  }
  cout << "]\n";
  cout << "dataFileName:" << dataFileName << "\n";*/
  planarZMilestones.push_back(
      PlanarZMilestoneInfo(numIndices1, numIndices2, offset1, offset2, offset3, atomIndices1, atomIndices2, dataFileName));
  endOnMiddleCrossing = endOnMiddleCrossingArg;
  
  
}

void SeekrForce::modifyPlanarZMilestone(int forceIndex, int numIndices1, int numIndices2, float offset1, 
                              float offset2, float offset3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2, bool endOnMiddleCrossingArg,
                              std::string dataFileName) 
{
  /*int i;
  planarZMilestone.numIndices1 = numIndices1;
  planarZMilestone.numIndices2 = numIndices2;
  planarZMilestone.offset1 = offset1;
  planarZMilestone.offset2 = offset2;
  planarZMilestone.offset3 = offset3;
  for (i=0; i<numIndices1; i++) {
    planarZMilestone.atomIndices1[i] = atomIndices1[i];
  }
  for (i=0; i<numIndices2; i++) {
    planarZMilestone.atomIndices2[i] = atomIndices2[i];
  }
  std::cout << "Modifying planarZ milestone. Indices: " << numIndices1 << ", " << numIndices2;
  std::cout << ". Radii: " << offset1 << ", " << offset2 << ", " << offset3 << ".\n";
  */
  planarZMilestones[forceIndex] = PlanarZMilestoneInfo(numIndices1, numIndices2, offset1, offset2, offset3, atomIndices1, atomIndices2, dataFileName);
  endOnMiddleCrossing = endOnMiddleCrossingArg;
  
}


ForceImpl* SeekrForce::createImpl() const {
    return new SeekrForceImpl(*this);
}

void SeekrForce::updateParametersInContext(Context& context) {
    dynamic_cast<SeekrForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
