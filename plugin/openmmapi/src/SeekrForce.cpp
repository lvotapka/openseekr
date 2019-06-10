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

int SeekrForce::getNumPlanarMilestones() const {
  return planarMilestones.size();
}

int SeekrForce::getNumPlanarAtomIndices() const {
  int total = 0;
  for(std::vector<PlanarMilestoneInfo>::const_iterator iter=planarMilestones.begin();	
      iter != planarMilestones.end(); ++iter) {
  total += iter->numIndices1;
  total += iter->numIndices2;
  }
  return total;
}

int SeekrForce::getPlanarNumIndices(int forceIndex, int molecule) const {
  if (molecule == 1) {
    return planarMilestones[forceIndex].numIndices1;
  } else if (molecule == 2) {
    return planarMilestones[forceIndex].numIndices2;
  } else {
    throw OpenMMException("Molecule index must be 1 or 2");
  }
}

float SeekrForce::getPlanarLength(int forceIndex, int milestone_id) const {
  if (milestone_id == 1) {
    return planarMilestones[forceIndex].length1;
  } else if (milestone_id == 2) {
    return planarMilestones[forceIndex].length2;
  } else if (milestone_id == 3) {
    return planarMilestones[forceIndex].length3;
  } else {
    throw OpenMMException("Molecule index must be 1, 2, or 3");
  }
}

void SeekrForce::getPlanarMilestoneAtoms(int forceIndex, int atomIndex, int& atom_id, int molecule) const {
cout << "ForceIndex: "<< forceIndex << " Atom Index: " << atomIndex << endl;
cout << "Get Planar Num Indicies: " << getPlanarNumIndices(forceIndex, molecule) << endl; 
  if (molecule == 1) {
    atom_id = planarMilestones[forceIndex].atomIndices1[atomIndex];
  } else if (molecule == 2) {
    atom_id = planarMilestones[forceIndex].atomIndices2[atomIndex];
  } else {
    throw OpenMMException("Molecule index must be 1 or 2");
  }
}

bool SeekrForce::getEndOnMiddleCrossing() const {
  return endOnMiddleCrossing;
}

std::string SeekrForce::getDataFileName(int forceIndex) const {
  return planarMilestones[forceIndex].dataFileName;
}

void SeekrForce::addPlanarMilestone(int numIndices1, int numIndices2, float length1, 
                              float length2, float length3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2, bool endOnMiddleCrossingArg,
                              std::string dataFileName) 
{
  /*cout << "Adding planar milestone. numIndices1: " << numIndices1 << " numIncices: ";
  cout << numIndices2 << " radii: (" << length1 << ", " << length2 << ", " << length3;
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
  planarMilestones.push_back(
      PlanarMilestoneInfo(numIndices1, numIndices2, length1, length2, length3, atomIndices1, atomIndices2, dataFileName));
  endOnMiddleCrossing = endOnMiddleCrossingArg;
  
  
}

void SeekrForce::modifyPlanarMilestone(int forceIndex, int numIndices1, int numIndices2, float length1, 
                              float length2, float length3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2, bool endOnMiddleCrossingArg,
                              std::string dataFileName) 
{
  /*int i;
  planarMilestone.numIndices1 = numIndices1;
  planarMilestone.numIndices2 = numIndices2;
  planarMilestone.length1 = length1;
  planarMilestone.length2 = length2;
  planarMilestone.length3 = length3;
  for (i=0; i<numIndices1; i++) {
    planarMilestone.atomIndices1[i] = atomIndices1[i];
  }
  for (i=0; i<numIndices2; i++) {
    planarMilestone.atomIndices2[i] = atomIndices2[i];
  }
  std::cout << "Modifying planar milestone. Indices: " << numIndices1 << ", " << numIndices2;
  std::cout << ". Radii: " << length1 << ", " << length2 << ", " << length3 << ".\n";
  */
  planarMilestones[forceIndex] = PlanarMilestoneInfo(numIndices1, numIndices2, length1, length2, length3, atomIndices1, atomIndices2, dataFileName);
  endOnMiddleCrossing = endOnMiddleCrossingArg;
  
}


ForceImpl* SeekrForce::createImpl() const {
    return new SeekrForceImpl(*this);
}

void SeekrForce::updateParametersInContext(Context& context) {
    dynamic_cast<SeekrForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
