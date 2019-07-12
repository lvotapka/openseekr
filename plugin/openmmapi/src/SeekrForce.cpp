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
  if (molecule == 1) {
    atom_id = planarZMilestones[forceIndex].atomIndices1[atomIndex];
  } else if (molecule == 2) {
    atom_id = planarZMilestones[forceIndex].atomIndices2[atomIndex];
  } else {
    throw OpenMMException("Molecule index must be 1 or 2");
  }
}

int SeekrForce::getNumSphericalMilestones() const {
  return sphericalMilestones.size();
}

int SeekrForce::getNumSphericalAtomIndices() const {
  int total = 0;
  for(std::vector<SphericalMilestoneInfo>::const_iterator iter=sphericalMilestones.begin();
      iter != sphericalMilestones.end(); ++iter) {
  total += iter->numIndices1;
  total += iter->numIndices2;
  }
  return total;
}

int SeekrForce::getSphericalNumIndices(int forceIndex, int molecule) const {
  if (molecule == 1) {
    return sphericalMilestones[forceIndex].numIndices1;
  } else if (molecule == 2) {
    return sphericalMilestones[forceIndex].numIndices2;
  } else {
    throw OpenMMException("Molecule index must be 1 or 2");
  }
}

float SeekrForce::getSphericalRadius(int forceIndex, int milestone_id) const {
  if (milestone_id == 1) {
    return sphericalMilestones[forceIndex].radius1;
  } else if (milestone_id == 2) {
    return sphericalMilestones[forceIndex].radius2;
  } else if (milestone_id == 3) {
    return sphericalMilestones[forceIndex].radius3;
  } else {
    throw OpenMMException("Molecule index must be 1, 2, or 3");
  }
}

void SeekrForce::getSphericalMilestoneAtoms(int forceIndex, int atomIndex, int& atom_id, int molecule) const {
  if (molecule == 1) {
    atom_id = sphericalMilestones[forceIndex].atomIndices1[atomIndex];
  } else if (molecule == 2) {
    atom_id = sphericalMilestones[forceIndex].atomIndices2[atomIndex];
  } else {
    throw OpenMMException("Molecule index must be 1 or 2");
  }
}

int SeekrForce::getNumRmsdMilestones() const {
  return rmsdMilestones.size();
}

int SeekrForce::getNumRmsdAtomIndices() const {
  int total = 0;
  for(std::vector<RmsdMilestoneInfo>::const_iterator iter=rmsdMilestones.begin();
      iter != rmsdMilestones.end(); ++iter) {
  total += iter->numIndices1;
  total += iter->numIndices2;
  }
  return total;
}

int SeekrForce::getRmsdNumIndices(int forceIndex, int molecule) const {
  if (molecule == 1) {
    return rmsdMilestones[forceIndex].numIndices1;
  } else if (molecule == 2) {
    return rmsdMilestones[forceIndex].numIndices2;
  } else {
    throw OpenMMException("Molecule index must be 1 or 2");
  }
}

float SeekrForce::getRmsdRadius(int forceIndex, int milestone_id) const {
  if (milestone_id == 1) {
    return rmsdMilestones[forceIndex].radius1;
  } else if (milestone_id == 2) {
    return rmsdMilestones[forceIndex].radius2;
  } else if (milestone_id == 3) {
    return rmsdMilestones[forceIndex].radius3;
  } else {
    throw OpenMMException("Molecule index must be 1, 2, or 3");
  }
}

void SeekrForce::getRmsdMilestoneAtoms(int forceIndex, int atomIndex, int& atom_id, int molecule) const {
  if (molecule == 1) {
    atom_id = rmsdMilestones[forceIndex].atomIndices1[atomIndex];
  } else if (molecule == 2) {
    atom_id = rmsdMilestones[forceIndex].atomIndices2[atomIndex];
  } else {
    throw OpenMMException("Molecule index must be 1 or 2");
  }
}

void SeekrForce::addPlanarZMilestone(int numIndices1, int numIndices2, float offset1, 
                              float offset2, float offset3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2) 
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
      PlanarZMilestoneInfo(numIndices1, numIndices2, offset1, offset2, offset3, atomIndices1, atomIndices2));
  //endOnMiddleCrossing = endOnMiddleCrossingArg;
  
  
}

void SeekrForce::modifyPlanarZMilestone(int forceIndex, int numIndices1, int numIndices2, float offset1, 
                              float offset2, float offset3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2) 
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
  planarZMilestones[forceIndex] = PlanarZMilestoneInfo(numIndices1, numIndices2, offset1, offset2, offset3, atomIndices1, atomIndices2);
  //endOnMiddleCrossing = endOnMiddleCrossingArg;
  
}


void SeekrForce::addSphericalMilestone(int numIndices1, int numIndices2, float radius1, 
                              float radius2, float radius3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2) 
{
  sphericalMilestones.push_back(
      SphericalMilestoneInfo(numIndices1, numIndices2, radius1, radius2, radius3, atomIndices1, atomIndices2));
  //endOnMiddleCrossing = endOnMiddleCrossingArg;
  
  
}

void SeekrForce::modifySphericalMilestone(int forceIndex, int numIndices1, int numIndices2, float radius1, 
                              float radius2, float radius3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2) 
{
  sphericalMilestones[forceIndex] = SphericalMilestoneInfo(numIndices1, numIndices2, radius1, radius2, radius3, atomIndices1, atomIndices2);
  //endOnMiddleCrossing = endOnMiddleCrossingArg;
  
}

void SeekrForce::addRmsdMilestone(int numIndices1, int numIndices2, float radius1, 
                              float radius2, float radius3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2) 
{
  rmsdMilestones.push_back(
      RmsdMilestoneInfo(numIndices1, numIndices2, radius1, radius2, radius3, atomIndices1, atomIndices2));
  //endOnMiddleCrossing = endOnMiddleCrossingArg;
  
  
}

void SeekrForce::modifyRmsdMilestone(int forceIndex, int numIndices1, int numIndices2, float radius1, 
                              float radius2, float radius3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2) 
{
  rmsdMilestones[forceIndex] = RmsdMilestoneInfo(numIndices1, numIndices2, radius1, radius2, radius3, atomIndices1, atomIndices2);
  //endOnMiddleCrossing = endOnMiddleCrossingArg;
  
}

bool SeekrForce::getEndOnMiddleCrossing() const {
  return endOnMiddleCrossing;
}

void SeekrForce::setEndOnMiddleCrossing(bool endOnMiddleCrossingArg) {
  endOnMiddleCrossing = endOnMiddleCrossingArg;
}

std::string SeekrForce::getDataFileName() const {
  return dataFileName;
}

void SeekrForce::setDataFileName(std::string dataFileNameArg) {
  dataFileName = dataFileNameArg;
}

std::string SeekrForce::getSaveStateFileName() const {
  return stateFileName;
}

void SeekrForce::setSaveStateFileName(std::string saveStateFileName) {
  stateFileName = saveStateFileName;
}

ForceImpl* SeekrForce::createImpl() const {
    return new SeekrForceImpl(*this);
}

void SeekrForce::updateParametersInContext(Context& context) {
    dynamic_cast<SeekrForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
