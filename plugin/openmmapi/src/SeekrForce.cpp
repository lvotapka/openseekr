/*
   Copyright 2018 by Lane Votapka
   All rights reserved
   
   -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

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

bool SeekrForce::getEndOnMiddleCrossing() const {
  return endOnMiddleCrossing;
}

std::string SeekrForce::getDataFileName(int forceIndex) const {
  return sphericalMilestones[forceIndex].dataFileName;
}

const string& SeekrForce::getSaveStateFileName() const {
    return saveStateFileName;
}

void SeekrForce::setSaveStateFileName(const string& fileName) {
    saveStateFileName = fileName;
}

void SeekrForce::addSphericalMilestone(int numIndices1, int numIndices2, float radius1, 
                              float radius2, float radius3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2, bool endOnMiddleCrossingArg,
                              std::string dataFileName) 
{
  /*cout << "Adding spherical milestone. numIndices1: " << numIndices1 << " numIncices: ";
  cout << numIndices2 << " radii: (" << radius1 << ", " << radius2 << ", " << radius3;
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
  sphericalMilestones.push_back(
      SphericalMilestoneInfo(numIndices1, numIndices2, radius1, radius2, radius3, atomIndices1, atomIndices2, dataFileName));
  endOnMiddleCrossing = endOnMiddleCrossingArg;
  
  
}

void SeekrForce::modifySphericalMilestone(int forceIndex, int numIndices1, int numIndices2, float radius1, 
                              float radius2, float radius3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2, bool endOnMiddleCrossingArg,
                              std::string dataFileName) 
{
  /*int i;
  sphericalMilestone.numIndices1 = numIndices1;
  sphericalMilestone.numIndices2 = numIndices2;
  sphericalMilestone.radius1 = radius1;
  sphericalMilestone.radius2 = radius2;
  sphericalMilestone.radius3 = radius3;
  for (i=0; i<numIndices1; i++) {
    sphericalMilestone.atomIndices1[i] = atomIndices1[i];
  }
  for (i=0; i<numIndices2; i++) {
    sphericalMilestone.atomIndices2[i] = atomIndices2[i];
  }
  std::cout << "Modifying spherical milestone. Indices: " << numIndices1 << ", " << numIndices2;
  std::cout << ". Radii: " << radius1 << ", " << radius2 << ", " << radius3 << ".\n";
  */
  sphericalMilestones[forceIndex] = SphericalMilestoneInfo(numIndices1, numIndices2, radius1, radius2, radius3, atomIndices1, atomIndices2, dataFileName);
  endOnMiddleCrossing = endOnMiddleCrossingArg;
  
}


ForceImpl* SeekrForce::createImpl() const {
    return new SeekrForceImpl(*this);
}

void SeekrForce::updateParametersInContext(Context& context) {
    dynamic_cast<SeekrForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
