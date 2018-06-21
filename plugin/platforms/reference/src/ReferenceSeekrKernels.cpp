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


#include "ReferenceSeekrKernels.h"
#include "SeekrForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"
#include <sstream>
#include <iostream>

using namespace SeekrPlugin;
using namespace OpenMM;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}
/*
static vector<RealVec>& extractMasses(ContextImpl& context) {
    vector<float>& masses;
    ReferencePlatform::PlatformData* system = reinterpret_cast<ReferencePlatform::PlatformData*>(contest.getSystem);
    
    return *((vector<float>*) masses);
} */

// Run some error collection routines

void checkAtomIndex(int numAtoms, const std::string& milestoneType, const int atomIndex) {
  bool error = false;
  if (atomIndex < 0) {
    error = true;
  }
  if (atomIndex >= numAtoms) {
    error = true;
  }
  
  if (error){
    std::stringstream m;
    m<<"Error: bad index for milestone of type: "<<milestoneType<<". atomIndex: "<<atomIndex;
    m<<"numAtoms: "<<numAtoms;
    throw OpenMMException(m.str());
  }
}

void checkSphericalMilestone(float radius1, float radius2, float radius3) {
  bool error = false;
  if (radius1 >= radius2) {
    error = true;
  }
  if (radius2 >= radius3) {
    error = true;
  }
  if (error){
    std::stringstream m;
    m<<"Error: bad radii for milestone of type: spherical milestone";
    m<<". radius1 (inner): "<<radius1<<" radius2 (middle): "<<radius2<< "radius3 (outer):"<<radius3;
    throw OpenMMException(m.str());
  }
}

void ReferenceCalcSeekrForceKernel::initialize(const System& system, const SeekrForce& force) {
    // Initialize bond parameters.
    int numAtoms = system.getNumParticles();
    std::string milestoneType = "spherical milestone";
    numSphericalMilestones = force.getNumSphericalMilestones();
    
    int currentAtomIndex1 = 0;
    int currentAtomIndex2 = 0;
    for (int i = 0; i < numSphericalMilestones; i++) {
      sphericalNumIndices1.push_back(force.getSphericalNumIndices(i, 1));
      sphericalNumIndices2.push_back(force.getSphericalNumIndices(i, 2));
      sphericalRadii1.push_back(force.getSphericalRadius(i, 1));
      sphericalRadii2.push_back(force.getSphericalRadius(i, 2));
      sphericalRadii3.push_back(force.getSphericalRadius(i, 3));
    
      checkSphericalMilestone(sphericalRadii1[i], sphericalRadii2[i], sphericalRadii3[i]); 
      sphericalAtomParamsBegin1.push_back(currentAtomIndex1);
      for(int j=0; j<sphericalNumIndices1[i]; j++) {
        int atom_id;
        force.getSphericalMilestoneAtoms(i, j, atom_id, 1);
        checkAtomIndex(numAtoms, milestoneType, atom_id);
        sphericalAtomIndices1.push_back(atom_id);
        currentAtomIndex1++;
      }
      sphericalAtomParamsEnd1.push_back(currentAtomIndex1);
      
      sphericalAtomParamsBegin2.push_back(currentAtomIndex2);
      for(int j=0; j<sphericalNumIndices2[i]; j++) {
        int atom_id;
        force.getSphericalMilestoneAtoms(i, j, atom_id, 2);
        checkAtomIndex(numAtoms, milestoneType, atom_id);
        sphericalAtomIndices2.push_back(atom_id);
        currentAtomIndex2++;
      }
      sphericalAtomParamsEnd2.push_back(currentAtomIndex2);
      crossed.push_back(0);
    }
}

double ReferenceCalcSeekrForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& pos = extractPositions(context);
    //vector<RealVec>& force = extractForces(context);
    //int numBonds = particle1.size();
    //double energy = 0;
    
    RealVec com1;
    RealVec com2;
    RealVec old_com1;
    RealVec old_com2;
    float totalmass1;
    float totalmass2;
    
    float mass;
    int i, j;
    int atomIndex;
    
    for (i=0; i < numSphericalMilestones; i++) {
        totalmass1 = 0.0;
        totalmass2 = 0.0;
        com1[0] = 0.0; com1[1] = 0.0; com1[2] = 0.0;
        com2[0] = 0.0; com2[1] = 0.0; com2[2] = 0.0;
        for (j=0; j<sphericalNumIndices1[i]; j++) {
            atomIndex = sphericalAtomIndices1[sphericalAtomParamsBegin1[i]+j];
            mass = 1.0; 
            com1 += pos[atomIndex]*mass;
            totalmass1 += mass;
        }
        com1 = com1 / totalmass1;
    
        for (j=0; j<sphericalNumIndices2[i]; j++) {
            atomIndex = sphericalAtomIndices2[sphericalAtomParamsBegin2[i]+j];
            mass = 1.0;
            com2 += pos[atomIndex]*mass;
            totalmass2 += mass;
        }
        com2 = com2 / totalmass2;
        
        RealVec delta = com2-com1;
        //RealVec old_delta = old_com2-old_com1;
        RealOpenMM r2 = delta.dot(delta);
        //RealOpenMM old_r2 = old_delta.dot(old_delta);
        RealOpenMM old_r2 = delta.dot(delta); // TODO: change this to old coordinates
        
        if (r2 < sphericalRadii1[i]*sphericalRadii1[i]) {
            crossed[i] = 1;
        } 
        else if ((r2 - sphericalRadii2[i]*sphericalRadii2[i])*(old_r2 - sphericalRadii2[i]*sphericalRadii2[i]) < 0.0) { // this will return true if the middle milestone is crossed
            crossed[i] = 2;
        }
        else if (r2 > sphericalRadii3[i]*sphericalRadii3[i]) {
            crossed[i] = 3;
        }
        
        // TODO: what to do with this calculation? 'crossed' is not put anywhere
    }
}

void ReferenceCalcSeekrForceKernel::copyParametersToContext(ContextImpl& context, const SeekrForce& force) {
    /*
    if (force.getNumBonds() != particle1.size())
        throw OpenMMException("updateParametersInContext: The number of Seekr bonds has changed");
    for (int i = 0; i < force.getSphericalNumBonds(); i++) {
        int p1, p2;
        force.getBondParameters(i, p1, p2, length[i], k[i]);
        if (p1 != particle1[i] || p2 != particle2[i])
            throw OpenMMException("updateParametersInContext: A particle index has changed");
    }
    */
}
