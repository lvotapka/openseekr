/*
   Copyright 2018 by Lane Votapka
   All rights reserved
 * -------------------------------------------------------------------------- *
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

#include "CudaSeekrKernels.h"
#include "CudaSeekrKernelSources.h"
#include "openmm/cuda/CudaContext.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/cuda/CudaBondedUtilities.h"
#include "openmm/cuda/CudaForceInfo.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>

using namespace SeekrPlugin;
using namespace OpenMM;
using namespace std;

#define CHECK_RESULT(result) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<errorMessage<<": "<<cu.getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

class CudaSeekrForceInfo : public CudaForceInfo {
public:
    CudaSeekrForceInfo(const SeekrForce& force) : force(force) {
    }
    
    
    int getNumParticleGroups() {
        //return force.getNumBonds();
        return 0;
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        
    }
    bool areGroupsIdentical(int group1, int group2) {
        return false;
    } 
private:
    const SeekrForce& force;
};

CudaCalcSeekrForceKernel::CudaCalcSeekrForceKernel(std::string name, const Platform& platform, CudaContext& cu,
                                                   const System& system) :
    //CalcSeekrForceKernel(name, platform), hasInitializedKernel(false), cu(cu), system(system), params(NULL)
    CalcSeekrForceKernel(name, platform), cu(cu), system(system)
{
    if (cu.getUseDoublePrecision()) {
        cout << "***\n";
        cout << "*** SeekrForce does not support double precision.\n";
        cout << "***" << endl;
        throw OpenMMException("SeekrForce does not support double precision");
    }
    numSphericalMilestones = 0;
    numSphericalAtomIndices = 0;
    
    sphericalNumIndices1 = nullptr;
    sphericalNumIndices2 = nullptr;
    sphericalRadii1 = nullptr;
    sphericalRadii2 = nullptr;
    sphericalRadii3 = nullptr;
    sphericalAtomIndices1 = nullptr;
    sphericalAtomIndices2 = nullptr;
    sphericalAtomBounds1 = nullptr;
    sphericalAtomBounds2 = nullptr;
    sphericalOldCom1 = nullptr;
    sphericalOldCom2 = nullptr;
    
    collectionReturnCode = nullptr;
    
    endSimulation = false;
    endOnMiddleCrossing = false;
    crossedStartingMilestone = false;
}

CudaCalcSeekrForceKernel::~CudaCalcSeekrForceKernel() {
    cu.setAsCurrent();
    
    delete sphericalNumIndices1;
    delete sphericalNumIndices2;
    delete sphericalRadii1;
    delete sphericalRadii2;
    delete sphericalRadii3;
    delete sphericalAtomIndices1;
    delete sphericalAtomIndices2;
    delete sphericalAtomBounds1;
    delete sphericalAtomBounds2;
    delete sphericalOldCom1;
    delete sphericalOldCom2;
    delete collectionReturnCode;
    //if (params != NULL)
    //    delete params;
}

void CudaCalcSeekrForceKernel::allocateMemory(const SeekrForce& force) {
  
  numSphericalMilestones = force.getNumSphericalMilestones();
  numSphericalAtomIndices = force.getNumSphericalAtomIndices();
  
  sphericalNumIndices1 =  CudaArray::create<int> (cu, numSphericalMilestones, "sphericalNumIndices1");
  sphericalNumIndices2 =  CudaArray::create<int> (cu, numSphericalMilestones, "sphericalNumIndices2");
  sphericalRadii1 =       CudaArray::create<float> (cu, numSphericalMilestones, "sphericalRadii1");
  sphericalRadii2 =       CudaArray::create<float> (cu, numSphericalMilestones, "sphericalRadii2");
  sphericalRadii3 =       CudaArray::create<float> (cu, numSphericalMilestones, "sphericalRadii3");
  sphericalAtomIndices1 = CudaArray::create<int> (cu, numSphericalAtomIndices, "sphericalAtomIndices1");
  sphericalAtomIndices2 = CudaArray::create<int> (cu, numSphericalAtomIndices, "sphericalAtomIndices2");
  sphericalAtomBounds1  = CudaArray::create<int2> (cu, numSphericalMilestones, "sphericalAtomBounds1");
  sphericalAtomBounds2  = CudaArray::create<int2> (cu, numSphericalMilestones, "sphericalAtomBounds2");
  sphericalOldCom1      = CudaArray::create<float4> (cu, numSphericalMilestones, "sphericalOldCom1");
  sphericalOldCom2      = CudaArray::create<float4> (cu, numSphericalMilestones, "sphericalOldCom2");
  collectionReturnCode  =  CudaArray::create<float> (cu, numSphericalMilestones, "collectionReturnCode");
  //collectionReturnCode =  CudaArray::create<int> (cu, numSphericalMilestones, "collectionReturnCode");
  
  // host memory
  
  h_sphericalNumIndices1  = std::vector<int> (numSphericalMilestones, 0);
  h_sphericalNumIndices2  = std::vector<int> (numSphericalMilestones, 0);
  h_sphericalRadii1       = std::vector<float> (numSphericalMilestones, 0.0);
  h_sphericalRadii2       = std::vector<float> (numSphericalMilestones, 0.0);
  h_sphericalRadii3       = std::vector<float> (numSphericalMilestones, 0.0);
  h_sphericalAtomIndices1 = std::vector<int> (numSphericalAtomIndices, 0);
  h_sphericalAtomIndices2 = std::vector<int> (numSphericalAtomIndices, 0);
  h_sphericalAtomBounds1  = std::vector<int2> (numSphericalMilestones, make_int2(-1, -1));
  h_sphericalAtomBounds2  = std::vector<int2> (numSphericalMilestones, make_int2(-1, -1));
  h_sphericalOldCom1      = std::vector<float4> (numSphericalMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
  h_sphericalOldCom2      = std::vector<float4> (numSphericalMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
  
  h_collectionReturnCode  = std::vector<float> (numSphericalMilestones, 0);
  //h_collectionReturnCode  = std::vector<int> (1, 0);
}
  //endSimulation = false;
  //endOnMiddleCrossing = force.getEndOnMiddleCrossing();
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
    m<<". radius1 (inner): "<<radius1<<" radius2 (middle): "<<radius2<<" radius3 (outer): "<<radius3;
    throw OpenMMException(m.str());
  }
}

void CudaCalcSeekrForceKernel::setupSphericalMilestones(const SeekrForce& force){
  int numAtoms = system.getNumParticles();
  std::string milestoneType = "spherical milestone";
  int currentAtomIndex1 = 0;
  int currentAtomIndex2 = 0;
  endOnMiddleCrossing = force.getEndOnMiddleCrossing();
  for (int i=0; i < numSphericalMilestones; i++) {
    
    int thisStart1 = currentAtomIndex1;
    int thisStart2 = currentAtomIndex2;
    // retrieve the spherical milestone parameters from the C++ api
    h_sphericalNumIndices1[i] = force.getSphericalNumIndices(i, 1); 
    h_sphericalNumIndices2[i] = force.getSphericalNumIndices(i, 2);
    h_sphericalRadii1[i] = force.getSphericalRadius(i, 1);
    h_sphericalRadii2[i] = force.getSphericalRadius(i, 2);
    h_sphericalRadii3[i] = force.getSphericalRadius(i, 3);
    
    // check the radii to make sure that they make sense
    checkSphericalMilestone(h_sphericalRadii1[i], h_sphericalRadii2[i], h_sphericalRadii3[i]); 
    
    // Retrieve the individual atoms from the API
    for(int j=0; j<h_sphericalNumIndices1[i]; j++) {
      int atom_id;
      force.getSphericalMilestoneAtoms(i, j, atom_id, 1);
      checkAtomIndex(numAtoms, milestoneType, atom_id);
      h_sphericalAtomIndices1[currentAtomIndex1] = atom_id;
      currentAtomIndex1++;
    }
    int thisEnd1 = currentAtomIndex1;
    h_sphericalAtomBounds1[i] = make_int2(thisStart1, thisEnd1);
    
    for(int j=0; j<h_sphericalNumIndices2[i]; j++) {
      int atom_id;
      force.getSphericalMilestoneAtoms(i, j, atom_id, 2);
      checkAtomIndex(numAtoms, milestoneType, atom_id);
      h_sphericalAtomIndices2[currentAtomIndex2] = atom_id;
      currentAtomIndex2++;
    }
    int thisEnd2 = currentAtomIndex2;
    h_sphericalAtomBounds2[i] = make_int2(thisStart2, thisEnd2);
    
    /*
    cout << "CUDA kernel: Setting up spherical milestone number: " << i << ". sphericalNumIndices1: ";
    cout << h_sphericalNumIndices1[i] << " sphericalNumIndices2: " << h_sphericalNumIndices2[i];
    cout << " sphericalRadii: (" << h_sphericalRadii1[i] << ", " << h_sphericalRadii2[i];
    cout << ", " << h_sphericalRadii3[i] << "). sphericalAtomBounds1: (";
    cout << h_sphericalAtomBounds1[i].x << ", " << h_sphericalAtomBounds1[i].y;
    cout << ") sphericalAtomIndices1: [";
    for (int j=h_sphericalAtomBounds1[i].x; j<h_sphericalAtomBounds1[i].y; j++) {
      cout << h_sphericalAtomIndices1[j] << " ";
    }
    cout << "] sphericalAtomBounds2: (";
    cout << h_sphericalAtomBounds2[i].x << ", " << h_sphericalAtomBounds2[i].y;
    cout << ") sphericalAtomIndices2: [";
    for (int j=h_sphericalAtomBounds2[i].x; j<h_sphericalAtomBounds2[i].y; j++) {
      cout << h_sphericalAtomIndices2[j] << " ";
    }
    cout << "]\n";*/
    dataFileNames.push_back(force.getDataFileName(i));
  }
}

void CudaCalcSeekrForceKernel::validateAndUpload() {
  cout << "Attempting to upload host arrays to device.\n";
  sphericalNumIndices1->upload(h_sphericalNumIndices1);
  sphericalNumIndices2->upload(h_sphericalNumIndices2);
  sphericalRadii1->upload(h_sphericalRadii1);
  sphericalRadii2->upload(h_sphericalRadii2);
  sphericalRadii3->upload(h_sphericalRadii3);
  sphericalAtomIndices1->upload(h_sphericalAtomIndices1);
  sphericalAtomIndices2->upload(h_sphericalAtomIndices2);
  sphericalAtomBounds1->upload(h_sphericalAtomBounds1);
  sphericalAtomBounds2->upload(h_sphericalAtomBounds2);
  sphericalOldCom1->upload(h_sphericalOldCom1);
  sphericalOldCom2->upload(h_sphericalOldCom2);
  cout << "Uploaded all host arrays to device.\n";
}

void CudaCalcSeekrForceKernel::initialize(const System& system, const SeekrForce& force) {
    cu.setAsCurrent();
    
    allocateMemory(force);
    setupSphericalMilestones(force);
    validateAndUpload();
    
    CUmodule module = cu.createModule(CudaSeekrKernelSources::vectorOps + CudaSeekrKernelSources::seekrForce);
    computeSphericalMilestonesKernel = cu.getKernel(module, "monitorSphericalMilestones");
    cu.addForce(new CudaSeekrForceInfo(force));
}

double CudaCalcSeekrForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (context.getTime() == 0.0) {
      endSimulation = false;
      h_sphericalOldCom1      = std::vector<float4> (numSphericalMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
      h_sphericalOldCom2      = std::vector<float4> (numSphericalMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
      sphericalOldCom1->upload(h_sphericalOldCom1);
      sphericalOldCom2->upload(h_sphericalOldCom2);
    }
    
    void* sphericalMilestoneArgs[] = {
      &cu.getPosq().getDevicePointer(),
      &cu.getVelm().getDevicePointer(),
      &sphericalNumIndices1->getDevicePointer(),
      &sphericalNumIndices2->getDevicePointer(),
      &sphericalRadii1->getDevicePointer(),
      &sphericalRadii2->getDevicePointer(),
      &sphericalRadii3->getDevicePointer(),
      &sphericalAtomIndices1->getDevicePointer(),
      &sphericalAtomIndices2->getDevicePointer(),
      &sphericalAtomBounds1->getDevicePointer(),
      &sphericalAtomBounds2->getDevicePointer(),
      &collectionReturnCode->getDevicePointer(),
      &sphericalOldCom1->getDevicePointer(),
      &sphericalOldCom2->getDevicePointer(),
      &numSphericalMilestones
    };
    cu.executeKernel(computeSphericalMilestonesKernel, sphericalMilestoneArgs, numSphericalMilestones);
    collectionReturnCode->download(h_collectionReturnCode); 
    //cout << "collectionReturnCode: " << h_collectionReturnCode[0] << endl;
    
    //cout << "step:" << context.getTime() << "\n";
    
    
    for (int i=0; i<numSphericalMilestones; i++) {
      if (h_collectionReturnCode[0] == 1) {
        if (endSimulation == false) {
          cout<<"Inner milestone crossed. Time:" << context.getTime() << " ps\n";
          endSimulation = true;
          ofstream datafile;
          datafile.open(dataFileNames[i], std::ios_base::app);
          if (crossedStartingMilestone == true) {
            datafile << "1 " << context.getTime() << "\n";
          } else {
            datafile << "1* " << context.getTime() << "\n";
          }
          crossedStartingMilestone = false;
          datafile.close();
        }
      } else if (h_collectionReturnCode[0] == 2) {
        if (endSimulation == false) {
          if (endOnMiddleCrossing == true) {
            cout<<"Middle milestone crossed. Time:" << context.getTime() << " ps\n";
            endSimulation = true;
            ofstream datafile;
            datafile.open(dataFileNames[i], std::ios_base::app);
            datafile << "2 " << context.getTime() << "\n";
            datafile.close();
          } else { // Then its the forward stage, so assert that this milestone is crossed
            if (crossedStartingMilestone == false) {
              crossedStartingMilestone = true;
              context.setTime(0.0);
            }
          }
        }
      } else if (h_collectionReturnCode[0] == 3) {
        if (endSimulation == false) {
          cout<<"Outer milestone crossed. Time:" << context.getTime() << " ps\n";
          endSimulation = true;
          ofstream datafile;
          datafile.open(dataFileNames[i], std::ios_base::app);
          if (crossedStartingMilestone == true) {
            datafile << "3 " << context.getTime() << "\n";
          } else {
            datafile << "3* " << context.getTime() << "\n";
          }
          crossedStartingMilestone = false;
          datafile.close();
        }
      }
    }
    return 0.0;
}

void CudaCalcSeekrForceKernel::copyParametersToContext(ContextImpl& context, const SeekrForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    //setupSphericalMilestones(force);
    //validateAndUpload();
    cout << "copyParametersToContext\n";
    
    cu.invalidateMolecules();
}
