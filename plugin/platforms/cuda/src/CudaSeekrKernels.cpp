/*
   Copyright 2019 by Andy Stokely
	All rights reserved
   PlanarZ Milestone C++ Code for the OpenMM Plugin SEEKR that links the seekrForce.cu CUDA code to 
   rest of the program


//path = /plugin/platforms/cuda/src/kernels
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

//TODO I just updated the "spherical", "Offset", etc names to their respective planarZ versions. I may have to go back and update the data structure types. (6-2-2019)

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
    numPlanarZMilestones = 0;
    numPlanarZAtomIndices = 0;
    
    planarZNumIndices1 = nullptr;
    planarZNumIndices2 = nullptr;
    planarZOffsets1 = nullptr;
    planarZOffsets2 = nullptr;
    planarZOffsets3 = nullptr;
    planarZAtomIndices1 = nullptr;
    planarZAtomIndices2 = nullptr; //TODO: ANDY-->RESOLVED
    planarZAtomBounds1 = nullptr;
    planarZAtomBounds2 = nullptr;
    planarZOldCom1 = nullptr;
    planarZOldCom2 = nullptr;
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
    planarZCollectionReturnCode = nullptr;
    
    endSimulation = false;
    endOnMiddleCrossing = false;
    crossedStartingMilestone = false;
}

CudaCalcSeekrForceKernel::~CudaCalcSeekrForceKernel() {
    cu.setAsCurrent();
    
    delete planarZNumIndices1;
    delete planarZNumIndices2;
    delete planarZOffsets1;
    delete planarZOffsets2;
    delete planarZOffsets3;
    delete planarZAtomIndices1;
    delete planarZAtomIndices2; //TODO: ANDY-->RESOLVED
    delete planarZAtomBounds1;
    delete planarZAtomBounds2;
    delete planarZOldCom1;
    delete planarZOldCom2;
    delete planarZCollectionReturnCode;
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
  
  numPlanarZMilestones = force.getNumPlanarZMilestones();
  numPlanarZAtomIndices = force.getNumPlanarZAtomIndices();
  
  planarZNumIndices1 =  CudaArray::create<int> (cu, numPlanarZMilestones, "planarZNumIndices1");
  planarZNumIndices2 =  CudaArray::create<int> (cu, numPlanarZMilestones, "planarZNumIndices2");
  planarZOffsets1 =       CudaArray::create<float> (cu, numPlanarZMilestones, "planarZOffsets1");
  planarZOffsets2 =       CudaArray::create<float> (cu, numPlanarZMilestones, "planarZOffsets2");
  planarZOffsets3 =       CudaArray::create<float> (cu, numPlanarZMilestones, "planarZOffsets3");
  planarZAtomIndices1 = CudaArray::create<int> (cu, numPlanarZAtomIndices, "planarZAtomIndices1");
  planarZAtomIndices2 = CudaArray::create<int> (cu, numPlanarZAtomIndices, "planarZAtomIndices2"); // TODO: ANDY
  planarZAtomBounds1  = CudaArray::create<int2> (cu, numPlanarZMilestones, "planarZAtomBounds1");
  planarZAtomBounds2  = CudaArray::create<int2> (cu, numPlanarZMilestones, "planarZAtomBounds2");
  planarZOldCom1      = CudaArray::create<float4> (cu, numPlanarZMilestones, "planarZOldCom1");
  planarZOldCom2      = CudaArray::create<float4> (cu, numPlanarZMilestones, "planarZOldCom2");
  planarZCollectionReturnCode  =  CudaArray::create<float> (cu, numPlanarZMilestones, "planarZCollectionReturnCode"); // TODO: ANDY watch out
  //planarZCollectionReturnCode =  CudaArray::create<int> (cu, numPlanarZMilestones, "planarZCollectionReturnCode");
  
  // host memory (CPU)
  
  h_planarZNumIndices1  = std::vector<int> (numPlanarZMilestones, 0);
  h_planarZNumIndices2  = std::vector<int> (numPlanarZMilestones, 0);
  h_planarZOffsets1       = std::vector<float> (numPlanarZMilestones, 0.0);
  h_planarZOffsets2       = std::vector<float> (numPlanarZMilestones, 0.0);
  h_planarZOffsets3       = std::vector<float> (numPlanarZMilestones, 0.0);
  h_planarZAtomIndices1 = std::vector<int> (numPlanarZAtomIndices, 0);
  h_planarZAtomIndices2 = std::vector<int> (numPlanarZAtomIndices, 0);
  h_planarZAtomBounds1  = std::vector<int2> (numPlanarZMilestones, make_int2(-1, -1));
  h_planarZAtomBounds2  = std::vector<int2> (numPlanarZMilestones, make_int2(-1, -1));
  h_planarZOldCom1      = std::vector<float4> (numPlanarZMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
  h_planarZOldCom2      = std::vector<float4> (numPlanarZMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
  
  h_planarZCollectionReturnCode  = std::vector<float> (numPlanarZMilestones, 0);
  //h_planarZCollectionReturnCode  = std::vector<int> (1, 0);


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

void checkPlanarZMilestone(float Offsets1, float Offsets2, float Offsets3) { // TODO: ANDY make one for planarZ 
  bool error = false;
  if (Offsets1 >= Offsets2) {
    error = true;
  }
  if (Offsets2 >= Offsets3) {
    error = true;
  }
  
  if (error){
    std::stringstream m;
    m<<"Error: bad Offsetss for milestone of type: planarZ milestone";
    m<<". Offsets1 (inner): "<<Offsets1<<" Offsets2 (middle): "<<Offsets2<<" Offsets3 (outer): "<<Offsets3;
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


void CudaCalcSeekrForceKernel::setupPlanarZMilestones(const SeekrForce& force){ //TODO: ANDY create your own DONE
  int numAtoms = system.getNumParticles();
  std::string milestoneType = "planarZ milestone";
  int currentAtomIndex1 = 0;
  int currentAtomIndex2 = 0;
  endOnMiddleCrossing = force.getEndOnMiddleCrossing();
  for (int i=0; i < numPlanarZMilestones; i++) {
    
    int thisStart1 = currentAtomIndex1;
    int thisStart2 = currentAtomIndex2;
    // retrieve the planarZ milestone parameters from the C++ api
    h_planarZNumIndices1[i] = force.getPlanarZNumIndices(i, 1); 
    h_planarZNumIndices2[i] = force.getPlanarZNumIndices(i, 2);
    h_planarZOffsets1[i] = force.getPlanarZOffset(i, 1);
    h_planarZOffsets2[i] = force.getPlanarZOffset(i, 2);
    h_planarZOffsets3[i] = force.getPlanarZOffset(i, 3);
    
    // check the Offsetss to make sure that they make sense
    checkPlanarZMilestone(h_planarZOffsets1[i], h_planarZOffsets2[i], h_planarZOffsets3[i]); 
    
    // Retrieve the individual atoms from the API
    for(int j=0; j<h_planarZNumIndices1[i]; j++) {
      int atom_id;
      force.getPlanarZMilestoneAtoms(i, j, atom_id, 1);
      checkAtomIndex(numAtoms, milestoneType, atom_id);
      h_planarZAtomIndices1[currentAtomIndex1] = atom_id;
      currentAtomIndex1++;
    }
    int thisEnd1 = currentAtomIndex1;
    h_planarZAtomBounds1[i] = make_int2(thisStart1, thisEnd1);
    
    for(int j=0; j<h_planarZNumIndices2[i]; j++) {
      int atom_id;
      force.getPlanarZMilestoneAtoms(i, j, atom_id, 2);
      checkAtomIndex(numAtoms, milestoneType, atom_id);
      h_planarZAtomIndices2[currentAtomIndex2] = atom_id;
      currentAtomIndex2++;
    }
    int thisEnd2 = currentAtomIndex2;
    h_planarZAtomBounds2[i] = make_int2(thisStart2, thisEnd2);
    
    /*
    cout << "CUDA kernel: Setting up planarZ milestone number: " << i << ". planarZNumIndices1: ";
    cout << h_planarZNumIndices1[i] << " planarZNumIndices2: " << h_planarZNumIndices2[i];
    cout << " planarZOffsets: (" << h_planarZOffsets1[i] << ", " << h_planarZOffsets2[i];
    cout << ", " << h_planarZOffsets3[i] << "). planarZAtomBounds1: (";
    cout << h_planarZAtomBounds1[i].x << ", " << h_planarZAtomBounds1[i].y;
    cout << ") planarZAtomIndices1: [";
    for (int j=h_planarZAtomBounds1[i].x; j<h_planarZAtomBounds1[i].y; j++) {
      cout << h_planarZAtomIndices1[j] << " ";
    }
    cout << "] planarZAtomBounds2: (";
    cout << h_planarZAtomBounds2[i].x << ", " << h_planarZAtomBounds2[i].y;
    cout << ") planarZAtomIndices2: [";
    for (int j=h_planarZAtomBounds2[i].x; j<h_planarZAtomBounds2[i].y; j++) {
      cout << h_planarZAtomIndices2[j] << " ";
    }
    cout << "]\n";*/
    dataFileNames.push_back(force.getSphericalDataFileName(i));
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
    dataFileNames.push_back(force.getPlanarZDataFileName(i));
  }
}

void CudaCalcSeekrForceKernel::validateAndUpload() {
  cout << "Attempting to upload host arrays to device.\n";
  planarZNumIndices1->upload(h_planarZNumIndices1);
  planarZNumIndices2->upload(h_planarZNumIndices2);
  planarZOffsets1->upload(h_planarZOffsets1);
  planarZOffsets2->upload(h_planarZOffsets2);
  planarZOffsets3->upload(h_planarZOffsets3);
  planarZAtomIndices1->upload(h_planarZAtomIndices1);
  planarZAtomIndices2->upload(h_planarZAtomIndices2); // TODO: ANDY-->RESOLVED
  planarZAtomBounds1->upload(h_planarZAtomBounds1);
  planarZAtomBounds2->upload(h_planarZAtomBounds2);
  planarZOldCom1->upload(h_planarZOldCom1);
  planarZOldCom2->upload(h_planarZOldCom2);
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
    setupPlanarZMilestones(force); // TODO:ANDY add yours here--RESOLVED
    setupSphericalMilestones(force);
    validateAndUpload();
    
    CUmodule module = cu.createModule(CudaSeekrKernelSources::vectorOps + CudaSeekrKernelSources::seekrForce);
    computePlanarZMilestonesKernel = cu.getKernel(module, "monitorPlanarZMilestones");
    computeSphericalMilestonesKernel = cu.getKernel(module, "monitorSphericalMilestones");
    cu.addForce(new CudaSeekrForceInfo(force));
}

double CudaCalcSeekrForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (context.getTime() == 0.0) {
      endSimulation = false;
      h_planarZOldCom1      = std::vector<float4> (numPlanarZMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
      h_planarZOldCom2      = std::vector<float4> (numPlanarZMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
      planarZOldCom1->upload(h_planarZOldCom1);
      planarZOldCom2->upload(h_planarZOldCom2);
       h_sphericalOldCom1      = std::vector<float4> (numSphericalMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
      h_sphericalOldCom2      = std::vector<float4> (numSphericalMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
      sphericalOldCom1->upload(h_sphericalOldCom1);
      sphericalOldCom2->upload(h_sphericalOldCom2);
    }
    
    void* planarZMilestoneArgs[] = {
      &cu.getPosq().getDevicePointer(),
      &cu.getVelm().getDevicePointer(),
      &planarZNumIndices1->getDevicePointer(),
      &planarZNumIndices2->getDevicePointer(),
      &planarZOffsets1->getDevicePointer(),
      &planarZOffsets2->getDevicePointer(),
      &planarZOffsets3->getDevicePointer(),
      &planarZAtomIndices1->getDevicePointer(),
      &planarZAtomIndices2->getDevicePointer(),
      &planarZAtomBounds1->getDevicePointer(),
      &planarZAtomBounds2->getDevicePointer(),
      &planarZCollectionReturnCode->getDevicePointer(),
      &planarZOldCom1->getDevicePointer(),
      &planarZOldCom2->getDevicePointer(),
      &numPlanarZMilestones
    };
    cu.executeKernel(computePlanarZMilestonesKernel, planarZMilestoneArgs, numPlanarZMilestones);
    planarZCollectionReturnCode->download(h_planarZCollectionReturnCode); 
    //cout << "planarZCollectionReturnCode: " << h_planarZCollectionReturnCode[0] << endl;
    
    //cout << "step:" << context.getTime() << "\n";
 
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
    
    
    for (int i=0; i<numPlanarZMilestones; i++) {
      if (h_planarZCollectionReturnCode[i] == 1) {
        if (endSimulation == false) { // if we haven't already crossed an ending milestone
          cout<<"Inner milestone crossed. Time:" << context.getTime() << " ps\n"; // output info to user
          endSimulation = true; // make sure we end after this point
          ofstream datafile; // open datafile for writing
          datafile.open(dataFileNames[i], std::ios_base::app); // append to file
          if (crossedStartingMilestone == true || endOnMiddleCrossing == true) { // if it's the reversal stage or we've crossed the starting milestone
            datafile << "1 " << context.getTime() << "\n";
          } else { // if its the forward stage and the starting milestone isn't crossed
            datafile << "1* " << context.getTime() << "\n";
          }
          crossedStartingMilestone = false; // reset whether we've crossed the starting milestone for the next simulation
          datafile.close(); // close data file
        }
      } else if (h_planarZCollectionReturnCode[i] == 2) {
        if (endSimulation == false) { // if we haven't already crossed an ending milestone
          if (endOnMiddleCrossing == true) { // then it's a reversal stage
            cout<<"Middle milestone crossed. Time:" << context.getTime() << " ps\n";
            endSimulation = true;
            ofstream datafile;
            datafile.open(dataFileNames[i], std::ios_base::app);
            datafile << "2 " << context.getTime() << "\n";
            datafile.close();
          } else { // Then its the forward stage, so assert that this milestone is crossed
            if (crossedStartingMilestone == false) { // as long as we haven't crossed once already
              crossedStartingMilestone = true; // we've crossed it once
              context.setTime(0.0); // reset the timer
            }
          }
        }
      } else if (h_planarZCollectionReturnCode[i] == 3) { // This should be a fairly identical procedure to inner milestone crossing above
        if (endSimulation == false) {
          cout<<"Outer milestone crossed. Time:" << context.getTime() << " ps\n";
          endSimulation = true;
          ofstream datafile;
          datafile.open(dataFileNames[i], std::ios_base::app);
          if (crossedStartingMilestone == true || endOnMiddleCrossing == true) {
            datafile << "3 " << context.getTime() << "\n";
          } else {
            datafile << "3* " << context.getTime() << "\n";
          }
          crossedStartingMilestone = false;
          datafile.close();
        }
      }
    }
    
    for (int i=0; i<numSphericalMilestones; i++) {
      if (h_collectionReturnCode[0] == 1) {
        if (endSimulation == false) { // if we haven't already crossed an ending milestone
          cout<<"Inner milestone crossed. Time:" << context.getTime() << " ps\n"; // output info to user
          endSimulation = true; // make sure we end after this point
          ofstream datafile; // open datafile for writing
          datafile.open(dataFileNames[i], std::ios_base::app); // append to file
          if (crossedStartingMilestone == true || endOnMiddleCrossing == true) { // if it's the reversal stage or we've crossed the starting milestone
            datafile << "1 " << context.getTime() << "\n";
          } else { // if its the forward stage and the starting milestone isn't crossed
            datafile << "1* " << context.getTime() << "\n";
          }
          crossedStartingMilestone = false; // reset whether we've crossed the starting milestone for the next simulation
          datafile.close(); // close data file
        }
      } else if (h_collectionReturnCode[0] == 2) {
        if (endSimulation == false) { // if we haven't already crossed an ending milestone
          if (endOnMiddleCrossing == true) { // then it's a reversal stage
            cout<<"Middle milestone crossed. Time:" << context.getTime() << " ps\n";
            endSimulation = true;
            ofstream datafile;
            datafile.open(dataFileNames[i], std::ios_base::app);
            datafile << "2 " << context.getTime() << "\n";
            datafile.close();
          } else { // Then its the forward stage, so assert that this milestone is crossed
            if (crossedStartingMilestone == false) { // as long as we haven't crossed once already
              crossedStartingMilestone = true; // we've crossed it once
              context.setTime(0.0); // reset the timer
            }
          }
        }
      } else if (h_collectionReturnCode[0] == 3) { // This should be a fairly identical procedure to inner milestone crossing above
        if (endSimulation == false) {
          cout<<"Outer milestone crossed. Time:" << context.getTime() << " ps\n";
          endSimulation = true;
          ofstream datafile;
          datafile.open(dataFileNames[i], std::ios_base::app);
          if (crossedStartingMilestone == true || endOnMiddleCrossing == true) {
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
    //setupPlanarZMilestones(force);
    //validateAndUpload();
    cout << "copyParametersToContext\n";
    
    cu.invalidateMolecules();
}























