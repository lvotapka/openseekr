/*
   Copyright 2019 by Andy Stokely
	All rights reserved
   Planar Milestone C++ Code for the OpenMM Plugin SEEKR that links the seekrForce.cu CUDA code to 
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

//TODO I just updated the "spherical", "Length", etc names to their respective planar versions. I may have to go back and update the data structure types. (6-2-2019)

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
    numPlanarMilestones = 0;
    numPlanarAtomIndices = 0;
    
    planarNumIndices1 = nullptr;
    planarNumIndices2 = nullptr;
    planarLengths1 = nullptr;
    planarLengths2 = nullptr;
    planarLengths3 = nullptr;
    planarAtomIndices1 = nullptr;
    planarAtomIndices2 = nullptr; //TODO: ANDY-->RESOLVED
    planarAtomBounds1 = nullptr;
    planarAtomBounds2 = nullptr;
    planarOldCom1 = nullptr;
    planarOldCom2 = nullptr;
    
    planarCollectionReturnCode = nullptr;
    
    endSimulation = false;
    endOnMiddleCrossing = false;
    crossedStartingMilestone = false;
}

CudaCalcSeekrForceKernel::~CudaCalcSeekrForceKernel() {
    cu.setAsCurrent();
    
    delete planarNumIndices1;
    delete planarNumIndices2;
    delete planarLengths1;
    delete planarLengths2;
    delete planarLengths3;
    delete planarAtomIndices1;
    delete planarAtomIndices2; //TODO: ANDY-->RESOLVED
    delete planarAtomBounds1;
    delete planarAtomBounds2;
    delete planarOldCom1;
    delete planarOldCom2;
    delete planarCollectionReturnCode;
    //if (params != NULL)
    //    delete params;
}

void CudaCalcSeekrForceKernel::allocateMemory(const SeekrForce& force) {
  
  numPlanarMilestones = force.getNumPlanarMilestones();
  numPlanarAtomIndices = force.getNumPlanarAtomIndices();
  
  planarNumIndices1 =  CudaArray::create<int> (cu, numPlanarMilestones, "planarNumIndices1");
  planarNumIndices2 =  CudaArray::create<int> (cu, numPlanarMilestones, "planarNumIndices2");
  planarLengths1 =       CudaArray::create<float> (cu, numPlanarMilestones, "planarLengths1");
  planarLengths2 =       CudaArray::create<float> (cu, numPlanarMilestones, "planarLengths2");
  planarLengths3 =       CudaArray::create<float> (cu, numPlanarMilestones, "planarLengths3");
  planarAtomIndices1 = CudaArray::create<int> (cu, numPlanarAtomIndices, "planarAtomIndices1");
  planarAtomIndices2 = CudaArray::create<int> (cu, numPlanarAtomIndices, "planarAtomIndices2"); // TODO: ANDY
  planarAtomBounds1  = CudaArray::create<int2> (cu, numPlanarMilestones, "planarAtomBounds1");
  planarAtomBounds2  = CudaArray::create<int2> (cu, numPlanarMilestones, "planarAtomBounds2");
  planarOldCom1      = CudaArray::create<float4> (cu, numPlanarMilestones, "planarOldCom1");
  planarOldCom2      = CudaArray::create<float4> (cu, numPlanarMilestones, "planarOldCom2");
  planarCollectionReturnCode  =  CudaArray::create<float> (cu, numPlanarMilestones, "planarCollectionReturnCode"); // TODO: ANDY watch out
  //planarCollectionReturnCode =  CudaArray::create<int> (cu, numPlanarMilestones, "planarCollectionReturnCode");
  
  // host memory (CPU)
  
  h_planarNumIndices1  = std::vector<int> (numPlanarMilestones, 0);
  h_planarNumIndices2  = std::vector<int> (numPlanarMilestones, 0);
  h_planarLengths1       = std::vector<float> (numPlanarMilestones, 0.0);
  h_planarLengths2       = std::vector<float> (numPlanarMilestones, 0.0);
  h_planarLengths3       = std::vector<float> (numPlanarMilestones, 0.0);
  h_planarAtomIndices1 = std::vector<int> (numPlanarAtomIndices, 0);
  h_planarAtomIndices2 = std::vector<int> (numPlanarAtomIndices, 0);
  h_planarAtomBounds1  = std::vector<int2> (numPlanarMilestones, make_int2(-1, -1));
  h_planarAtomBounds2  = std::vector<int2> (numPlanarMilestones, make_int2(-1, -1));
  h_planarOldCom1      = std::vector<float4> (numPlanarMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
  h_planarOldCom2      = std::vector<float4> (numPlanarMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
  
  h_planarCollectionReturnCode  = std::vector<float> (numPlanarMilestones, 0);
  //h_planarCollectionReturnCode  = std::vector<int> (1, 0);
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

void checkPlanarMilestone(float Lengths1, float Lengths2, float Lengths3) { // TODO: ANDY make one for planar 
  bool error = false;
  if (Lengths1 >= Lengths2) {
    error = true;
  }
  if (Lengths2 >= Lengths3) {
    error = true;
  }
  
  if (error){
    std::stringstream m;
    m<<"Error: bad Lengthss for milestone of type: planar milestone";
    m<<". Lengths1 (inner): "<<Lengths1<<" Lengths2 (middle): "<<Lengths2<<" Lengths3 (outer): "<<Lengths3;
    throw OpenMMException(m.str());
  }
}

void CudaCalcSeekrForceKernel::setupPlanarMilestones(const SeekrForce& force){ //TODO: ANDY create your own DONE
  int numAtoms = system.getNumParticles();
  std::string milestoneType = "planar milestone";
  int currentAtomIndex1 = 0;
  int currentAtomIndex2 = 0;
  endOnMiddleCrossing = force.getEndOnMiddleCrossing();
  for (int i=0; i < numPlanarMilestones; i++) {
    
    int thisStart1 = currentAtomIndex1;
    int thisStart2 = currentAtomIndex2;
    // retrieve the planar milestone parameters from the C++ api
    h_planarNumIndices1[i] = force.getPlanarNumIndices(i, 1); 
    h_planarNumIndices2[i] = force.getPlanarNumIndices(i, 2);
    h_planarLengths1[i] = force.getPlanarLength(i, 1);
    h_planarLengths2[i] = force.getPlanarLength(i, 2);
    h_planarLengths3[i] = force.getPlanarLength(i, 3);
    
    // check the Lengthss to make sure that they make sense
    checkPlanarMilestone(h_planarLengths1[i], h_planarLengths2[i], h_planarLengths3[i]); 
    
    // Retrieve the individual atoms from the API
    for(int j=0; j<h_planarNumIndices1[i]; j++) {
      int atom_id;
      force.getPlanarMilestoneAtoms(i, j, atom_id, 1);
      checkAtomIndex(numAtoms, milestoneType, atom_id);
      h_planarAtomIndices1[currentAtomIndex1] = atom_id;
      currentAtomIndex1++;
    }
    int thisEnd1 = currentAtomIndex1;
    h_planarAtomBounds1[i] = make_int2(thisStart1, thisEnd1);
    
    for(int j=0; j<h_planarNumIndices2[i]; j++) {
      int atom_id;
      force.getPlanarMilestoneAtoms(i, j, atom_id, 2);
      checkAtomIndex(numAtoms, milestoneType, atom_id);
      h_planarAtomIndices2[currentAtomIndex2] = atom_id;
      currentAtomIndex2++;
    }
    int thisEnd2 = currentAtomIndex2;
    h_planarAtomBounds2[i] = make_int2(thisStart2, thisEnd2);
    
    /*
    cout << "CUDA kernel: Setting up planar milestone number: " << i << ". planarNumIndices1: ";
    cout << h_planarNumIndices1[i] << " planarNumIndices2: " << h_planarNumIndices2[i];
    cout << " planarLengths: (" << h_planarLengths1[i] << ", " << h_planarLengths2[i];
    cout << ", " << h_planarLengths3[i] << "). planarAtomBounds1: (";
    cout << h_planarAtomBounds1[i].x << ", " << h_planarAtomBounds1[i].y;
    cout << ") planarAtomIndices1: [";
    for (int j=h_planarAtomBounds1[i].x; j<h_planarAtomBounds1[i].y; j++) {
      cout << h_planarAtomIndices1[j] << " ";
    }
    cout << "] planarAtomBounds2: (";
    cout << h_planarAtomBounds2[i].x << ", " << h_planarAtomBounds2[i].y;
    cout << ") planarAtomIndices2: [";
    for (int j=h_planarAtomBounds2[i].x; j<h_planarAtomBounds2[i].y; j++) {
      cout << h_planarAtomIndices2[j] << " ";
    }
    cout << "]\n";*/
    dataFileNames.push_back(force.getDataFileName(i));
  }
}

void CudaCalcSeekrForceKernel::validateAndUpload() {
  cout << "Attempting to upload host arrays to device.\n";
  planarNumIndices1->upload(h_planarNumIndices1);
  planarNumIndices2->upload(h_planarNumIndices2);
  planarLengths1->upload(h_planarLengths1);
  planarLengths2->upload(h_planarLengths2);
  planarLengths3->upload(h_planarLengths3);
  planarAtomIndices1->upload(h_planarAtomIndices1);
  planarAtomIndices2->upload(h_planarAtomIndices2); // TODO: ANDY-->RESOLVED
  planarAtomBounds1->upload(h_planarAtomBounds1);
  planarAtomBounds2->upload(h_planarAtomBounds2);
  planarOldCom1->upload(h_planarOldCom1);
  planarOldCom2->upload(h_planarOldCom2);
  cout << "Uploaded all host arrays to device.\n";
}

void CudaCalcSeekrForceKernel::initialize(const System& system, const SeekrForce& force) {
    cu.setAsCurrent();
    
    allocateMemory(force);
    setupPlanarMilestones(force); // TODO:ANDY add yours here--RESOLVED
    validateAndUpload();
    
    CUmodule module = cu.createModule(CudaSeekrKernelSources::vectorOps + CudaSeekrKernelSources::seekrForce);
    computePlanarMilestonesKernel = cu.getKernel(module, "monitorPlanarMilestones");
    cu.addForce(new CudaSeekrForceInfo(force));
}

double CudaCalcSeekrForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (context.getTime() == 0.0) {
      endSimulation = false;
      h_planarOldCom1      = std::vector<float4> (numPlanarMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
      h_planarOldCom2      = std::vector<float4> (numPlanarMilestones, make_float4(-9.0e5, -9.0e5, -9.0e5, -9.0e5));
      planarOldCom1->upload(h_planarOldCom1);
      planarOldCom2->upload(h_planarOldCom2);
    }
    
    void* planarMilestoneArgs[] = {
      &cu.getPosq().getDevicePointer(),
      &cu.getVelm().getDevicePointer(),
      &planarNumIndices1->getDevicePointer(),
      &planarNumIndices2->getDevicePointer(),
      &planarLengths1->getDevicePointer(),
      &planarLengths2->getDevicePointer(),
      &planarLengths3->getDevicePointer(),
      &planarAtomIndices1->getDevicePointer(),
      &planarAtomIndices2->getDevicePointer(),
      &planarAtomBounds1->getDevicePointer(),
      &planarAtomBounds2->getDevicePointer(),
      &planarCollectionReturnCode->getDevicePointer(),
      &planarOldCom1->getDevicePointer(),
      &planarOldCom2->getDevicePointer(),
      &numPlanarMilestones
    };
    cu.executeKernel(computePlanarMilestonesKernel, planarMilestoneArgs, numPlanarMilestones);
    planarCollectionReturnCode->download(h_planarCollectionReturnCode); 
    //cout << "planarCollectionReturnCode: " << h_planarCollectionReturnCode[0] << endl;
    
    //cout << "step:" << context.getTime() << "\n";
    
    
    for (int i=0; i<numPlanarMilestones; i++) {
      if (h_planarCollectionReturnCode[i] == 1) {
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
      } else if (h_planarCollectionReturnCode[i] == 2) {
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
      } else if (h_planarCollectionReturnCode[i] == 3) { // This should be a fairly identical procedure to inner milestone crossing above
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
    //setupPlanarMilestones(force);
    //validateAndUpload();
    cout << "copyParametersToContext\n";
    
    cu.invalidateMolecules();
}























