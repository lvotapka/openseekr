#ifndef CUDA_SEEKR_KERNELS_H_
#define CUDA_SEEKR_KERNELS_H_

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

#include "SeekrKernels.h"
#include "openmm/kernels.h"
#include "openmm/System.h"
#include "openmm/cuda/CudaContext.h"
#include "openmm/cuda/CudaArray.h"


namespace SeekrPlugin {

/**
 * This kernel is invoked by SeekrForce to calculate the forces acting on the system and the energy of the system.
 */
class CudaCalcSeekrForceKernel : public CalcSeekrForceKernel {
public:
    CudaCalcSeekrForceKernel(std::string name, const OpenMM::Platform& platform, OpenMM::CudaContext& cu, const OpenMM::System& system); /*:
            CalcSeekrForceKernel(name, platform), hasInitializedKernel(false), cu(cu), system(system), params(NULL) {
    }*/
    
    ~CudaCalcSeekrForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the SeekrForce this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const SeekrForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the SeekrForce to copy the parameters from
     */
    void copyParametersToContext(OpenMM::ContextImpl& context, const SeekrForce& force);
private:
    //int numBonds;
    
    CUfunction computeSphericalMilestonesKernel;
    
    int numSphericalMilestones;
    int numSphericalAtomIndices;
    
    OpenMM::CudaArray* sphericalNumIndices1;
    std::vector<int> h_sphericalNumIndices1;
    
    OpenMM::CudaArray* sphericalNumIndices2;
    std::vector<int> h_sphericalNumIndices2;
    
    OpenMM::CudaArray* sphericalRadii1;
    std::vector<float> h_sphericalRadii1;
    
    OpenMM::CudaArray* sphericalRadii2;
    std::vector<float> h_sphericalRadii2;
    
    OpenMM::CudaArray* sphericalRadii3;
    std::vector<float> h_sphericalRadii3;
    
    OpenMM::CudaArray* sphericalAtomIndices1;
    std::vector<int> h_sphericalAtomIndices1;
    
    OpenMM::CudaArray* sphericalAtomIndices2;
    std::vector<int> h_sphericalAtomIndices2;
    
    OpenMM::CudaArray* sphericalAtomBounds1;
    std::vector<int2> h_sphericalAtomBounds1;
    
    OpenMM::CudaArray* sphericalAtomBounds2;
    std::vector<int2> h_sphericalAtomBounds2;
    
    OpenMM::CudaArray* sphericalOldCom1;
    std::vector<float4> h_sphericalOldCom1;
    
    OpenMM::CudaArray* sphericalOldCom2;
    std::vector<float4> h_sphericalOldCom2;
    
    OpenMM::CudaArray* collectionReturnCode;
    std::vector<float> h_collectionReturnCode;
    //std::vector<int> h_collectionReturnCode;
    
    bool endSimulation;
    bool endOnMiddleCrossing;
    std::vector<std::string> dataFileNames;
    
    bool hasInitializedKernel;
    OpenMM::CudaContext& cu;
    const OpenMM::System& system;
    OpenMM::CudaArray* params;
    
    void allocateMemory(const SeekrForce& force);
    void setupSphericalMilestones(const SeekrForce& force);
    void validateAndUpload();
};

} // namespace SeekrPlugin

#endif /*CUDA_SEEKR_KERNELS_H_*/
