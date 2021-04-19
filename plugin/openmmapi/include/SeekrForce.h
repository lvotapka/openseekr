#ifndef OPENMM_SEEKRFORCE_H_
#define OPENMM_SEEKRFORCE_H_

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

#include "openmm/Context.h"
#include "openmm/Force.h"
#include <vector>
#include "internal/windowsExportSeekr.h"

namespace SeekrPlugin {

/**
 * This class implements milestoning in OpenMM
 */

class OPENMM_EXPORT_SEEKR SeekrForce : public OpenMM::Force {
public:
    /**
     * Create a SeekrForce.
     */
    SeekrForce();
    
    int getNumSphericalMilestones() const;
    
    int getNumSphericalAtomIndices() const;
    
    int getSphericalNumIndices(int forceIndex, int molecule) const;
    
    float getSphericalRadius(int forceIndex, int milestone_id) const;
    
    void getSphericalMilestoneAtoms(int forceIndex, int atomIndex, int& atom_id, int molecule) const;
    
    bool getEndOnMiddleCrossing() const;
    
    std::string getDataFileName(int forceIndex) const;
    
    /**
     * Get the file name that the integrator writes positions/velocities to
     */
    const std::string& getSaveStateFileName() const;
    
    /**
     * Set the file name that the integrator would write positions/velocities to
     *
     * @param fileName    the string of the state file name upon crossing
     */
    void setSaveStateFileName(const std::string& fileName);
    
    void addSphericalMilestone(int numIndices1, int numIndices2, float radius1, 
                              float radius2, float radius3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2, bool endOnMiddleCrossingArg,
                              std::string dataFileName);
                              
    void modifySphericalMilestone(int forceIndex, int numIndices1, int numIndices2, float radius1, 
                              float radius2, float radius3, std::vector<int> atomIndices1,
                              std::vector<int> atomIndices2, bool endOnMiddleCrossingArg,
                              std::string dataFileName);
    
    /**
     * Update the per-bond parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setBondParameters() to modify this object's parameters, then call updateParametersInState()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-bond parameters.  The set of particles involved
     * in a bond cannot be changed, nor can new bonds be added.
     */
    void updateParametersInContext(OpenMM::Context& context);
    /**
     * Returns true if the force uses periodic boundary conditions and false otherwise. Your force should implement this
     * method appropriately to ensure that `System.usesPeriodicBoundaryConditions()` works for all systems containing
     * your force.
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
    
    
protected:
    OpenMM::ForceImpl* createImpl() const;
private:
    class SphericalMilestoneInfo;
    
    class SphericalMilestoneInfo {
      public:
          int numIndices1, numIndices2;
          float radius1, radius2, radius3;
          std::vector<int> atomIndices1;
          std::vector<int> atomIndices2;
          std::string dataFileName;
          
          SphericalMilestoneInfo() {
            numIndices1 = numIndices2 = 0;
            radius1 = radius2 = radius3 = 0.0;
            // vectors?
          }
          
          SphericalMilestoneInfo(int numIndices1, int numIndices2, float radius1, float radius2, float radius3, std::vector<int> atomIndices1, std::vector<int> atomIndices2, std::string dataFileName) : 
                                 numIndices1(numIndices1), numIndices2(numIndices2), radius1(radius1),
                                 radius2(radius2), radius3(radius3), atomIndices1(atomIndices1), 
                                 atomIndices2(atomIndices2), dataFileName(dataFileName) {
          }
    };
    
    std::vector<SphericalMilestoneInfo> sphericalMilestones;
    bool endOnMiddleCrossing;
    std::string saveStateFileName;
    
};

} // namespace SeekrPlugin

#endif /*OPENMM_SEEKRFORCE_H_*/
