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

#include "SeekrForceProxy.h"
#include "SeekrForce.h"
#include "openmm/serialization/SerializationNode.h"
#include <sstream>

using namespace SeekrPlugin;
using namespace OpenMM;
using namespace std;

SeekrForceProxy::SeekrForceProxy() : SerializationProxy("SeekrForce") {
}

// This module allows everything to be converted to XML and back again

void SeekrForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const SeekrForce& force = *reinterpret_cast<const SeekrForce*>(object);
    SerializationNode& sphericalMilestones = node.createChildNode("SphericalMilestones");
    
    for (int i = 0; i < force.getNumSphericalMilestones(); i++) {
      std::vector<int> atomIndices1;
      std::vector<int> atomIndices2;
      SerializationNode& sphericalMilestone = sphericalMilestones.createChildNode("SphericalMilestone");
      SerializationNode& atomIndicesNode = sphericalMilestone.setIntProperty("numIndices1", force.getSphericalNumIndices(i,1)).setIntProperty("numIndices2", force.getSphericalNumIndices(i,2)).setDoubleProperty("radius1", force.getSphericalRadius(i,1)).setDoubleProperty("radius2", force.getSphericalRadius(i,2)).setDoubleProperty("radius3", force.getSphericalRadius(i,3)).setDoubleProperty("endOnMiddleCrossing", force.getEndOnMiddleCrossing()); //.setStringProperty("dataFileName", force.getDataFileName(i));
      SerializationNode& atomIndicesNode1 = atomIndicesNode.createChildNode("atomIndex1");
      for (int j = 0; j < force.getSphericalNumIndices(i,1); j++) {
        force.getSphericalMilestoneAtoms(i, j, atomIndices1[j], 1);
        atomIndicesNode1.setIntProperty("index", atomIndices1[j]);
      }
      SerializationNode& atomIndicesNode2 = atomIndicesNode.createChildNode("atomIndex2");
      for (int j = 0; j < force.getSphericalNumIndices(i,2); j++) {
        force.getSphericalMilestoneAtoms(i, j, atomIndices2[j], 2);
        atomIndicesNode2.setIntProperty("index", atomIndices2[j]);
      }
    }
}

void* SeekrForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    SeekrForce* force = new SeekrForce();
    try {
        const SerializationNode& sphericalMilestones = node.getChildNode("SphericalMilestones");
        
        for (int i = 0; i < sphericalMilestones.getChildren().size(); i++) {
          const SerializationNode& sphericalMilestone = sphericalMilestones.getChildNode("SphericalMilestone");
          std::vector<int> atomIndices1;
          std::vector<int> atomIndices2;
          int numIndices1 = sphericalMilestone.getIntProperty("numIndices1");
          int numIndices2 = sphericalMilestone.getIntProperty("numIndices2");
          float radius1 = sphericalMilestone.getDoubleProperty("radius1");
          float radius2 = sphericalMilestone.getDoubleProperty("radius2");
          float radius3 = sphericalMilestone.getDoubleProperty("radius3");
          const SerializationNode& atomIndexNode1 = sphericalMilestone.getChildNode("atomIndex1");
          const SerializationNode& atomIndexNode2 = sphericalMilestone.getChildNode("atomIndex2");
          bool endOnMiddleCrossing = sphericalMilestone.getBoolProperty("endOnMiddleCrossing");
          std::string dataFileName = "/tmp/test.txt"; //sphericalMilestone.getStringProperty("dataFileName");
          
          for (int j = 0; j < numIndices1; j++) {
            atomIndices1.push_back(atomIndexNode1.getChildren()[j].getIntProperty("index"));
          }
          
          for (int j = 0; j < numIndices2; j++) {
            atomIndices2.push_back(atomIndexNode2.getChildren()[j].getIntProperty("index"));
          }
          
          force->addSphericalMilestone(numIndices1, numIndices2, radius1, radius2, radius3, atomIndices1, atomIndices2, endOnMiddleCrossing, dataFileName);
          
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
