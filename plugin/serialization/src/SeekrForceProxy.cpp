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
#include <iostream>

using namespace SeekrPlugin;
using namespace OpenMM;
using namespace std;


SeekrForceProxy::SeekrForceProxy() : SerializationProxy("SeekrForce") {
}


// This module allows everything to be converted to XML and back again

void SeekrForceProxy::serialize(const void* object, SerializationNode& node) const {
    cout << "mark10\n";
    node.setIntProperty("version", 1);
    const SeekrForce& force = *reinterpret_cast<const SeekrForce*>(object);
    SerializationNode& planarZMilestones = node.createChildNode("PlanarZMilestones");
    
	cout << "getNumPlanarZMilestones: " << force.getNumPlanarZMilestones() << endl; 

    for (int i = 0; i < force.getNumPlanarZMilestones(); i++) {
      //std::vector<int> atomIndices1;
      int atomIndex1;
      //std::vector<int> atomIndices2;
      int atomIndex2;
      SerializationNode& planarZMilestone = planarZMilestones.createChildNode("PlanarZMilestone");
      SerializationNode& atomIndicesNode = planarZMilestone.setIntProperty("numIndices1", force.getPlanarZNumIndices(i,1)).setIntProperty("numIndices2", force.getPlanarZNumIndices(i,2)).setDoubleProperty("offset1", force.getPlanarZOffset(i,1)).setDoubleProperty("offset2", force.getPlanarZOffset(i,2)).setDoubleProperty("offset3", force.getPlanarZOffset(i,3)).setDoubleProperty("endOnMiddleCrossing", force.getEndOnMiddleCrossing()); //
      planarZMilestone.setStringProperty("dataFileName", force.getPlanarZDataFileName(i));
      cout << "mark50\n";
      SerializationNode& atomIndicesNode1 = atomIndicesNode.createChildNode("atomIndex1");
      cout << "mark60\n";
      for (int j = 0; j < force.getPlanarZNumIndices(i,1); j++) {
		cout << "mark65" << endl;
        //force.getPlanarZMilestoneAtoms(i, j, atomIndices1[j], 1);
        force.getPlanarZMilestoneAtoms(i, j, atomIndex1, 1);
        cout << "mark70\n";
        //atomIndicesNode1.setIntProperty("index", atomIndices1[j]);
        atomIndicesNode1.setIntProperty("index", atomIndex1);
      }
      cout << "mark75\n";
      SerializationNode& atomIndicesNode2 = atomIndicesNode.createChildNode("atomIndex2");
      for (int j = 0; j < force.getPlanarZNumIndices(i,2); j++) {
        //force.getPlanarZMilestoneAtoms(i, j, atomIndices2[j], 2);
        force.getPlanarZMilestoneAtoms(i, j, atomIndex2, 2);
        //atomIndicesNode2.setIntProperty("index", atomIndices2[j]);
        atomIndicesNode2.setIntProperty("index", atomIndex2);
      }
    }
    
    SerializationNode& sphericalMilestones = node.createChildNode("SphericalMilestones");
    
	cout << "getNumSphericalMilestones: " << force.getNumSphericalMilestones() << endl; 

    for (int i = 0; i < force.getNumSphericalMilestones(); i++) {
      //std::vector<int> atomIndices1;
      int atomIndex1;
      //std::vector<int> atomIndices2;
      int atomIndex2;
      SerializationNode& sphericalMilestone = sphericalMilestones.createChildNode("SphericalMilestone");
      SerializationNode& atomIndicesNode = sphericalMilestone.setIntProperty("numIndices1", force.getSphericalNumIndices(i,1)).setIntProperty("numIndices2", force.getSphericalNumIndices(i,2)).setDoubleProperty("offset1", force.getSphericalRadius(i,1)).setDoubleProperty("offset2", force.getSphericalRadius(i,2)).setDoubleProperty("offset3", force.getSphericalRadius(i,3)).setDoubleProperty("endOnMiddleCrossing", force.getEndOnMiddleCrossing()); //
      sphericalMilestone.setStringProperty("dataFileName", force.getSphericalDataFileName(i));
      cout << "mark50\n";
      SerializationNode& atomIndicesNode1 = atomIndicesNode.createChildNode("atomIndex1");
      cout << "mark60\n";
      for (int j = 0; j < force.getSphericalNumIndices(i,1); j++) {
		cout << "mark65" << endl;
        //force.getSphericalMilestoneAtoms(i, j, atomIndices1[j], 1);
        force.getSphericalMilestoneAtoms(i, j, atomIndex1, 1);
        cout << "mark70\n";
        //atomIndicesNode1.setIntProperty("index", atomIndices1[j]);
        atomIndicesNode1.setIntProperty("index", atomIndex1);
      }
      cout << "mark75\n";
      SerializationNode& atomIndicesNode2 = atomIndicesNode.createChildNode("atomIndex2");
      for (int j = 0; j < force.getSphericalNumIndices(i,2); j++) {
        //force.getSphericalMilestoneAtoms(i, j, atomIndices2[j], 2);
        force.getSphericalMilestoneAtoms(i, j, atomIndex2, 2);
        //atomIndicesNode2.setIntProperty("index", atomIndices2[j]);
        atomIndicesNode2.setIntProperty("index", atomIndex2);
      }
    }
    cout << "mark99\n";
}

void* SeekrForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    cout << "mark101\n";
    SeekrForce* force = new SeekrForce();
    try {
        const SerializationNode& planarZMilestones = node.getChildNode("PlanarZMilestones");
        
        for (int i = 0; i < planarZMilestones.getChildren().size(); i++) {
          const SerializationNode& planarZMilestone = planarZMilestones.getChildNode("PlanarZMilestone");
          std::vector<int> atomIndices1;
          std::vector<int> atomIndices2;
          int numIndices1 = planarZMilestone.getIntProperty("numIndices1");
          int numIndices2 = planarZMilestone.getIntProperty("numIndices2");
          float offset1 = planarZMilestone.getDoubleProperty("offset1");
          float offset2 = planarZMilestone.getDoubleProperty("offset2");
          float offset3 = planarZMilestone.getDoubleProperty("offset3");
          const SerializationNode& atomIndexNode1 = planarZMilestone.getChildNode("atomIndex1");
          const SerializationNode& atomIndexNode2 = planarZMilestone.getChildNode("atomIndex2");
          bool endOnMiddleCrossing = planarZMilestone.getBoolProperty("endOnMiddleCrossing");
          cout << "mark150\n";
          //TODO: this is a problem
          //std::string dataFileName = "/tmp/test.txt";
          const std::string dataFileName = planarZMilestone.getStringProperty("dataFileName"); // This might not work, resort to above if necessary
          cout << "mark160\n";
          for (int j = 0; j < numIndices1; j++) {
            cout << "long command: " << atomIndexNode1.getChildren()[j].getIntProperty("index");
            atomIndices1.push_back(atomIndexNode1.getChildren()[j].getIntProperty("index"));
          }
          cout << "mark170\n";
          for (int j = 0; j < numIndices2; j++) {
            atomIndices2.push_back(atomIndexNode2.getChildren()[j].getIntProperty("index"));
          }

          force->addPlanarZMilestone(numIndices1, numIndices2, offset1, offset2, offset3, atomIndices1, atomIndices2, endOnMiddleCrossing, dataFileName);
          
        }
    }
    
     catch (...) {
        delete force;
        throw;
    }
    cout << "mark200\n";
    return force;
    
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
          cout << "mark150\n";
          //TODO: this is a problem
          //std::string dataFileName = "/tmp/test.txt";
          const std::string dataFileName = sphericalMilestone.getStringProperty("dataFileName"); // This might not work, resort to above if necessary
          cout << "mark160\n";
          for (int j = 0; j < numIndices1; j++) {
            cout << "long command: " << atomIndexNode1.getChildren()[j].getIntProperty("index");
            atomIndices1.push_back(atomIndexNode1.getChildren()[j].getIntProperty("index"));
          }
          cout << "mark170\n";
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
    cout << "mark200\n";
    return force;
}


//TODO: Get this thing to work...made it all the way to mark170




