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
    SerializationNode& planarMilestones = node.createChildNode("PlanarMilestones");
    
	cout << "getNumPlanarMilestones: " << force.getNumPlanarMilestones() << endl; 

    for (int i = 0; i < force.getNumPlanarMilestones(); i++) {
      //std::vector<int> atomIndices1;
      int atomIndex1;
      //std::vector<int> atomIndices2;
      int atomIndex2;
      SerializationNode& planarMilestone = planarMilestones.createChildNode("PlanarMilestone");
      SerializationNode& atomIndicesNode = planarMilestone.setIntProperty("numIndices1", force.getPlanarNumIndices(i,1)).setIntProperty("numIndices2", force.getPlanarNumIndices(i,2)).setDoubleProperty("length1", force.getPlanarLength(i,1)).setDoubleProperty("length2", force.getPlanarLength(i,2)).setDoubleProperty("length3", force.getPlanarLength(i,3)).setDoubleProperty("endOnMiddleCrossing", force.getEndOnMiddleCrossing()); //
      planarMilestone.setStringProperty("dataFileName", force.getDataFileName(i));
      cout << "mark50\n";
      SerializationNode& atomIndicesNode1 = atomIndicesNode.createChildNode("atomIndex1");
      cout << "mark60\n";
      for (int j = 0; j < force.getPlanarNumIndices(i,1); j++) {
		cout << "mark65" << endl;
        //force.getPlanarMilestoneAtoms(i, j, atomIndices1[j], 1);
        force.getPlanarMilestoneAtoms(i, j, atomIndex1, 1);
        cout << "mark70\n";
        //atomIndicesNode1.setIntProperty("index", atomIndices1[j]);
        atomIndicesNode1.setIntProperty("index", atomIndex1);
      }
      cout << "mark75\n";
      SerializationNode& atomIndicesNode2 = atomIndicesNode.createChildNode("atomIndex2");
      for (int j = 0; j < force.getPlanarNumIndices(i,2); j++) {
        //force.getPlanarMilestoneAtoms(i, j, atomIndices2[j], 2);
        force.getPlanarMilestoneAtoms(i, j, atomIndex2, 2);
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
        const SerializationNode& planarMilestones = node.getChildNode("PlanarMilestones");
        
        for (int i = 0; i < planarMilestones.getChildren().size(); i++) {
          const SerializationNode& planarMilestone = planarMilestones.getChildNode("PlanarMilestone");
          std::vector<int> atomIndices1;
          std::vector<int> atomIndices2;
          int numIndices1 = planarMilestone.getIntProperty("numIndices1");
          int numIndices2 = planarMilestone.getIntProperty("numIndices2");
          float length1 = planarMilestone.getDoubleProperty("length1");
          float length2 = planarMilestone.getDoubleProperty("length2");
          float length3 = planarMilestone.getDoubleProperty("length3");
          const SerializationNode& atomIndexNode1 = planarMilestone.getChildNode("atomIndex1");
          const SerializationNode& atomIndexNode2 = planarMilestone.getChildNode("atomIndex2");
          bool endOnMiddleCrossing = planarMilestone.getBoolProperty("endOnMiddleCrossing");
          cout << "mark150\n";
          //TODO: this is a problem
          //std::string dataFileName = "/tmp/test.txt";
          const std::string dataFileName = planarMilestone.getStringProperty("dataFileName"); // This might not work, resort to above if necessary
          cout << "mark160\n";
          for (int j = 0; j < numIndices1; j++) {
            cout << "long command: " << atomIndexNode1.getChildren()[j].getIntProperty("index");
            atomIndices1.push_back(atomIndexNode1.getChildren()[j].getIntProperty("index"));
          }
          cout << "mark170\n";
          for (int j = 0; j < numIndices2; j++) {
            atomIndices2.push_back(atomIndexNode2.getChildren()[j].getIntProperty("index"));
          }

          force->addPlanarMilestone(numIndices1, numIndices2, length1, length2, length3, atomIndices1, atomIndices2, endOnMiddleCrossing, dataFileName);
          
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




