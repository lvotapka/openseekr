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
#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace SeekrPlugin;
using namespace OpenMM;
using namespace std;


extern "C" void registerSeekrSerializationProxies();

void testPlanarZSerialization() {
    // Create a Force.
	
    std::vector<int> atomIndices1 {23, 34};
    std::vector<int> atomIndices2 {45, 56, 67};
    SeekrForce force;
    std::string testDataFileName = "/tmp/test.dat";
    force.addPlanarZMilestone(2, 3, 2.71, 3.01, 3.14, atomIndices1, atomIndices2, false, testDataFileName);

    /*
    force.addBond(0, 1, 1.0, 2.0);
    force.addBond(0, 2, 2.0, 2.1);
    force.addBond(2, 3, 3.0, 2.2);
    force.addBond(5, 1, 4.0, 2.3);
    */

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<SeekrForce>(&force, "Force", buffer);
    SeekrForce* copy = XmlSerializer::deserialize<SeekrForce>(buffer);

    // Compare the two forces to see if they are identical.
	
    SeekrForce& force2 = *copy;

    ASSERT_EQUAL(force.getPlanarZNumIndices(0,1), force2.getPlanarZNumIndices(0,1));
    ASSERT_EQUAL(force.getPlanarZNumIndices(0,2), force2.getPlanarZNumIndices(0,2));
    ASSERT_EQUAL(force.getPlanarZOffset(0,1), force2.getPlanarZOffset(0,1));
    ASSERT_EQUAL(force.getPlanarZOffset(0,2), force2.getPlanarZOffset(0,2));
    ASSERT_EQUAL(force.getPlanarZOffset(0,3), force2.getPlanarZOffset(0,3))
    for (int i = 0; i < force.getPlanarZNumIndices(0,1); i++) {
        int atomIndex, atomIndexCopy;
        force.getPlanarZMilestoneAtoms(0, i, atomIndex, 1);
        force2.getPlanarZMilestoneAtoms(0, i, atomIndexCopy, 1);
        ASSERT_EQUAL(atomIndex, atomIndexCopy);
    }
    for (int i = 0; i < force.getPlanarZNumIndices(0,2); i++) {
        int atomIndex, atomIndexCopy;
        force.getPlanarZMilestoneAtoms(0, i, atomIndex, 2);
        force2.getPlanarZMilestoneAtoms(0, i, atomIndexCopy, 2);
        ASSERT_EQUAL(atomIndex, atomIndexCopy);
    }
}

void testSphericalSerialization() {
    // Create a Force.
    std::vector<int> atomIndices1 {23, 34};
    std::vector<int> atomIndices2 {45, 56, 67};
    SeekrForce force;
    std::string testDataFileName = "/tmp/test.dat";
    force.addSphericalMilestone(2, 3, 2.71, 3.01, 3.14, atomIndices1, atomIndices2, false, testDataFileName);
    /*
    force.addBond(0, 1, 1.0, 2.0);
    force.addBond(0, 2, 2.0, 2.1);
    force.addBond(2, 3, 3.0, 2.2);
    force.addBond(5, 1, 4.0, 2.3);
    */

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<SeekrForce>(&force, "Force", buffer);
    SeekrForce* copy = XmlSerializer::deserialize<SeekrForce>(buffer);

    // Compare the two forces to see if they are identical.

    SeekrForce& force2 = *copy;
    ASSERT_EQUAL(force.getSphericalNumIndices(0,1), force2.getSphericalNumIndices(0,1));
    ASSERT_EQUAL(force.getSphericalNumIndices(0,2), force2.getSphericalNumIndices(0,2));
    ASSERT_EQUAL(force.getSphericalRadius(0,1), force2.getSphericalRadius(0,1));
    ASSERT_EQUAL(force.getSphericalRadius(0,2), force2.getSphericalRadius(0,2));
    ASSERT_EQUAL(force.getSphericalRadius(0,3), force2.getSphericalRadius(0,3))
    for (int i = 0; i < force.getSphericalNumIndices(0,1); i++) {
        int atomIndex, atomIndexCopy;
        force.getSphericalMilestoneAtoms(0, i, atomIndex, 1);
        force2.getSphericalMilestoneAtoms(0, i, atomIndexCopy, 1);
        ASSERT_EQUAL(atomIndex, atomIndexCopy);
    }
    for (int i = 0; i < force.getSphericalNumIndices(0,2); i++) {
        int atomIndex, atomIndexCopy;
        force.getSphericalMilestoneAtoms(0, i, atomIndex, 2);
        force2.getSphericalMilestoneAtoms(0, i, atomIndexCopy, 2);
        ASSERT_EQUAL(atomIndex, atomIndexCopy);
    }
}



int main() {
    try {
        registerSeekrSerializationProxies(); 
        testPlanarZSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    
    try {
        registerSeekrSerializationProxies();
        testSphericalSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
