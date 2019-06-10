/*
   Copyright 2019 by Andy Stokely
	All rights reserved
   Planar Milestone CUDA Code for the OpenMM Plugin SEEKR
*/

//path = /plugin/platforms/cuda/src/kernels

#define ELEM_SWAP(a,b) { int t=(a);(a)=(b);(b)=t; }

// TODO: these arguments might need to be arranged with arrays first, scalars last...

extern "C" __global__ void monitorPlanarMilestones(
                            const real4* __restrict__ posq,             // positions and charges
                            const mixed4* __restrict__ velm,             // velocities and masses
                            //const float* __restrict__ masses,           // the masses of all of the atoms
                            const int* __restrict__ numIndices1,         // number of atoms in receptor
                            const int* __restrict__ numIndices2,         // number of atoms in ligand
                            const float* __restrict__ length1,           // length of inner planar milestone
                            const float* __restrict__ length2,           // length of middle planar milestone ANDY
                            const float* __restrict__ length3,           // length of outer planar milestone
                            const int* __restrict__ atomIndices1,       // atom indices of receptor
                            const int* __restrict__ atomIndices2,       // atom indices of ligand
                            const int2* __restrict__ atomBounds1,
                            const int2* __restrict__ atomBounds2,
                            float* __restrict__ returncode,              // whether the milestone was crossed: 0 = uncrossed, 1 = crossed inner, 2 = crossed outer
                            float4* __restrict__ old_com1,              // Keeps track of the previous timestep's receptor COM to determine if the middle milestone was crossed
                            float4* __restrict__ old_com2,              // old ligand COM
                            const int numPlanarMilestones) {           // length of outer planar milestone

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //													DIDN'T CHANGE ANY OF THIS CODE
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int index=blockIdx.x*blockDim.x+threadIdx.x; index<numPlanarMilestones; index+=blockDim.x*gridDim.x) { // Dedicate a given GPU to the milestone object
        //int atomid;
        returncode[index] = 0; // Initialize a nominal message to return to the CPU regarding milestone crossing
        
        real4 com1; // center of mass of the 'receptor'
        real4 com2; // center of mass of the 'ligand'
        float totalmass1 = 0.0; // total mass of the 'receptor' selection
        float totalmass2 = 0.0; // total mass of the 'ligand' selection
        float mass = 0.0; // temporary placeholder for masses of particles
        int i;
        int atomIndex; // temporary placeholder for atom index
        com1.x = com1.y = com1.z = 0.0;
        com2.x = com2.y = com2.z = 0.0;  // initialize the center of mass positions to 0,0,0
        
        for (i=atomBounds1[index].x; i<atomBounds1[index].y; i++) { // loop through this milestone set's atom indices
            atomIndex = atomIndices1[i]; // obtain this atom's index
            mass = 1.0 / velm[atomIndex].w; // extract the mass from the velm array
            com1 += posq[atomIndex]*mass; // add this atom's contribution to the center of mass
            totalmass1 += mass; // increment the selection's total mass by the particle's mass
        }
        com1 = com1 / totalmass1; // divide by total mass of all particles to get the normalized center of mass
        
        for (i=atomBounds2[index].x; i<atomBounds2[index].y; i++) { // see above for loop: applies to the ligand
            atomIndex = atomIndices2[i];
            mass = 1.0 / velm[atomIndex].w;
            com2 += posq[atomIndex]*mass;
            totalmass2 += mass;
        }
        com2 = com2 / totalmass2; 
        //old_com2 = old_com2 / totalmass2;
        
        // TODO: find a less hacky way to deal with this
        if (old_com1[index].x == -9.0e5) { // then this is the first step, so initialize old_posq to equal posq 
          returncode[index] = 4; //TODO **NO CHANGE***
          old_com1[index].x = com1.x; old_com1[index].y = com1.y; old_com1[index].z = com1.z; 
          old_com2[index].x = com2.x; old_com2[index].y = com2.y; old_com2[index].z = com2.z; 
        }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                
        real delta = com2.z - com1.z; // ANDY: com2.z - com1.z TODO **DONE** 
        float old_delta = old_com2[index].z - old_com1[index].z; // ANDY: use .z here too TODO **DONE** 
       
        old_com1[index].x = com1.x; old_com1[index].y = com1.y; old_com1[index].z = com1.z; // save current COM to be the next timestep's old_COM TODO **NO CHANGE**
        old_com2[index].x = com2.x; old_com2[index].y = com2.y; old_com2[index].z = com2.z; // TODO **NO CHANGE***
        
        if (delta < length1) { // crossed inner milestone ANDY: change length1 to z1 TODO **DONE**
            returncode[index] = 1;			
        } 
        else if ((delta - length2[index]*(old_delta - length2[index]) < 0) { // This will return true if the particle has crossed the middle milestone since the last timestep ANDY TODO **DONE**
            returncode[index] = 2;
        }
        else if (DELTA > length3[index]) { // crossed outer milestone ANDY TODO **DONE**
            returncode[index] = 3;
        
        }
        
        //returncode[index] = length2[index]*length2[index];
    }
}
0


























