/*
   Copyright 2018 by Lane Votapka
   All rights reserved
*/

#define ELEM_SWAP(a,b) { int t=(a);(a)=(b);(b)=t; }

// TODO: these arguments might need to be arranged with arrays first, scalars last...

extern "C" __global__ void monitorSphericalMilestones(
                            const real4* __restrict__ posq,             // positions and charges
                            const mixed4* __restrict__ velm,             // velocities and masses
                            //const float* __restrict__ masses,           // the masses of all of the atoms
                            const int* __restrict__ numIndices1,         // number of atoms in receptor
                            const int* __restrict__ numIndices2,         // number of atoms in ligand
                            const float* __restrict__ radius1,           // radius of inner spherical milestone
                            const float* __restrict__ radius2,           // radius of middle spherical milestone
                            const float* __restrict__ radius3,           // radius of outer spherical milestone
                            const int* __restrict__ atomIndices1,       // atom indices of receptor
                            const int* __restrict__ atomIndices2,       // atom indices of ligand
                            const int2* __restrict__ atomBounds1,
                            const int2* __restrict__ atomBounds2,
                            float* __restrict__ returncode,              // whether the milestone was crossed: 0 = uncrossed, 1 = crossed inner, 2 = crossed outer
                            float4* __restrict__ old_com1,
                            float4* __restrict__ old_com2,
                            const int numSphericalMilestones) {           // radius of outer spherical milestone
    
    for (int index=blockIdx.x*blockDim.x+threadIdx.x; index<numSphericalMilestones; index+=blockDim.x*gridDim.x) {
        //int atomid;
        returncode[index] = 0;
        
        real4 com1;
        real4 com2;
        float totalmass1 = 0.0;
        float totalmass2 = 0.0;
        float mass = 0.0;
        int i;
        int atomIndex;
        com1.x = com1.y = com1.z = 0.0;
        com2.x = com2.y = com2.z = 0.0;  // initialize the center of mass positions to 0,0,0
        
        for (i=atomBounds1[index].x; i<atomBounds1[index].y; i++) {
            atomIndex = atomIndices1[i];
            mass = 1.0 / velm[atomIndex].w;
            com1 += posq[atomIndex]*mass;
            totalmass1 += mass;
        }
        com1 = com1 / totalmass1;
        
        for (i=atomBounds2[index].x; i<atomBounds2[index].y; i++) {
            atomIndex = atomIndices2[i];
            mass = 1.0 / velm[atomIndex].w;
            com2 += posq[atomIndex]*mass;
            totalmass2 += mass;
        }
        com2 = com2 / totalmass2;
        //old_com2 = old_com2 / totalmass2;
        
        if (old_com1[index].x == -9.0e5) { // then this is the first step, so initialize old_posq to equal posq
          returncode[index] = 4;
          old_com1[index].x = com1.x; old_com1[index].y = com1.y; old_com1[index].z = com1.z;
          old_com2[index].x = com2.x; old_com2[index].y = com2.y; old_com2[index].z = com2.z;
        }
                
        real4 delta = com2 - com1;
        float4 old_delta = old_com2[index] - old_com1[index];
        real distSquared = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
        real old_distSquared = old_delta.x*old_delta.x + old_delta.y*old_delta.y + old_delta.z*old_delta.z;
        //real r = SQRT(distSquared);
        
        old_com1[index].x = com1.x; old_com1[index].y = com1.y; old_com1[index].z = com1.z;
        old_com2[index].x = com2.x; old_com2[index].y = com2.y; old_com2[index].z = com2.z;
        
        if (distSquared < radius1[index]*radius1[index]) { // crossed inner milestone
            returncode[index] = 1;
        } 
        else if ((distSquared - radius2[index]*radius2[index])*(old_distSquared - radius2[index]*radius2[index]) < 0) { // This will return true if the particle has crossed the middle milestone since the last timestep
            returncode[index] = 2;
        }
        else if (distSquared > radius3[index]*radius3[index]) { // crossed outer milestone
            returncode[index] = 3;
        
        }
        
        if (isnan(distSquared) == true) {
            returncode[index] = 5;
        }
        
        //returncode[index] = distSquared;
    }
}

/*
extern "C" __global__ void computeDistRest(
                            const real4* __restrict__ posq,             // positions and charges
                            const int2* __restrict__ atomIndices,       // pair of atom indices
                            const float4* __restrict__ distanceBounds,  // r1, r2, r3, r4
                            const float* __restrict__ forceConstants,   // k
                            int* __restrict__ indexToGlobal,            // array of indices into global arrays
                            float* __restrict__ energies,               // global array of restraint energies
                            float3* __restrict__ forceBuffer,           // temporary buffer to hold the force
                            const int numRestraints) {
    for (int index=blockIdx.x*blockDim.x+threadIdx.x; index<numRestraints; index+=blockDim.x*gridDim.x) {
        // get my global index
        const int globalIndex = indexToGlobal[index];

        // get the distances
        const float r1 = distanceBounds[index].x;
        const float r2 = distanceBounds[index].y;
        const float r3 = distanceBounds[index].z;
        const float r4 = distanceBounds[index].w;

        // get the force constant
        const float k = forceConstants[index];

        // get atom indices and compute distance
        int atomIndexA = atomIndices[index].x;
        int atomIndexB = atomIndices[index].y;
        real4 delta = posq[atomIndexA] - posq[atomIndexB];
        real distSquared = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
        real r = SQRT(distSquared);

        // compute force and energy
        float energy = 0.0;
        float dEdR = 0.0;
        float diff = 0.0;
        float diff2 = 0.0;
        float3 f;

        if(r < r1) {
            energy = k * (r - r1) * (r1 - r2) + 0.5 * k * (r1 - r2) * (r1 - r2);
            dEdR = k * (r1 - r2);
        }
        else if(r < r2) {
            diff = r - r2;
            diff2 = diff * diff;
            energy = 0.5 * k * diff2;
            dEdR = k * diff;
        }
        else if(r < r3) {
            dEdR = 0.0;
            energy = 0.0;
        }
        else if(r < r4) {
            diff = r - r3;
            diff2 = diff * diff;
            energy = 0.5 * k * diff2;
            dEdR = k * diff;
        }
        else {
            energy = k * (r - r4) * (r4 - r3) + 0.5 * k * (r4 - r3) * (r4 - r3);
            dEdR = k * (r4 - r3);
        }

        // store force into local buffer
        if (r > 0) {
            f.x = delta.x * dEdR / r;
            f.y = delta.y * dEdR / r;
            f.z = delta.z * dEdR / r;
        } else {
            f.x = 0.0;
            f.y = 0.0;
            f.z = 0.0;
        }
        forceBuffer[index] = f;

        // store energy into global buffer
        energies[globalIndex] = energy;
    }
}
*/
