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

#include "CudaSeekrKernelSources.h"

using namespace SeekrPlugin;
using namespace std;

const string CudaSeekrKernelSources::seekrForce = "/*\n"
"   Copyright 2019 by Andy Stokely and Lane Votapka\n"
"	All rights reserved\n"
"   PlanarZ Milestone CUDA Code for the OpenMM Plugin SEEKR\n"
"*/\n"
"\n"
"//path = /plugin/platforms/cuda/src/kernels\n"
"\n"
"#define ELEM_SWAP(a,b) { int t=(a);(a)=(b);(b)=t; }\n"
"\n"
"// TODO: these arguments might need to be arranged with arrays first, scalars last...\n"
"\n"
"extern \"C\" __global__ void monitorPlanarZMilestones(\n"
"                            const real4* __restrict__ posq,             // positions and charges\n"
"                            const mixed4* __restrict__ velm,             // velocities and masses\n"
"                            //const float* __restrict__ masses,           // the masses of all of the atoms\n"
"                            const int* __restrict__ numIndices1,         // number of atoms in receptor\n"
"                            const int* __restrict__ numIndices2,         // number of atoms in ligand\n"
"                            const float* __restrict__ offset1,           // offset of inner planarZ milestone\n"
"                            const float* __restrict__ offset2,           // offset of middle planarZ milestone ANDY\n"
"                            const float* __restrict__ offset3,           // offset of outer planarZ milestone\n"
"                            const int* __restrict__ atomIndices1,       // atom indices of receptor\n"
"                            const int* __restrict__ atomIndices2,       // atom indices of ligand\n"
"                            const int2* __restrict__ atomBounds1,\n"
"                            const int2* __restrict__ atomBounds2,\n"
"                            float* __restrict__ returncode,              // whether the milestone was crossed: 0 = uncrossed, 1 = crossed inner, 2 = crossed outer\n"
"                            float4* __restrict__ old_com1,              // Keeps track of the previous timestep's receptor COM to determine if the middle milestone was crossed\n"
"                            float4* __restrict__ old_com2,              // old ligand COM\n"
"                            const int numPlanarZMilestones) {           // offset of outer planarZ milestone\n"
"\n"
"    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n"
"   //													DIDN'T CHANGE ANY OF THIS CODE\n"
"  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n"
"\n"
"    for (int index=blockIdx.x*blockDim.x+threadIdx.x; index<numPlanarZMilestones; index+=blockDim.x*gridDim.x) { // Dedicate a given GPU to the milestone object\n"
"        //int atomid;\n"
"        returncode[index] = 0; // Initialize a nominal message to return to the CPU regarding milestone crossing\n"
"        \n"
"        real4 com1; // center of mass of the 'receptor'\n"
"        real4 com2; // center of mass of the 'ligand'\n"
"        float totalmass1 = 0.0; // total mass of the 'receptor' selection\n"
"        float totalmass2 = 0.0; // total mass of the 'ligand' selection\n"
"        float mass = 0.0; // temporary placeholder for masses of particles\n"
"        int i;\n"
"        int atomIndex; // temporary placeholder for atom index\n"
"        com1.x = com1.y = com1.z = 0.0;\n"
"        com2.x = com2.y = com2.z = 0.0;  // initialize the center of mass positions to 0,0,0\n"
"        \n"
"        for (i=atomBounds1[index].x; i<atomBounds1[index].y; i++) { // loop through this milestone set's atom indices\n"
"            atomIndex = atomIndices1[i]; // obtain this atom's index\n"
"            mass = 1.0 / velm[atomIndex].w; // extract the mass from the velm array\n"
"            com1 += posq[atomIndex]*mass; // add this atom's contribution to the center of mass\n"
"            totalmass1 += mass; // increment the selection's total mass by the particle's mass\n"
"        }\n"
"        com1 = com1 / totalmass1; // divide by total mass of all particles to get the normalized center of mass\n"
"        \n"
"        for (i=atomBounds2[index].x; i<atomBounds2[index].y; i++) { // see above for loop: applies to the ligand\n"
"            atomIndex = atomIndices2[i];\n"
"            mass = 1.0 / velm[atomIndex].w;\n"
"            com2 += posq[atomIndex]*mass;\n"
"            totalmass2 += mass;\n"
"        }\n"
"        com2 = com2 / totalmass2; \n"
"        //old_com2 = old_com2 / totalmass2;\n"
"        \n"
"        // TODO: find a less hacky way to deal with this\n"
"        if (old_com1[index].x == -9.0e5) { // then this is the first step, so initialize old_posq to equal posq \n"
"          returncode[index] = 4; //TODO **NO CHANGE***\n"
"          old_com1[index].x = com1.x; old_com1[index].y = com1.y; old_com1[index].z = com1.z; \n"
"          old_com2[index].x = com2.x; old_com2[index].y = com2.y; old_com2[index].z = com2.z; \n"
"        }\n"
"\n"
"	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n"
"	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\n"
"\n"
"\n"
"                \n"
"        real delta = com2.z - com1.z; // ANDY: com2.z - com1.z TODO **DONE** \n"
"        float old_delta = old_com2[index].z - old_com1[index].z; // ANDY: use .z here too TODO **DONE** \n"
"       \n"
"        old_com1[index].x = com1.x; old_com1[index].y = com1.y; old_com1[index].z = com1.z; // save current COM to be the next timestep's old_COM TODO **NO CHANGE**\n"
"        old_com2[index].x = com2.x; old_com2[index].y = com2.y; old_com2[index].z = com2.z; // TODO **NO CHANGE***\n"
"        \n"
"        if (delta < offset1) { // crossed inner milestone ANDY: change offset1 to z1 TODO **DONE**\n"
"            returncode[index] = 1;			\n"
"        } \n"
"        else if ((delta - offset2[index]*(old_delta - offset2[index]) < 0) { // This will return true if the particle has crossed the middle milestone since the last timestep ANDY TODO **DONE**\n"
"            returncode[index] = 2;\n"
"        }\n"
"        else if (DELTA > offset3[index]) { // crossed outer milestone ANDY TODO **DONE**\n"
"            returncode[index] = 3;\n"
"        \n"
"        }\n"
"        \n"
"        //returncode[index] = offset2[index]*offset2[index];\n"
"    }\n"
"}\n"
"\n"
"extern \"C\" __global__ void monitorSphericalMilestones(\n"
"                            const real4* __restrict__ posq,             // positions and charges\n"
"                            const mixed4* __restrict__ velm,             // velocities and masses\n"
"                            //const float* __restrict__ masses,           // the masses of all of the atoms\n"
"                            const int* __restrict__ numIndices1,         // number of atoms in receptor\n"
"                            const int* __restrict__ numIndices2,         // number of atoms in ligand\n"
"                            const float* __restrict__ radius1,           // radius of inner spherical milestone\n"
"                            const float* __restrict__ radius2,           // radius of middle spherical milestone\n"
"                            const float* __restrict__ radius3,           // radius of outer spherical milestone\n"
"                            const int* __restrict__ atomIndices1,       // atom indices of receptor\n"
"                            const int* __restrict__ atomIndices2,       // atom indices of ligand\n"
"                            const int2* __restrict__ atomBounds1,\n"
"                            const int2* __restrict__ atomBounds2,\n"
"                            float* __restrict__ returncode,              // whether the milestone was crossed: 0 = uncrossed, 1 = crossed inner, 2 = crossed outer\n"
"                            float4* __restrict__ old_com1,\n"
"                            float4* __restrict__ old_com2,\n"
"                            const int numSphericalMilestones) {           // radius of outer spherical milestone\n"
"    \n"
"    for (int index=blockIdx.x*blockDim.x+threadIdx.x; index<numSphericalMilestones; index+=blockDim.x*gridDim.x) {\n"
"        //int atomid;\n"
"        returncode[index] = 0;\n"
"        \n"
"        real4 com1;\n"
"        real4 com2;\n"
"        float totalmass1 = 0.0;\n"
"        float totalmass2 = 0.0;\n"
"        float mass = 0.0;\n"
"        int i;\n"
"        int atomIndex;\n"
"        com1.x = com1.y = com1.z = 0.0;\n"
"        com2.x = com2.y = com2.z = 0.0;  // initialize the center of mass positions to 0,0,0\n"
"        \n"
"        for (i=atomBounds1[index].x; i<atomBounds1[index].y; i++) {\n"
"            atomIndex = atomIndices1[i];\n"
"            mass = 1.0 / velm[atomIndex].w;\n"
"            com1 += posq[atomIndex]*mass;\n"
"            totalmass1 += mass;\n"
"        }\n"
"        com1 = com1 / totalmass1;\n"
"        \n"
"        for (i=atomBounds2[index].x; i<atomBounds2[index].y; i++) {\n"
"            atomIndex = atomIndices2[i];\n"
"            mass = 1.0 / velm[atomIndex].w;\n"
"            com2 += posq[atomIndex]*mass;\n"
"            totalmass2 += mass;\n"
"        }\n"
"        com2 = com2 / totalmass2;\n"
"        //old_com2 = old_com2 / totalmass2;\n"
"        \n"
"        if (old_com1[index].x == -9.0e5) { // then this is the first step, so initialize old_posq to equal posq\n"
"          returncode[index] = 4;\n"
"          old_com1[index].x = com1.x; old_com1[index].y = com1.y; old_com1[index].z = com1.z;\n"
"          old_com2[index].x = com2.x; old_com2[index].y = com2.y; old_com2[index].z = com2.z;\n"
"        }\n"
"                \n"
"        real4 delta = com2 - com1;\n"
"        float4 old_delta = old_com2[index] - old_com1[index];\n"
"        real distSquared = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;\n"
"        real old_distSquared = old_delta.x*old_delta.x + old_delta.y*old_delta.y + old_delta.z*old_delta.z;\n"
"        //real r = SQRT(distSquared);\n"
"        \n"
"        old_com1[index].x = com1.x; old_com1[index].y = com1.y; old_com1[index].z = com1.z;\n"
"        old_com2[index].x = com2.x; old_com2[index].y = com2.y; old_com2[index].z = com2.z;\n"
"        \n"
"        if (distSquared < radius1[index]*radius1[index]) { // crossed inner milestone\n"
"            returncode[index] = 1;\n"
"        } \n"
"        else if ((distSquared - radius2[index]*radius2[index])*(old_distSquared - radius2[index]*radius2[index]) < 0) { // This will return true if the particle has crossed the middle milestone since the last timestep\n"
"            returncode[index] = 2;\n"
"        }\n"
"        else if (distSquared > radius3[index]*radius3[index]) { // crossed outer milestone\n"
"            returncode[index] = 3;\n"
"        \n"
"        }\n"
"        \n"
"        //returncode[index] = radius2[index]*radius2[index];\n"
"    }\n"
"}\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"";
const string CudaSeekrKernelSources::vectorOps = "/**\n"
" * This file defines vector operations to simplify code elsewhere.\n"
"*/\n"
"\n"
"// Versions of make_x() that take a single value and set all components to that.\n"
"\n"
"inline __device__ int2 make_int2(int a) {\n"
"    return make_int2(a, a);\n"
"}\n"
"\n"
"inline __device__ int3 make_int3(int a) {\n"
"    return make_int3(a, a, a);\n"
"}\n"
"\n"
"inline __device__ int4 make_int4(int a) {\n"
"    return make_int4(a, a, a, a);\n"
"}\n"
"\n"
"inline __device__ float2 make_float2(float a) {\n"
"    return make_float2(a, a);\n"
"}\n"
"\n"
"inline __device__ float3 make_float3(float a) {\n"
"    return make_float3(a, a, a);\n"
"}\n"
"\n"
"inline __device__ float4 make_float4(float a) {\n"
"    return make_float4(a, a, a, a);\n"
"}\n"
"\n"
"inline __device__ double2 make_double2(double a) {\n"
"    return make_double2(a, a);\n"
"}\n"
"\n"
"inline __device__ double3 make_double3(double a) {\n"
"    return make_double3(a, a, a);\n"
"}\n"
"\n"
"inline __device__ double4 make_double4(double a) {\n"
"    return make_double4(a, a, a, a);\n"
"}\n"
"\n"
"// Negate a vector.\n"
"\n"
"inline __device__ int2 operator-(int2 a) {\n"
"    return make_int2(-a.x, -a.y);\n"
"}\n"
"\n"
"inline __device__ int3 operator-(int3 a) {\n"
"    return make_int3(-a.x, -a.y, -a.z);\n"
"}\n"
"\n"
"inline __device__ int4 operator-(int4 a) {\n"
"    return make_int4(-a.x, -a.y, -a.z, -a.w);\n"
"}\n"
"\n"
"inline __device__ float2 operator-(float2 a) {\n"
"    return make_float2(-a.x, -a.y);\n"
"}\n"
"\n"
"inline __device__ float3 operator-(float3 a) {\n"
"    return make_float3(-a.x, -a.y, -a.z);\n"
"}\n"
"\n"
"inline __device__ float4 operator-(float4 a) {\n"
"    return make_float4(-a.x, -a.y, -a.z, -a.w);\n"
"}\n"
"\n"
"inline __device__ double2 operator-(double2 a) {\n"
"    return make_double2(-a.x, -a.y);\n"
"}\n"
"\n"
"inline __device__ double3 operator-(double3 a) {\n"
"    return make_double3(-a.x, -a.y, -a.z);\n"
"}\n"
"\n"
"inline __device__ double4 operator-(double4 a) {\n"
"    return make_double4(-a.x, -a.y, -a.z, -a.w);\n"
"}\n"
"\n"
"// Add two vectors.\n"
"\n"
"inline __device__ int2 operator+(int2 a, int2 b) {\n"
"    return make_int2(a.x+b.x, a.y+b.y);\n"
"}\n"
"\n"
"inline __device__ int3 operator+(int3 a, int3 b) {\n"
"    return make_int3(a.x+b.x, a.y+b.y, a.z+b.z);\n"
"}\n"
"\n"
"inline __device__ int4 operator+(int4 a, int4 b) {\n"
"    return make_int4(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);\n"
"}\n"
"\n"
"inline __device__ float2 operator+(float2 a, float2 b) {\n"
"    return make_float2(a.x+b.x, a.y+b.y);\n"
"}\n"
"\n"
"inline __device__ float3 operator+(float3 a, float3 b) {\n"
"    return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);\n"
"}\n"
"\n"
"inline __device__ float4 operator+(float4 a, float4 b) {\n"
"    return make_float4(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);\n"
"}\n"
"\n"
"inline __device__ double2 operator+(double2 a, double2 b) {\n"
"    return make_double2(a.x+b.x, a.y+b.y);\n"
"}\n"
"\n"
"inline __device__ double3 operator+(double3 a, double3 b) {\n"
"    return make_double3(a.x+b.x, a.y+b.y, a.z+b.z);\n"
"}\n"
"\n"
"inline __device__ double4 operator+(double4 a, double4 b) {\n"
"    return make_double4(a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w);\n"
"}\n"
"\n"
"// Subtract two vectors.\n"
"\n"
"inline __device__ int2 operator-(int2 a, int2 b) {\n"
"    return make_int2(a.x-b.x, a.y-b.y);\n"
"}\n"
"\n"
"inline __device__ int3 operator-(int3 a, int3 b) {\n"
"    return make_int3(a.x-b.x, a.y-b.y, a.z-b.z);\n"
"}\n"
"\n"
"inline __device__ int4 operator-(int4 a, int4 b) {\n"
"    return make_int4(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);\n"
"}\n"
"\n"
"inline __device__ float2 operator-(float2 a, float2 b) {\n"
"    return make_float2(a.x-b.x, a.y-b.y);\n"
"}\n"
"\n"
"inline __device__ float3 operator-(float3 a, float3 b) {\n"
"    return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);\n"
"}\n"
"\n"
"inline __device__ float4 operator-(float4 a, float4 b) {\n"
"    return make_float4(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);\n"
"}\n"
"\n"
"inline __device__ double2 operator-(double2 a, double2 b) {\n"
"    return make_double2(a.x-b.x, a.y-b.y);\n"
"}\n"
"\n"
"inline __device__ double3 operator-(double3 a, double3 b) {\n"
"    return make_double3(a.x-b.x, a.y-b.y, a.z-b.z);\n"
"}\n"
"\n"
"inline __device__ double4 operator-(double4 a, double4 b) {\n"
"    return make_double4(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);\n"
"}\n"
"\n"
"// Multiply two vectors.\n"
"\n"
"inline __device__ int2 operator*(int2 a, int2 b) {\n"
"    return make_int2(a.x*b.x, a.y*b.y);\n"
"}\n"
"\n"
"inline __device__ int3 operator*(int3 a, int3 b) {\n"
"    return make_int3(a.x*b.x, a.y*b.y, a.z*b.z);\n"
"}\n"
"\n"
"inline __device__ int4 operator*(int4 a, int4 b) {\n"
"    return make_int4(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);\n"
"}\n"
"\n"
"inline __device__ float2 operator*(float2 a, float2 b) {\n"
"    return make_float2(a.x*b.x, a.y*b.y);\n"
"}\n"
"\n"
"inline __device__ float3 operator*(float3 a, float3 b) {\n"
"    return make_float3(a.x*b.x, a.y*b.y, a.z*b.z);\n"
"}\n"
"\n"
"inline __device__ float4 operator*(float4 a, float4 b) {\n"
"    return make_float4(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);\n"
"}\n"
"\n"
"inline __device__ double2 operator*(double2 a, double2 b) {\n"
"    return make_double2(a.x*b.x, a.y*b.y);\n"
"}\n"
"\n"
"inline __device__ double3 operator*(double3 a, double3 b) {\n"
"    return make_double3(a.x*b.x, a.y*b.y, a.z*b.z);\n"
"}\n"
"\n"
"inline __device__ double4 operator*(double4 a, double4 b) {\n"
"    return make_double4(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);\n"
"}\n"
"\n"
"// Divide two vectors.\n"
"\n"
"inline __device__ int2 operator/(int2 a, int2 b) {\n"
"    return make_int2(a.x/b.x, a.y/b.y);\n"
"}\n"
"\n"
"inline __device__ int3 operator/(int3 a, int3 b) {\n"
"    return make_int3(a.x/b.x, a.y/b.y, a.z/b.z);\n"
"}\n"
"\n"
"inline __device__ int4 operator/(int4 a, int4 b) {\n"
"    return make_int4(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);\n"
"}\n"
"\n"
"inline __device__ float2 operator/(float2 a, float2 b) {\n"
"    return make_float2(a.x/b.x, a.y/b.y);\n"
"}\n"
"\n"
"inline __device__ float3 operator/(float3 a, float3 b) {\n"
"    return make_float3(a.x/b.x, a.y/b.y, a.z/b.z);\n"
"}\n"
"\n"
"inline __device__ float4 operator/(float4 a, float4 b) {\n"
"    return make_float4(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);\n"
"}\n"
"\n"
"inline __device__ double2 operator/(double2 a, double2 b) {\n"
"    return make_double2(a.x/b.x, a.y/b.y);\n"
"}\n"
"\n"
"inline __device__ double3 operator/(double3 a, double3 b) {\n"
"    return make_double3(a.x/b.x, a.y/b.y, a.z/b.z);\n"
"}\n"
"\n"
"inline __device__ double4 operator/(double4 a, double4 b) {\n"
"    return make_double4(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);\n"
"}\n"
"\n"
"// += operator\n"
"\n"
"inline __device__ void operator+=(int2& a, int2 b) {\n"
"    a.x += b.x; a.y += b.y;\n"
"}\n"
"\n"
"inline __device__ void operator+=(int3& a, int3 b) {\n"
"    a.x += b.x; a.y += b.y; a.z += b.z;\n"
"}\n"
"\n"
"inline __device__ void operator+=(int4& a, int4 b) {\n"
"    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;\n"
"}\n"
"\n"
"inline __device__ void operator+=(float2& a, float2 b) {\n"
"    a.x += b.x; a.y += b.y;\n"
"}\n"
"\n"
"inline __device__ void operator+=(float3& a, float3 b) {\n"
"    a.x += b.x; a.y += b.y; a.z += b.z;\n"
"}\n"
"\n"
"inline __device__ void operator+=(float4& a, float4 b) {\n"
"    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;\n"
"}\n"
"\n"
"inline __device__ void operator+=(double2& a, double2 b) {\n"
"    a.x += b.x; a.y += b.y;\n"
"}\n"
"\n"
"inline __device__ void operator+=(double3& a, double3 b) {\n"
"    a.x += b.x; a.y += b.y; a.z += b.z;\n"
"}\n"
"\n"
"inline __device__ void operator+=(double4& a, double4 b) {\n"
"    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;\n"
"}\n"
"\n"
"// -= operator\n"
"\n"
"inline __device__ void operator-=(int2& a, int2 b) {\n"
"    a.x -= b.x; a.y -= b.y;\n"
"}\n"
"\n"
"inline __device__ void operator-=(int3& a, int3 b) {\n"
"    a.x -= b.x; a.y -= b.y; a.z -= b.z;\n"
"}\n"
"\n"
"inline __device__ void operator-=(int4& a, int4 b) {\n"
"    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;\n"
"}\n"
"\n"
"inline __device__ void operator-=(float2& a, float2 b) {\n"
"    a.x -= b.x; a.y -= b.y;\n"
"}\n"
"\n"
"inline __device__ void operator-=(float3& a, float3 b) {\n"
"    a.x -= b.x; a.y -= b.y; a.z -= b.z;\n"
"}\n"
"\n"
"inline __device__ void operator-=(float4& a, float4 b) {\n"
"    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;\n"
"}\n"
"\n"
"inline __device__ void operator-=(double2& a, double2 b) {\n"
"    a.x -= b.x; a.y -= b.y;\n"
"}\n"
"\n"
"inline __device__ void operator-=(double3& a, double3 b) {\n"
"    a.x -= b.x; a.y -= b.y; a.z -= b.z;\n"
"}\n"
"\n"
"inline __device__ void operator-=(double4& a, double4 b) {\n"
"    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;\n"
"}\n"
"\n"
"// *= operator\n"
"\n"
"inline __device__ void operator*=(int2& a, int2 b) {\n"
"    a.x *= b.x; a.y *= b.y;\n"
"}\n"
"\n"
"inline __device__ void operator*=(int3& a, int3 b) {\n"
"    a.x *= b.x; a.y *= b.y; a.z *= b.z;\n"
"}\n"
"\n"
"inline __device__ void operator*=(int4& a, int4 b) {\n"
"    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;\n"
"}\n"
"\n"
"inline __device__ void operator*=(float2& a, float2 b) {\n"
"    a.x *= b.x; a.y *= b.y;\n"
"}\n"
"\n"
"inline __device__ void operator*=(float3& a, float3 b) {\n"
"    a.x *= b.x; a.y *= b.y; a.z *= b.z;\n"
"}\n"
"\n"
"inline __device__ void operator*=(float4& a, float4 b) {\n"
"    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;\n"
"}\n"
"\n"
"inline __device__ void operator*=(double2& a, double2 b) {\n"
"    a.x *= b.x; a.y *= b.y;\n"
"}\n"
"\n"
"inline __device__ void operator*=(double3& a, double3 b) {\n"
"    a.x *= b.x; a.y *= b.y; a.z *= b.z;\n"
"}\n"
"\n"
"inline __device__ void operator*=(double4& a, double4 b) {\n"
"    a.x *= b.x; a.y *= b.y; a.z *= b.z; a.w *= b.w;\n"
"}\n"
"\n"
"// /= operator\n"
"\n"
"inline __device__ void operator/=(int2& a, int2 b) {\n"
"    a.x /= b.x; a.y /= b.y;\n"
"}\n"
"\n"
"inline __device__ void operator/=(int3& a, int3 b) {\n"
"    a.x /= b.x; a.y /= b.y; a.z /= b.z;\n"
"}\n"
"\n"
"inline __device__ void operator/=(int4& a, int4 b) {\n"
"    a.x /= b.x; a.y /= b.y; a.z /= b.z; a.w /= b.w;\n"
"}\n"
"\n"
"inline __device__ void operator/=(float2& a, float2 b) {\n"
"    a.x /= b.x; a.y /= b.y;\n"
"}\n"
"\n"
"inline __device__ void operator/=(float3& a, float3 b) {\n"
"    a.x /= b.x; a.y /= b.y; a.z /= b.z;\n"
"}\n"
"\n"
"inline __device__ void operator/=(float4& a, float4 b) {\n"
"    a.x /= b.x; a.y /= b.y; a.z /= b.z; a.w /= b.w;\n"
"}\n"
"\n"
"inline __device__ void operator/=(double2& a, double2 b) {\n"
"    a.x /= b.x; a.y /= b.y;\n"
"}\n"
"\n"
"inline __device__ void operator/=(double3& a, double3 b) {\n"
"    a.x /= b.x; a.y /= b.y; a.z /= b.z;\n"
"}\n"
"\n"
"inline __device__ void operator/=(double4& a, double4 b) {\n"
"    a.x /= b.x; a.y /= b.y; a.z /= b.z; a.w /= b.w;\n"
"}\n"
"\n"
"// Multiply a vector by a constant.\n"
"\n"
"inline __device__ int2 operator*(int2 a, int b) {\n"
"    return make_int2(a.x*b, a.y*b);\n"
"}\n"
"\n"
"inline __device__ int3 operator*(int3 a, int b) {\n"
"    return make_int3(a.x*b, a.y*b, a.z*b);\n"
"}\n"
"\n"
"inline __device__ int4 operator*(int4 a, int b) {\n"
"    return make_int4(a.x*b, a.y*b, a.z*b, a.w*b);\n"
"}\n"
"\n"
"inline __device__ int2 operator*(int a, int2 b) {\n"
"    return make_int2(a*b.x, a*b.y);\n"
"}\n"
"\n"
"inline __device__ int3 operator*(int a, int3 b) {\n"
"    return make_int3(a*b.x, a*b.y, a*b.z);\n"
"}\n"
"\n"
"inline __device__ int4 operator*(int a, int4 b) {\n"
"    return make_int4(a*b.x, a*b.y, a*b.z, a*b.w);\n"
"}\n"
"\n"
"inline __device__ float2 operator*(float2 a, float b) {\n"
"    return make_float2(a.x*b, a.y*b);\n"
"}\n"
"\n"
"inline __device__ float3 operator*(float3 a, float b) {\n"
"    return make_float3(a.x*b, a.y*b, a.z*b);\n"
"}\n"
"\n"
"inline __device__ float4 operator*(float4 a, float b) {\n"
"    return make_float4(a.x*b, a.y*b, a.z*b, a.w*b);\n"
"}\n"
"\n"
"inline __device__ float2 operator*(float a, float2 b) {\n"
"    return make_float2(a*b.x, a*b.y);\n"
"}\n"
"\n"
"inline __device__ float3 operator*(float a, float3 b) {\n"
"    return make_float3(a*b.x, a*b.y, a*b.z);\n"
"}\n"
"\n"
"inline __device__ float4 operator*(float a, float4 b) {\n"
"    return make_float4(a*b.x, a*b.y, a*b.z, a*b.w);\n"
"}\n"
"\n"
"inline __device__ double2 operator*(double2 a, double b) {\n"
"    return make_double2(a.x*b, a.y*b);\n"
"}\n"
"\n"
"inline __device__ double3 operator*(double3 a, double b) {\n"
"    return make_double3(a.x*b, a.y*b, a.z*b);\n"
"}\n"
"\n"
"inline __device__ double4 operator*(double4 a, double b) {\n"
"    return make_double4(a.x*b, a.y*b, a.z*b, a.w*b);\n"
"}\n"
"\n"
"inline __device__ double2 operator*(double a, double2 b) {\n"
"    return make_double2(a*b.x, a*b.y);\n"
"}\n"
"\n"
"inline __device__ double3 operator*(double a, double3 b) {\n"
"    return make_double3(a*b.x, a*b.y, a*b.z);\n"
"}\n"
"\n"
"inline __device__ double4 operator*(double a, double4 b) {\n"
"    return make_double4(a*b.x, a*b.y, a*b.z, a*b.w);\n"
"}\n"
"\n"
"// Divide a vector by a constant.\n"
"\n"
"inline __device__ int2 operator/(int2 a, int b) {\n"
"    return make_int2(a.x/b, a.y/b);\n"
"}\n"
"\n"
"inline __device__ int3 operator/(int3 a, int b) {\n"
"    return make_int3(a.x/b, a.y/b, a.z/b);\n"
"}\n"
"\n"
"inline __device__ int4 operator/(int4 a, int b) {\n"
"    return make_int4(a.x/b, a.y/b, a.z/b, a.w/b);\n"
"}\n"
"\n"
"inline __device__ float2 operator/(float2 a, float b) {\n"
"    float scale = 1.0f/b;\n"
"    return a*scale;\n"
"}\n"
"\n"
"inline __device__ float3 operator/(float3 a, float b) {\n"
"    float scale = 1.0f/b;\n"
"    return a*scale;\n"
"}\n"
"\n"
"inline __device__ float4 operator/(float4 a, float b) {\n"
"    float scale = 1.0f/b;\n"
"    return a*scale;\n"
"}\n"
"\n"
"inline __device__ double2 operator/(double2 a, double b) {\n"
"    double scale = 1.0/b;\n"
"    return a*scale;\n"
"}\n"
"\n"
"inline __device__ double3 operator/(double3 a, double b) {\n"
"    double scale = 1.0/b;\n"
"    return a*scale;\n"
"}\n"
"\n"
"inline __device__ double4 operator/(double4 a, double b) {\n"
"    double scale = 1.0/b;\n"
"    return a*scale;\n"
"}\n"
"\n"
"// *= operator (multiply vector by constant)\n"
"\n"
"inline __device__ void operator*=(int2& a, int b) {\n"
"    a.x *= b; a.y *= b;\n"
"}\n"
"\n"
"inline __device__ void operator*=(int3& a, int b) {\n"
"    a.x *= b; a.y *= b; a.z *= b;\n"
"}\n"
"\n"
"inline __device__ void operator*=(int4& a, int b) {\n"
"    a.x *= b; a.y *= b; a.z *= b; a.w *= b;\n"
"}\n"
"\n"
"inline __device__ void operator*=(float2& a, float b) {\n"
"    a.x *= b; a.y *= b;\n"
"}\n"
"\n"
"inline __device__ void operator*=(float3& a, float b) {\n"
"    a.x *= b; a.y *= b; a.z *= b;\n"
"}\n"
"\n"
"inline __device__ void operator*=(float4& a, float b) {\n"
"    a.x *= b; a.y *= b; a.z *= b; a.w *= b;\n"
"}\n"
"\n"
"inline __device__ void operator*=(double2& a, double b) {\n"
"    a.x *= b; a.y *= b;\n"
"}\n"
"\n"
"inline __device__ void operator*=(double3& a, double b) {\n"
"    a.x *= b; a.y *= b; a.z *= b;\n"
"}\n"
"\n"
"inline __device__ void operator*=(double4& a, double b) {\n"
"    a.x *= b; a.y *= b; a.z *= b; a.w *= b;\n"
"}\n"
"\n"
"// Dot product\n"
"\n"
"inline __device__ float dot(float3 a, float3 b) {\n"
"    return a.x*b.x+a.y*b.y+a.z*b.z;\n"
"}\n"
"\n"
"inline __device__ double dot(double3 a, double3 b) {\n"
"    return a.x*b.x+a.y*b.y+a.z*b.z;\n"
"}\n"
"\n"
"// Cross product\n"
"\n"
"inline __device__ float3 cross(float3 a, float3 b) {\n"
"    return make_float3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);\n"
"}\n"
"\n"
"inline __device__ float3 cross(float4 a, float4 b) {\n"
"    return make_float3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);\n"
"}\n"
"\n"
"inline __device__ double3 cross(double3 a, double3 b) {\n"
"    return make_double3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);\n"
"}\n"
"\n"
"inline __device__ double3 cross(double4 a, double4 b) {\n"
"    return make_double3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);\n"
"}\n"
"\n"
"// Normalize a vector\n"
"\n"
"inline __device__ float2 normalize(float2 a) {\n"
"    return a*rsqrtf(a.x*a.x+a.y*a.y);\n"
"}\n"
"\n"
"inline __device__ float3 normalize(float3 a) {\n"
"    return a*rsqrtf(a.x*a.x+a.y*a.y+a.z*a.z);\n"
"}\n"
"\n"
"inline __device__ float4 normalize(float4 a) {\n"
"    return a*rsqrtf(a.x*a.x+a.y*a.y+a.z*a.z+a.w*a.w);\n"
"}\n"
"\n"
"inline __device__ double2 normalize(double2 a) {\n"
"    return a*rsqrt(a.x*a.x+a.y*a.y);\n"
"}\n"
"\n"
"inline __device__ double3 normalize(double3 a) {\n"
"    return a*rsqrt(a.x*a.x+a.y*a.y+a.z*a.z);\n"
"}\n"
"\n"
"inline __device__ double4 normalize(double4 a) {\n"
"    return a*rsqrt(a.x*a.x+a.y*a.y+a.z*a.z+a.w*a.w);\n"
"}\n"
"\n"
"// Strip off the fourth component of a vector.\n"
"\n"
"inline __device__ float3 trimTo3(float4 v) {\n"
"    return make_float3(v.x, v.y, v.z);\n"
"}\n"
"\n"
"inline __device__ double3 trimTo3(double4 v) {\n"
"    return make_double3(v.x, v.y, v.z);\n"
"}\n"
"";