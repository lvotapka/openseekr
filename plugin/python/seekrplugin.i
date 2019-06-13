/*
   Copyright 2018 by Lane Votapka
   All rights reserved
*/

%module seekrplugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include "SeekrForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import simtk.openmm as mm
import simtk.unit as unit
%}

/*
 * Add units to function outputs.
*/

%pythonappend SeekrPlugin::SeekrForce::getSphericalNumIndices(int forceIndex, int molecule) const %{
    
%}

%pythonappend SeekrPlugin::SeekrForce::getSphericalRadius(int forceIndex, int milestone_id) const %{
    
%}

%pythonappend SeekrPlugin::SeekrForce::getSphericalMilestoneAtoms(int forceIndex, int atomIndex, int atom_id, int molecule) const %{
    
%}

%pythonappend SeekrPlugin::SeekrForce::getEndOnMiddleCrossing() const %{

%}

%pythonappend SeekrPlugin::SeekrForce::getSphericalDataFileName(int forceIndex) const %{

%}

%pythonappend SeekrPlugin::SeekrForce::getPlanarZNumIndices(int forceIndex, int molecule) const %{
    
%}

%pythonappend SeekrPlugin::SeekrForce::getPlanarZOffset(int forceIndex, int milestone_id) const %{
    
%}

%pythonappend SeekrPlugin::SeekrForce::getPlanarZMilestoneAtoms(int forceIndex, int atomIndex, int atom_id, int molecule) const %{
    
%}

%pythonappend SeekrPlugin::SeekrForce::getEndOnMiddleCrossing() const %{

%}

%pythonappend SeekrPlugin::SeekrForce::getPlanarZDataFileName(int forceIndex) const %{

%}


namespace SeekrPlugin {

class SeekrForce : public OpenMM::Force {
public:
    SeekrForce();
    
    int getSphericalNumIndices(int forceIndex, int molecule) const;
    
    float getSphericalRadius(int forceIndex, int milestone_id) const;
    
    void getSphericalMilestoneAtoms(int forceIndex, int atomIndex, int& index, int molecule) const;
    
    bool getEndOnMiddleCrossing() const;
    
    std::string getSphericalDataFileName(int foceIndex) const;
    //void getDataFileName(int forceIndex) const;
    
    void updateParametersInContext(OpenMM::Context& context);
    
    
    /*
     * The reference parameters to this function are output values.
     * Marking them as such will cause swig to return a tuple.
    */
    
    void addSphericalMilestone(int numIndices1, int numIndices2, float radius1,
          float radius2, float radius3, std::vector<int>atomIndices1, std::vector<int>atomIndices2,
          bool argEndOnMiddleCrossing, std::string dataFileName);
    
    void modifySphericalMilestone(int forceIndex, int numIndices1, int numIndices2, float radius1,
          float radius2, float radius3, std::vector<int>atomIndices1, std::vector<int>atomIndices2,
          bool argEndOnMiddleCrossing, std::string dataFileName);
    int getPlanarZNumIndices(int forceIndex, int molecule) const;
    
    float getPlanarZOffset(int forceIndex, int milestone_id) const;
    
    void getPlanarZMilestoneAtoms(int forceIndex, int atomIndex, int& index, int molecule) const;
    
    bool getEndOnMiddleCrossing() const;
    
    std::string getPlanarZDataFileName(int foceIndex) const;
    //void getDataFileName(int forceIndex) const;
    
    void updateParametersInContext(OpenMM::Context& context);
    
    
    /*
     * The reference parameters to this function are output values.
     * Marking them as such will cause swig to return a tuple.
    */
    
    void addPlanarZMilestone(int numIndices1, int numIndices2, float offset1,
          float offset2, float offset3, std::vector<int>atomIndices1, std::vector<int>atomIndices2,
          bool argEndOnMiddleCrossing, std::string dataFileName);
    
    void modifyPlanarZMilestone(int forceIndex, int numIndices1, int numIndices2, float offset1,
          float offset2, float offset3, std::vector<int>atomIndices1, std::vector<int>atomIndices2,
          bool argEndOnMiddleCrossing, std::string dataFileName);
};

}
