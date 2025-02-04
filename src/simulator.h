#ifndef SIMULATOR_H
#define SIMULATOR_H

// #ifndef _SIMULATOR_H_
// #define _SIMULATOR_H_

#include <list>
#include "landscape.h"
#include "individual.h"
#include <Rcpp.h>  // <-- ADD THIS LINE



using namespace std;

/* ----------------------
  Using random library from Agner Fog
  -----------------------
*/
#include "randomc/randomc.h"            // define classes for random number generators
#include "stocc/stocc.h"              // define random library classes
// define which random number generator to base random library on:

#include <time.h>
/* ----------------------
  End of include for random library from Agner Fog
  -----------------------
*/


struct TSimParam
{
 Mat_DP* land;           // matrix with landscape, each cell having a habitat quality between 0 and 1
 int initpopulation;     // initial population size
 int nsteps;             // number of steps in simulation
 int hrsize;             // home range size
 double birthrate;       // fecundity per individual (see CalculateFledglings for stochastic/determinitisc options)
 int breedingage;        // age of first breeding
 double survival;
    // In stochastic simulations survival is the annual survival and takes values in the interval ]0,1[
    // In deterministic simulations, survival is the maximum life span and takes any integer value >=1
    // The number of fledglings are also stochastic or deterministic dependent on the survival parameter
 double distanceweight;
    // Weight in the fitness (parameter c in the manuscript) of the distance of a cell to the home range center
 double dispersaldistance;
    // Depends on dispersal model:
    // Dispersal Model 0: not used
    // Dispersal Mode 1: radius of the dispersal kernel.
    // Dispersal Mode 2: number of random-walk steps
 int dispersalmode;
    // Dispersal mode can take the following values
    // 0: Global dispersal
    // 1: Local dispersal, habitat search in a local kernel
    // 2: Random walk
 double sinkavoidance;
    // Probability of avoiding sink habitats (habitat quality 0) during dispersal. Takes values between 0 and 1.
 double neighavoidance;
    // Probability of avoiding neighbors during dispersal. Takes values between 0 and 1.
 double sinkmortality;
    // Probability (per dispersal step) of dying in a sink habitat (habitat quality 0, eg. roads)
 Rcpp::List disturbance_matrices; 
    
 string filename;
};

class TSimulator
{
 private:
        long double Growth(int x, long double nx);
        void OutputGeneration();
        void OutputParameters();
        // data members
        TPopulation population;  //population of settlers
        TLandscape* landscape;
        int nsteps;
        unsigned int hrsize;
        double birthrate;
        unsigned int breedingage;
        double survival;
        long int initpopulation;
        double distanceweight;
        double dispersaldistance;
        int dispersalmode;
        int step;
        double sinkavoidance;
        double neighavoidance;
        double sinkmortality;
        string filename;
        Rcpp::List disturbance_matrices; 
        double optimalfitness;
 public:
        TSimulator(const TSimParam&);
        ~TSimulator();
        void Step();
        unsigned int GetHomeRangeSize() {return hrsize;}
        double GetDistanceWeight() {return distanceweight;}
        double GetBirthRate() {return birthrate;}
        double GetSurvival() {return survival;}
        unsigned int GetBreedingAge() {return breedingage;}
        TLandscape* GetLandscape() {return landscape;}
        double GetOptimalFitness() {return optimalfitness;}
        const char* GetFileName() {return filename.c_str();}
        long GetPopulationSize() {return population.size();}
        int GetDispersalMode() {return dispersalmode;}
        double GetDispersalDistance() {return dispersaldistance;}
        double GetSinkAvoidance() {return sinkavoidance;}
        double GetNeighAvoidance() {return neighavoidance;}
        double GetSinkMortality() {return sinkmortality;}
        Rcpp::List GetDisturbanceMatrices() { return disturbance_matrices; }  // <-- Getter function
        int GetStep() {return step;}
		int GetNSteps() {return nsteps;}
        StochasticLib1* sto;
};
#endif
