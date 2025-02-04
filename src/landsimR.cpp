#include <Rcpp.h>
using namespace Rcpp;
#include "simulator.h"
#include "nrtypes.h"

// indmodelseC is the function called from R with the respective parameters

//' @useDynLib PopIn
//' @importFrom Rcpp sourceCpp

// [[Rcpp::export]]
List indmodelseC(
      SEXP land_r,
      int nrow,
      int ncol,
      List disturbance_matrices,  // <-- Fixed: default to empty list
      int n_steps=20,
      int init_population=10, 
      int hr_size=1,
      double birth_rate=2.0,
      int breeding_age=1,
      double survival=0.4,
      double distance_weight=0.001,
      double dispersal_distance=5.0,
      int dispersal_mode=2,
      double sink_avoidance=0.5,
      double neigh_avoidance=1.0,
      double sink_mortality=0.7,
      const char* file_name="Res_imse")
{

// land_r is the input landscape that is obtained from R, which is a matrix
// the original matrix is transposed before being passed to Cpp
// so that the end landscape is equal to the original matrix
// the dimensions of the matrix (nrow and ncol) are also defined in R (before 
// transposing the original matrix)

  std::vector<double> landc = as< std::vector<double> >(land_r);

  double land[landc.size()];
  for(int i=0;i<landc.size();i++)
    land[i]=landc[i];
  

// Stores the simulation parameters in the object param
TSimParam param;

   param.land = new Mat_DP(land,nrow,ncol); 
      //stores the landscape passed from R
   param.initpopulation = init_population;
   param.nsteps = n_steps;
   param.hrsize = hr_size;
   param.birthrate = birth_rate;
   param.breedingage = breeding_age;
   param.survival = survival;
   param.distanceweight = distance_weight;
   param.dispersaldistance = dispersal_distance;
   param.dispersalmode = dispersal_mode;
   param.sinkavoidance = sink_avoidance;
   param.neighavoidance = neigh_avoidance;
   param.sinkmortality = sink_mortality;
   param.disturbance_matrices = disturbance_matrices;  
   param.filename = file_name;
      // name the file that is created with the results
      // of the simulation when indmodel.se is run

   TSimulator simulator(param);  // creates and starts simulation

   List popsizehist(n_steps + 1);
      // vector which stores population sizes at each time step
   popsizehist[0]=simulator.GetPopulationSize(); // stores first value of population size
   
   

   for (int i=1; i<=n_steps; i++)
   {
      simulator.Step();  // executes a step of the simulation
      popsizehist[i]=simulator.GetPopulationSize(); // stores the value of population size
   }

  IntegerVector ts = seq_len(n_steps+1);
      // creates a vector with the time steps

  List output = List::create(_["popsize"]= NumericVector(popsizehist.begin(), 
                                                         popsizehist.end()), 
                                                         _["n_steps"]= ts);

  return (output); 
       
   }
