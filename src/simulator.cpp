#include <cstdlib>
#include <numeric>
#include <fstream>
#include <functional>
#include "simulator.h"


#include "randomc/mersenne.cpp"         // members of class TRandomMersenne
#define STOC_BASE CRandomMersenne
#include "stocc/stoc1.cpp"            // random library source code
#include "randomc/userintf.cpp"
#include <sys/time.h>


// Constructor of TSimulator (it is run when the object is first created): initializes and starts the simulation

TSimulator::TSimulator(const TSimParam& param)
{
 struct timeval time; 
 gettimeofday(&time,NULL);
	
 // microsecond has 1 000 000
 // Assuming you did not need quite that accuracy
 // Also do not assume the system clock has that accuracy.
 // srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
 // The trouble here is that the seed will repeat every
 // 24 days or so.
	
 // If you use 100 (rather than 1000) the seed repeats every 248 days.
	
 // Do not make the MISTAKE of using just the tv_usec
 // This will mean your seed repeats every second.
	
 long int seed=(time.tv_sec * 100) + (time.tv_usec / 100);	
 //long int seed = time(0);      // random seed
 sto=new StochasticLib1(seed); // make instance of random library

    
 //stores the parameter values of the simulation in local variables to the object
 nsteps=param.nsteps;
 initpopulation=param.initpopulation;
 hrsize=param.hrsize;
 birthrate=param.birthrate;
 breedingage=param.breedingage;
 survival=param.survival;
 distanceweight=param.distanceweight;
 dispersaldistance=param.dispersaldistance;
 dispersalmode=param.dispersalmode;
 neighavoidance=param.neighavoidance;
 sinkavoidance=param.sinkavoidance;
 sinkmortality=param.sinkmortality;
    
 filename=param.filename;

 // the first step of the simulation is run here so the step counter is set to 1
 step=1;
    
 // creates a landscape object based on the input landscape in param.land
 landscape = new TLandscape(this,param.land);

 // writes in the output file (with name filename) the value of the parameters
 OutputParameters();

 // creates an optimal landscape with the same of the size of the landscape of the simulation
 // in the optimal landscpae all cell have habitat affinity 1
 TLandscape landscapeopt(this,1,param.land->nrows(),param.land->ncols());
 optimalfitness = landscapeopt.CalculateOptimalFitness();
 // optimalfitness is the fitness of an individual with the best possible home range in an empty and uniform landscape
 // it corresponds to the normalizing value Phi in the paper

 // starts at the center
 TCell mothercell(param.land->nrows()/2,param.land->ncols()/2);
 
 // creates the initial population of individuals (ninitpopulation objects of the class TIndividual) and stores them in the
 // population list of individuals
 for (int n=0; n<initpopulation; n++)
   population.push_back(TIndividual(this,mothercell));
 
 // setlles the home range of each individual in the initial population
 for_each(population.begin(),population.end(),
          mem_fun_ref(&TIndividual::SettleHomeRange));
 
//kill floaters (individuals with empty home ranges)
 population.erase(remove_if(population.begin(),population.end(),
                            mem_fun_ref(&TIndividual::HasEmptyHomeRange)),
                  population.end());
 
 //writes in the output file the list of the cells of each individual and the ages of each individual
 OutputGeneration();
}


// Step: executes one step of the simulation

void TSimulator::Step()
{
 TPopulation popjuv;      //list of juveniles

 step++;                  //increases step counter
    
 //increases ages for each individual and produces juveniles (stores in popjuv)
 for (TPopulation::iterator i = population.begin(); i!=population.end(); i++)
      i->ApplyBreeding(popjuv);
    
 /*cout << "j:" << popjuv.size() << ' ';*/
    
 //kill adults randomly and remove them from the population (adult mortality)
 population.erase(remove_if(population.begin(),population.end(),
                            mem_fun_ref(&TIndividual::ApplyMortality)),
                  population.end());
    
 landscape->Update(population);  //actualize matrix of free cells opened by adult mortality
    
    
 //kill juveniles randomly and remove them from the population (first stage of juvenile mortality)
 popjuv.erase(remove_if(popjuv.begin(),popjuv.end(),
                            mem_fun_ref(&TIndividual::ApplyMortality)),
                  popjuv.end());
 
 //settle the home-range of each juvenile
 for_each(popjuv.begin(),popjuv.end(),
          mem_fun_ref(&TIndividual::SettleHomeRange));
 
 //kill juveniles without home-range (floaters)
 popjuv.erase(remove_if(popjuv.begin(),popjuv.end(),
                            mem_fun_ref(&TIndividual::HasEmptyHomeRange)),
                  popjuv.end());
 
 //insert juveniles into the adult population
 population.insert(population.end(), popjuv.begin(), popjuv.end());
 
 
 // apply spatially autocorrelated mortality
 population.erase(
   remove_if(population.begin(), population.end(),
             [this](TIndividual& ind) { return ind.ApplySpatialMortality(this->step); }),
             population.end());
 
 
 
 landscape->Update(population);  //actualize matrix of free cells opened by spatial mortality
 
 
 /*cout << population.size() << ' '; */
 
 //writes in the output file the list of the cells of each individual and the ages of each individual
 OutputGeneration();
}


// OutputParameters: writes the parameters of the simulation in a file

void TSimulator::OutputParameters()
{
    // opens a file with name filename for writing
    ofstream os(filename.c_str());
    
    // writes each parameter into the file in a format readable by Mathematica
    os << "simoptions = \n{";
    os << "HomeRangeSize -> " << hrsize << ", ";
    os << "Fecundity -> " << birthrate << ", ";
    os << "BreedingAge -> " << breedingage << ", ";
    os << "Survival -> " << survival << ", ";
    os << "InitialPopulation -> " << initpopulation << ", ";
    os << "DistanceWeight -> " << distanceweight << ", ";
    os << "DispersalDistance -> " << dispersaldistance << ", ";
    os << "DispersalMode -> " << dispersalmode;
    os << "};\n";
    os << "landscape = \n" << *landscape << ";\n";
    // landscape is a matrix with the habitat; each cell has the value of the habitat affinity for that cell
    os << "hrmaphist = ageshist = Table[Null,{" << (nsteps+1) << "}];\n";
    os << "popsize = Table[Null,{" << (nsteps+1) << "}];\n";
    // hrmaphist(ory) is a list of list of the cells of each individual at each step of the simulation
    // agehist is a list of the individuals ages at each step of the simulation
}

// OutputGeneration: writes the individuals alive at a given simulation step into a file

void TSimulator::OutputGeneration()
{
 // opens a file with name filename for writing
 ofstream os(filename.c_str(),ios_base::app);
    
 // writes a list of the list of home range cells of each individual in a Mathematica format
 os << "hrmaphist[[" << step << "]]=\n{";

 for (TPopulation::iterator i=population.begin();
      i!=population.end(); i++)
    {
    i->OutputHomeRange(os);
    if (i != --population.end())
      os << ",\n";
    }
 os << "};\n";

    
 // writes a list of the individual ages in a Mathematica format   
 os << "ageshist[[" << step << "]]=\n{";
 for (TPopulation::iterator i=population.begin();
      i!=population.end(); i++)
   {
   os << i->GetAge();
   if (i != --population.end())
     os << ", ";
   }
 os << "};\n";
    
 os << "popsize[[" << step << "]]=\n";
 os << population.size() << ";\n";
}


// Destructor of TSimulator (it is run when the object is destroyed): releases allocated memory

TSimulator::~TSimulator()
{
 delete landscape;
}

