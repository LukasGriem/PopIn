
#include "individual.h"
#include "simulator.h"



// Constructor of TIndividual: it is run when an individual is created

TIndividual::TIndividual(TSimulator* simulatorIn, TCell& hrcentermotherIn)
{
 // store the parameter of the constructor function in the object
 simulator = simulatorIn;
 hrcentermother = hrcentermotherIn;
 
 // the age and the number of offspring of the individual are initialized to zero
 age = 0;
 fledglings = 0;
}


// SettleHomeRange: Selects a home-range in the landscape for an individual

void TIndividual::SettleHomeRange()
{
 TLandscape* landscape = simulator->GetLandscape();
    
 if (!landscape->PlaceHomeRange(homerange,hrcentermother))   // tries to setlle a home-range
   homerange.clear();                                   // if not sucessful clears home-range
 else hrcenter = landscape->HomeRangeCenter(homerange); //else calculates the center of the HR
	
}


// OutputHomeRange: Writes the homerange cells in a file

void TIndividual::OutputHomeRange(ostream& os)
{
 //#define os cout
 os << homerange;
}


// ApplyMortality: Determines whether the individual dies, either stochastically or deterministically
// Returns true if the individual dies and false if the individual survives

bool TIndividual::ApplyMortality()
{
 if (simulator->GetSurvival()>1.0) //Deterministic simulation: deterministic death
   {
   if(age < simulator->GetSurvival()) //if age smaller than maximum age
     return false;                    //then individual survives
   return true;                       //otherwise dies
   }
 else                               //Stochastic simulation: stochastic death
   {                                //if Bernouli event with probability simulator->GetSurvival() is true
   if (simulator->sto->Bernoulli(simulator->GetSurvival()))
      return false;                 //then individual survives
    return true;                    //otherwise dies
   }
}

// ApplySpatialMortality: Kills all individuals at specific locations
// Returns true if the individual dies and false if the individual survives
bool TIndividual::ApplySpatialMortality(int step) {
  if (!simulator) return false;
  
  Rcpp::List disturbance_matrices = simulator->GetDisturbanceMatrices(); // retrieve matrix
  
  if (step >= disturbance_matrices.size()) return false;  // Prevent out-of-bounds error
  Rcpp::NumericMatrix disturbance_matrix = disturbance_matrices[step];
  
  int x = hrcenter.x;
  int y = hrcenter.y;
  
  if (x < 0 || x >= disturbance_matrix.nrow() || y < 0 || y >= disturbance_matrix.ncol()) {
    return false;  // Out of bounds, individual survives
  }
  
  double extinction_prob = disturbance_matrix(x, y);
  return (R::runif(0, 1) < extinction_prob);
}





// ApplyBreeding: Produces the fledglings of an individual in a given year

void TIndividual::ApplyBreeding(TPopulation& popjuv)
{
 age++;  // Increase the age of the individual (a reproductive season has happened)
 CalculateFledglings();  // Calculates the number of fledglings based on the home range

       // Store the fledglings in the list popjuv
 for (int n=0; n<fledglings; n++)
   popjuv.push_back(TIndividual(simulator,hrcenter));
}


// CalculateFledglings: calculates the number of fledgings based on the home range

void TIndividual::CalculateFledglings()
{
 if (age >= simulator->GetBreedingAge())  //if age is greater than breeding age
   {
   TLandscape* land = simulator->GetLandscape();
   double d = 0;

   for (THomeRange::iterator i=homerange.begin();  //sums the fitness (energy yield) over all the cells of the home range
        i!=homerange.end(); i++)
     d+=land->EvaluatePoint(*i,hrcenter);
   // d = simulator->GetBirthRate();                // even if an individual ends up on a sink habitat at the end of the dispersal steps, it can still reproduce
	   
   d *= simulator->GetBirthRate();           // multiplies fitness by fecundity (b0 in the model)
												  // this was active in the initial version
   d /= simulator->GetOptimalFitness();           // normalizes by the optimal home-range fitness (Phi in the model)
   if (simulator->GetSurvival()>=1.0)             //Deterministic simulation
      fledglings = iround(d);                     //The number of fledglings equals the fecundity
   else
      fledglings = simulator->sto->Poisson(d);    //The number of fledglings is a Poisson with mean equal to fecundity
   }
 else fledglings = 0;    //if individual has not reached breeding age
}


// iround: Helper function that rounds a real number to the nearest integer
int iround(double x)
{
    double dum;
    if (fabs(modf(x,&dum))==0.5)
    {
        if (int(floor(x))%2==0)
            return floor(x);
        else return ceil(x);
    }
    return floor(x+.5);
}

 
