#include <fstream>
#include <vector>
#include <numeric>
#include <deque>
#include "landscape.h"
#include "simulator.h"
#include "individual.h"
#include <iterator>


// TCell == operator overload: tests whether two cells are equal

bool operator==(const TCell& c1, const TCell& c2)
{
 return ((c1.x==c2.x)&&(c1.y==c2.y));
}


// Tcell < operator overload: a cell is considered smaller than other if it has a smaller x value
// in case equal x values, than the cell is considered smaller if it has a smaller y value

bool operator<(const TCell& c1, const TCell& c2)
{
 if (c1.x!=c2.x)
   return (c1.x < c2.x);
 else
   return (c1.y < c2.y);
}


// TCell << operator overload: writes a cell in a stream or file in a Mathematica format

ostream& operator<<(ostream& os, const TCell& c)
{
 return os << '{' << c.x << ',' << c.y << '}';
}

// THomeRange << operator overload: writes a home range in a stream or a file in a Mathematica format

ostream& operator<<(ostream& os, const THomeRange& hr)
{
 os << '{';
 std::ostream_iterator<TCell,char> oo(os,",");
 copy(hr.begin(),--hr.end(),oo);
 os<<hr.back();
 os << '}';
 return os;
}

// TLandscape << operator overload: writes the landscape in a stream or a file in a Mathematica format

ostream& operator<<(ostream& s, const TLandscape& land)
{
 const Mat_DP& mat = land.GetLandscapeMatrix();

 s << "{";
 for (int i=0; i<mat.nrows(); i++)
   {
   s << "{";
   for (int j=0; j<mat.ncols(); j++)
      {
      s << mat[i][j];
      if (j!=mat.ncols()-1) s << ", ";
      }
   s << "}";
   if (i!=mat.nrows()-1) s << ",\n ";
   }
 s << "}";

 return s;
}


// TLandscape constructor(it is run when the object is first created): stores the parameters of the landscape
// Takes as input a matrix of real numbers as the landscape

TLandscape::TLandscape(TSimulator* simulatorIn, Mat_DP* land)
{
 xmax = land->nrows();
 ymax = land->ncols();
 mland = *land;
 mfree = mland;
 simulator = simulatorIn;
}

// Alternative TLandscape constructor(it is run when the object is first created): stores the parameters of the landscape
// and creates a uniform landscape where all cells have habaffty value.
// Takes as input the size of the landscape and the habitat affinity value.

TLandscape::TLandscape(TSimulator* simulatorIn, double habaffty,
                       int xmaxIn, int ymaxIn)
{
 xmax = xmaxIn;
 ymax = ymaxIn;
 mland = Mat_DP(habaffty,xmax,ymax);
 mfree = mland;
 simulator = simulatorIn;
}

// TLandscape destructor (it is run when the object is eliminated)
TLandscape::~TLandscape()
{
}


// Update: updates the matrix of free cells (mfree) by assigning -1 to all ocuppied cells

void TLandscape::Update(const TPopulation& population)
{
  mfree = mland;
  for(TPopulation::const_iterator individual=population.begin(); individual!=population.end(); individual++)
	  { //for all individuals of the population
		const THomeRange& homerange = individual->GetHomeRange();
		for (THomeRange::const_iterator i=homerange.begin(); i!=homerange.end();i++) // for all cells of the home range
			mfree[i->x][i->y]=-1;  // assign -1 to the occupied cell
	  }
}


// ChooseStartingPoint: choses the first cell of the home range based on the mothercell
// Stores the first cell in start cell and returns true if successful

bool TLandscape::ChooseStartingPoint(TCell& startcell, TCell& mothercell)
{
    if (simulator->GetStep()==1) //Global dispersal for initial population
        return ChooseStartingPointMode0(startcell);
    switch(simulator->GetDispersalMode())
    {
        case 0: return ChooseStartingPointMode0(startcell);  // Global dispersal
        case 1: return ChooseStartingPointMode1(startcell, mothercell); // Local dispersal, habitat search in a local kernel
        case 2: return ChooseStartingPointMode2(startcell, mothercell); // Local dispersal, random walk
    }
    return true;
}

// Global dispersal mode
// Chooses starting point for the home range based on global dispersal
bool TLandscape::ChooseStartingPointMode0(TCell& startcell)
{
 double maxaffty = mfree.max();  // calculates the maximum available affinity in the matrix
 if (maxaffty<0)                 // if matrix if full it is not possible to choose start point
   return false;

 // creates vector to store available cells in the lanscape
 vector<TCell> vlandmax(mland.nrows() * mland.ncols());

 int ncells=0;

 // identifies available cells in the landscape with maximum affinity
 for (int i=0; i<xmax; i++)
   for (int j=0; j<ymax; j++)
      if (mfree[i][j]==maxaffty)
         {
         vlandmax[ncells]=TCell(i,j);
         ncells++;
         }

 // generates a random number between 0 and ncells - 1
 int start = simulator->sto->IRandom(0,ncells-1);
 // selects a random cell from the vector of available cells with maximum affinity
 startcell = vlandmax[start];
 return true;
}

// Local dispersal mode: local habitat choice in a kernel (NEW: SELECTS BEST AVAILABLE HABITAT IN KERNEL)
// Chooses starting point for the home range based on local dispersal from mother cell
bool TLandscape::ChooseStartingPointMode1(TCell& startcell,
                                          TCell& mothercell)
{
  int r = simulator->GetDispersalDistance();
  int rsq = SQR(r);                 // stores r^2 in variable rsq
  
  vector<TCell> vlandmax(4 * rsq);  // Creates vector to store dispersal kernel
  int ncells = 0;
  
  double local_maxaffty = -1; // Initialize local maximum affinity
  
  for (int i = MAX(mothercell.x - r, 0); i < MIN(mothercell.x + r + 1, xmax); i++)
    for (int j = MAX(mothercell.y - r, 0); j < MIN(mothercell.y + r + 1, ymax); j++)
      if (SQR(i - mothercell.x) + SQR(j - mothercell.y) <= rsq) // Check if cell is in circle
        local_maxaffty = MAX(local_maxaffty, mfree[i][j]); // Update local max affinity
      
      if (local_maxaffty < 0) // If no cells have affinity > -1
        return false;
      
      for (int i = MAX(mothercell.x - r, 0); i < MIN(mothercell.x + r + 1, xmax); i++)
        for (int j = MAX(mothercell.y - r, 0); j < MIN(mothercell.y + r + 1, ymax); j++)
          if (mfree[i][j] == local_maxaffty && // Check if the cell matches local max affinity
              SQR(i - mothercell.x) + SQR(j - mothercell.y) <= rsq) // and lies in the circle
          {
            vlandmax[ncells] = TCell(i, j); // Add the cell to the dispersal kernel
            ncells++;
          }
          
          if (ncells > 0)
          {
            int start = simulator->sto->IRandom(0, ncells - 1); // Randomly select a cell
            startcell = vlandmax[start]; // Set the chosen cell as the start cell
            return true;
          }
          else
            return false;
}
// 
// Local dispersal mode: random walk
// Chooses starting point for the home range based on random-walk from mother cell
bool TLandscape::ChooseStartingPointMode2(TCell& startcell,
                                          TCell& mothercell)
{
    int r=simulator->GetDispersalDistance();
    
    TCell cell=mothercell; //current cell in dispersal
    
    for (int walkstep=0; walkstep<r; walkstep++)   // at each dispersal step
    {
        TCell newcell;
        //creates vector with neigbor cells of current cell
        std::vector<TCell> neigh;
        if (cell.x > 0) neigh.push_back(TCell(cell.x-1,cell.y));
        if (cell.y > 0) neigh.push_back(TCell(cell.x,cell.y-1));
        if (cell.x < xmax-1) neigh.push_back(TCell(cell.x+1,cell.y));
        if (cell.y < ymax-1) neigh.push_back(TCell(cell.x,cell.y+1));
        
        //creates vector with cumulative probabilities for each neighbor cell
        std::vector<double> neighprob;
        double cumprob=0;
        for (vector<TCell>::iterator i = neigh.begin(); i!=neigh.end(); i++)
         {
             double prob=1.0;   // default probability of dispersing to a neighbour
             if ((mland[i->x][i->y]==0)&&(mfree[i->x][i->y]<0)) // if sink habitat and occupied
                 prob=1-MAX(simulator->GetSinkAvoidance(),simulator->GetNeighAvoidance());  //decreases prob. by maximum of sinkavoidance and neighavoidance
             else if (mland[i->x][i->y]==0) // if sink habitat only
                 prob=1-simulator->GetSinkAvoidance();
             else if (mfree[i->x][i->y]<0) // if occupied only
                 prob=1-simulator->GetNeighAvoidance();
             cumprob+=prob;
             neighprob.push_back(cumprob);
         }
        
        // chooses a cell with its probabilty (multinomial) and updates the current cell
		if (cumprob==0) // if all probabilities are zero then cannot find starting point
		    return false;
		else {
			double urand=simulator->sto->Random();
			for (int k=neigh.size()-1; k>=0; k--)
				if (urand <= neighprob[k]/cumprob)
					newcell=neigh[k];
		}
		
		// if in a sink cell apply sink dispersal mortality and move again
        if ((mland[newcell.x][newcell.y]==0)&&(mland[cell.x][cell.y]!=0))
        {
			if (simulator->sto->Random()<=simulator->GetSinkMortality())
				return false;
            
            int dirx = newcell.x-cell.x;
            int diry = newcell.y-cell.y;
            if (dirx>0)
            {
                if (newcell.x < xmax-1) newcell.x++;
                else newcell.x--;
            }
            if (dirx<0)
            {
                if (newcell.x > 0) newcell.x--;
                else newcell.x++;
            }
            if (diry>0)
            {
                if (newcell.y < ymax-1) newcell.y++;
                else newcell.y--;
            }
            if (diry<0)
            {
                if (newcell.y > 0) newcell.y--;
                else newcell.y++;
            }
            
        }
        cell=newcell;
    }
    startcell=cell;

    if (mfree[cell.x][cell.y]<0)  // if cell occupied
	//if ((mfree[cell.x][cell.y]<0)||(mland[cell.x][cell.y]==0)) //if cell occupied or sink
        return false;                                           // dispersal unsuccessul
    else return true;
}


// PlaceHomeRange: places the home range in the landscape starting the dispersal in the mother cell

bool TLandscape::PlaceHomeRange(THomeRange& homerange,
                                TCell& hrcentermother)
{
 TCell start;

 while (ChooseStartingPoint(start, hrcentermother))   // while it is possible to find a starting cell for the HR expansion
    {
    mfree[start.x][start.y]=-1;
    homerange.push_back(start);                       // Stores the starting cell in the home range
    if (ExpandHomeRange(homerange))                   // If it is possible to exand home range to its desired size
      return true;                                    // return succes
    else homerange.clear();                           // else fails - it seems we also need to release back the cells to mfree
                                                      // and it is not being done, possible bug for local dispersal with HR>1
    }

 return false;    // if it was not possible to find a start cell fails
}

//---------------------------------------------------------------------------

bool TLandscape::ExpandHomeRange(THomeRange& homerange)
{
 TNeighbors neighbors;

 while (homerange.size() < simulator->GetHomeRangeSize())
   {
   CalculateNeighbors(homerange, neighbors);
   if (neighbors.empty())
      return false;
   TCell pt = ChoosePoint(homerange, neighbors);
   mfree[pt.x][pt.y]=-1;
   homerange.push_back(pt);
   }
 return true;
}

//---------------------------------------------------------------------------

void TLandscape::CalculateNeighbors(THomeRange& homerange,
                                    TNeighbors& neighbors)
{
 int neighdiff[8][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1},
                        {0, 1}, {1, -1}, {1, 0}, {1, 1}};

 const TCell& last = homerange.back();
 TNeighbors newneigh;
 for (int i=0; i<8; i++)
   {
   int x = last.x + neighdiff[i][0];
   int y = last.y + neighdiff[i][1];

   if ((x < xmax) && (x >= 0) &&
       (y < ymax) && (y >= 0))
      if (mfree[x][y] >= 0)
         newneigh.push_back(TCell(x,y));
   }
 newneigh.sort();
 neighbors.merge(newneigh);
 neighbors.unique();
 homerange.sort();
 set_difference(neighbors.begin(),neighbors.end(),
                homerange.begin(),homerange.end(),
                neighbors.begin());
}

//---------------------------------------------------------------------------

TCell TLandscape::ChoosePoint(const THomeRange& homerange,
                              TNeighbors& neighbors)
{
 TCell ctr = HomeRangeCenter(homerange);

 TNeighbors::iterator i;
 double max=-10000.;
 deque<TNeighbors::iterator> bestneigh;
 double val;
 for (i=neighbors.begin(); i!=neighbors.end(); ++i)
   {
   val = EvaluatePoint(*i,ctr);
   if (val>max)
      {
      max = val;
      bestneigh.clear();
      bestneigh.push_back(i);
      }
   else if (val==max)
      bestneigh.push_back(i);
   }
 int rndneigh = simulator->sto->IRandom(0,bestneigh.size()-1);
 neighbors.erase(bestneigh[rndneigh]);
 return *(bestneigh[rndneigh]);
}

//---------------------------------------------------------------------------

TCell TLandscape::HomeRangeCenter (const THomeRange& homerange)
{
 THomeRange::const_iterator i;
 double x=0;
 double y=0;
 double puse=0;

 if (homerange.size()==1)
	 return *homerange.begin();
	
 for (i=homerange.begin(); i!=homerange.end(); ++i)
   {
   puse += mland[i->x][i->y];
   x += i->x * mland[i->x][i->y];
   y += i->y * mland[i->x][i->y];
   }
 if (homerange.empty())
   return TCell(0,0);
 return TCell(round(x/puse),round(y/puse));
}

//---------------------------------------------------------------------------
double TLandscape::EvaluatePoint (const TCell& pt, const TCell& ctr)
{
 double distw = simulator->GetDistanceWeight() * Distance(pt, ctr);

 return mland[pt.x][pt.y]*(1-distw);
}

//---------------------------------------------------------------------------
double TLandscape::CalculateOptimalFitness ()
{
 TCell start(xmax/2,ymax/2);
 THomeRange homerange;

 mfree[start.x][start.y]=-1;
 homerange.push_back(start);
 ExpandHomeRange(homerange);
 TCell hrcenter = HomeRangeCenter(homerange);
 double d = 0;
 for (THomeRange::iterator i=homerange.begin();
      i!=homerange.end(); i++)
         d+=EvaluatePoint(*i,hrcenter);
 return d;
}

//---------------------------------------------------------------------------
double Distance(const TCell& c1, const TCell& c2)
{
 return SQR(c1.x-c2.x)+SQR(c1.y-c2.y);
}

