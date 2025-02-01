#ifndef _LANDSCAPE_H_
#define _LANDSCAPE_H_

#include <list>

#include "nrtypes.h"

using namespace std;

class TSimulator;
class TIndividual;
typedef list<TIndividual> TPopulation;

struct TCell
{
   int x;
   int y;
   TCell(int xIn, int yIn): x(xIn), y(yIn) {};
   TCell() {};
};

bool operator==(const TCell& c1, const TCell& c2);
bool operator<(const TCell& c1, const TCell& c2);
ostream& operator<<(ostream& s, const TCell& c);

typedef list<TCell> THomeRange;
typedef list<TCell> TNeighbors;

ostream& operator<<(ostream& s, const THomeRange& hr);

class TLandscape
{
 public:
   TLandscape(TSimulator*, Mat_DP*);
   TLandscape(TSimulator*, double, int, int);
   ~TLandscape();
   Mat_DP& GetLandscapeMatrix() {return mland;}
   const Mat_DP& GetLandscapeMatrix() const {return mland;}
   bool PlaceHomeRange(THomeRange&, TCell&);
   void Update(const TPopulation&);
   double EvaluatePoint (const TCell& pt, const TCell& ctr);
   TCell HomeRangeCenter(const THomeRange&);
   double CalculateOptimalFitness();
   int xmax;
   int ymax;
 private:
   bool ChooseStartingPoint(TCell&, TCell&);
   bool ChooseStartingPointMode0(TCell&);
   bool ChooseStartingPointMode1(TCell&, TCell&);
   bool ChooseStartingPointMode2(TCell&, TCell&);
   bool ExpandHomeRange(THomeRange&);
   void CalculateNeighbors(THomeRange&, TNeighbors&);
   TCell ChoosePoint(const THomeRange& homerange, TNeighbors& neighbors);


   Mat_DP mland;        // a matrix of the landscape with each cell having an affinity value (between 0 and 1)
   Mat_DP mfree;        // a matrix where the occupied cells have value -1, only free cells have the correct affinity
   TSimulator* simulator;
};

ostream& operator<<(ostream& s, const TLandscape& land);

double Distance(const TCell& c1, const TCell& c2);

#endif
