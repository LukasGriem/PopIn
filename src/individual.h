#ifndef _INDIVIDUAL_H_
#define _INDIVIDUAL_H_

#include "landscape.h"

class TSimulator;

class TIndividual
{
 public:
         TIndividual(TSimulator*, TCell&);
         THomeRange& GetHomeRange() {return homerange;}
         const THomeRange& GetHomeRange() const {return homerange;}
         unsigned int GetAge() {return age;}
         bool ApplyMortality();
         bool ApplySpatialMortality(int step);
         void ApplyBreeding(TPopulation& popjuv);
         bool HasEmptyHomeRange() {return homerange.empty();}
         void OutputHomeRange(ostream&);
         void SettleHomeRange();
 private:
         unsigned int age;
         int fledglings;
         THomeRange homerange;
         TCell hrcenter;
         TCell hrcentermother;
         TSimulator* simulator;
         void CalculateFledglings();
};

int iround(double x);

#endif
