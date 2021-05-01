//
//  SolverFunctions.hpp
//  cs6241 Zika HJB Solver
//
//  Created by Mallory Gaspard on 4/22/21.
//

#ifndef SolverFunctions_hpp
#define SolverFunctions_hpp

#include <stdio.h>
//#include "boost/multi_array.hpp"
#include "helper.h"

typedef boost::multi_array<double, 4> multiarray;

double runningCost (double aAWF, double aAWM, double aR, double aTau);

void optimalValue (multiarray *aValueFunction, multiarray *aFemaleWControls, multiarray *aMaleWControls, int aGridNumR, int aGridNumFW, int aGridNumMW, int aTimeIndex);

void initializeArray (multiarray *aArray, const double aIllegalValue);

void writeToFile(multiarray *aValueFunction, multiarray *aFemaleWControls, multiarray *aMaleWControls, int aSelectedTimeslice);

void zikaHJBSolver();

#endif /* SolverFunctions_hpp */
