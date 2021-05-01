//
//  helper.h
//  cs6241 Zika HJB Solver
//
//  Created by Mallory Gaspard on 4/22/21.
//

//Libraries
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <tuple>
#include <utility>
//Header Files
#include "boost/multi_array.hpp"

#ifndef helper_h
#define helper_h

//Model parameters
extern double gAlpha; //proportion of population that are female. 
extern double gK; //carrying capacity
extern double gRate; // r
extern double gAlleleEffect; //b - strength of the allele effect
extern double gDelta; // constant death rate
extern double gBitingPenalty; 
//Grid formation
extern int gRegularMaxNodes; //max number of regular mosquitoes
extern int gFemaleWMaxNodes; //max number of female w-type
extern int gMaleWMaxNodes; //max number of male w-type
extern int gNumMosquitoesPerNode; //number of mosquitoes per grid node 
extern int gNt; //max number of timeslices taken
extern double gTerminalT; //terminal time
extern double gDT; //dt
extern double gH;
extern double gTau; 
extern double gTau; //semi-lagrangian step

//Control specifications
extern double gAWMMax; //max males that can be released
extern double gAWFMax;
extern double gInfinity;
extern double gControlInfty;

extern double gFemaleWControls[]; //array to hold control values tested for female w type
extern double gMaleWControls[]; //array to hold control values tested for male w type

#endif /* helper_h */
