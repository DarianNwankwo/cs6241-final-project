//
//  main.cpp
//  cs6241 Zika HJB Solver
//
//  Created by Mallory Gaspard on 4/22/21.
//

#include <iostream>
// Libraries
#include <iostream>
#include <fstream>
#include <cmath>
//#include <mpi.h>
//Project-specific Header files
#include "boost/multi_array.hpp"
#include "helper.h"
#include "SolverFunctions.hpp"

//Declaring parameter values
double gAlpha = 0.5;
double gK = 1000; //carrying capacity
double gRate = 0.15; // r
double gAlleleEffect = 2; //b - strength of the allele effect
double gDelta = 0.0025; // constant death rate
double gBitingPenalty = 2; 
//Grid formation - Note: All of these are initial values to test / placeholders
int gRegularMaxNodes = 70; //max node number of regular mosquitoes
int gFemaleWMaxNodes = 70; //max node number of female w-type
int gMaleWMaxNodes = 70; //max node number of male w-type
int gNumMosquitoesPerNode = 10;

double gTerminalT = 100; //terminal time
int gNt = 100; //max number of timeslices taken
//double gDT = gTerminalT / gNt; //dt
double gDT = 1; //dt
double gTau = gDT; 
double gH = 1;

double gAWMMax = 4; //max males that can be released
double gAWFMax = 2;
double gInfinity = gNumMosquitoesPerNode * 1000;
double gControlInfty = 20;

double gFemaleWControls[2] = {0, gAWFMax}; //array to hold control values tested for female w type
double gMaleWControls[2] = {0, gAWMMax}; //array to hold control values tested for male w type

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Running Zika Solver:\n";
    
    zikaHJBSolver();
    
    std::cout << "Finished running Zika Solver:\n";
    return 0;
}
