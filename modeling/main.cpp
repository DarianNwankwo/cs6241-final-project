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
double gAlpha = 0.75;
double gK = 10000; //carrying capacity
double gRate = 0.5; // r
double gAlleleEffect = 2; //b - strength of the allele effect
double gDelta = 100; // constant death rate

//Grid formation - Note: All of these are initial values to test / placeholders
int gRegularMaxNodes = 50; //max node number of regular mosquitoes
int gFemaleWMaxNodes = 50; //max node number of female w-type
int gMaleWMaxNodes = 50; //max node number of male w-type
int gNumMosquitoesPerNode = 100;

double gTerminalT = 1; //terminal time
int gNt = 100; //max number of timeslices taken
double gDT = gTerminalT / gNt; //dt

double gAWMMax = 10; //max males that can be released
double gAWFMax = 10;
double gInfinity = 1e7;

double gFemaleWControls[2] = {0, gAWFMax}; //array to hold control values tested for female w type
double gMaleWControls[2] = {0, gAWMMax}; //array to hold control values tested for male w type

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Running Zika Solver:\n";
    
    zikaHJBSolver();
    
    std::cout << "Finished running Zika Solver:\n";
    return 0;
}
