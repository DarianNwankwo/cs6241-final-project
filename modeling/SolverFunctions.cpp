//
//  SolverFunctions.cpp
//  cs6241 Zika HJB Solver
//
//  Created by Mallory Gaspard on 4/22/21.
//
#include <iostream>
#include <fstream>
#include <cmath>
#include "SolverFunctions.hpp"
#include "helper.h"
#include "boost/multi_array.hpp"

using namespace std;

typedef boost::multi_array<double, 4> multiarray;

//Running cost
double runningCost (double aAWF, double aAWM, double aR, double aTau){
    double rCostVal;
    
    //compute integral approximation
    rCostVal = (aAWF + aAWM + aR * gAlpha * gBitingPenalty) * aTau; 
    
    return rCostVal;
}

//SL interpolation scheme 
void optimalValue (multiarray *aValueFunction, multiarray *aFemaleWControls, multiarray *aMaleWControls, int aGridNumR, int aGridNumFW, int aGridNumMW, int aTimeIndex){
    //takes in number of R, FW, and MW type mosquitoes, carries out trilinear interpolation with four corners of the grid cell to determine value function at
    
    //Plan: Bilinear on upper face and lower face of interp cube, then linear interpolation between the two faces to get value at the point in the cube.
    
    //Mallory Note: For now, I'd prefer to NOT do interpolation in time as well. With that said, we need to think about the grid / timescales to make sure they make sense. Also, be careful with typecasting!
    
    //Current number of mosquitoes of each type at the current  timestep.
    double currentNumR = gNumMosquitoesPerNode * aGridNumR;
    double currentNumFW = gNumMosquitoesPerNode * aGridNumFW;
    double currentNumMW = gNumMosquitoesPerNode * aGridNumMW;
    double currentNumF = currentNumR * gAlpha + currentNumFW;
    double currentNumM = currentNumR * (1 - gAlpha) + currentNumMW;
    
    //Initialize placeholders for optimal controls and value function values
    double optimalAWF;
    double optimalAWM;
    double optimalValueFunction = gInfinity;
    
    //About the very artificial boundary conditions: Severely penalizes having the maximum amount of mosquitoes of any type. Easiest thing to do for now, ideal thing would be to carve out a boundary based on the maximum values of each type and what you'd need to ensure that the number of those mosquitoes decrease from there, but I would prefer not to deal with cut cells in my SL scheme right now 
    if(currentNumR == 0){
        //No R-type mosquitoes --> no Zika transmission --> zero value function
        optimalAWF = 0;
        optimalAWM = 0;
        optimalValueFunction = 0;
    }//end first bc if
    
    else if ((currentNumFW == gFemaleWMaxNodes * gNumMosquitoesPerNode) && (currentNumR != 0)){
        //FOR NOW: If at the max, not going to add anything of any type
        optimalAWF = gControlInfty;
        optimalAWM = gControlInfty;
        optimalValueFunction = gInfinity;
    }//end FW = maxFw

    else if ((currentNumMW == gMaleWMaxNodes * gNumMosquitoesPerNode) && (currentNumR != 0)){
        optimalAWF = gControlInfty;
        optimalAWM = gControlInfty;
        optimalValueFunction = gInfinity;
    }//end MW = maxFw
    
    else if (currentNumR == gRegularMaxNodes * gNumMosquitoesPerNode){
        optimalAWF = gControlInfty;
        optimalAWM = gControlInfty;
        optimalValueFunction = gInfinity;
    }//end R = Rmax
    
    /*else if ((currentNumMW == 0) && (currentNumR != 0) && (currentNumFW != gFemaleWMaxNodes * gNumMosquitoesPerNode) && (currentNumR != gNumMosquitoesPerNode * gRegularMaxNodes)){
        optimalAWF = gControlInfty;
        optimalAWM = gControlInfty;
        optimalValueFunction = gInfinity;
    }//end MW = 0
    
    else if ((currentNumFW == 0) && (currentNumR != 0)){
        optimalAWF = gControlInfty;
        optimalAWM = gControlInfty;
        optimalValueFunction = gInfinity;
    }//end FW = 0 */
    
    else {
    //LOOP OVER ADMISSIBLE CONTROL VALUES
        for (short int fwControlCandidateIndex = 0; fwControlCandidateIndex < 2; fwControlCandidateIndex++ ){
            
            double currentWFControl = gFemaleWControls[fwControlCandidateIndex];
            
            for (short int mwControlCandidateIndex = 0; mwControlCandidateIndex < 2; mwControlCandidateIndex++ ){
                
                int nextTimeIndex = aTimeIndex + 1;
                
                //get current controls to test
                double currentWMControl = gMaleWControls[mwControlCandidateIndex];
                //placeholder for current controls and vf value
                double currentValueFunction;
                
                //Euler step to advance mosquito populations one step forward in time
                double nextR = currentNumR + gDT * (gRate * (1 - ((currentNumF + currentNumM) / gK)) * ((gAlpha * currentNumR * (1 - gAlpha) * currentNumR) / (gAlleleEffect + currentNumM)) - gDelta * currentNumR);
                double nextFW = currentNumFW + gDT * (gRate * gAlpha *  (1 - ((currentNumF + currentNumM) / gK)) * (currentNumFW * (currentNumM / (gAlleleEffect + currentNumM))) - gDelta * currentNumFW + currentWFControl);
                double nextMW = currentNumMW + gDT * (gRate * (1 - gAlpha) * (1 - ((currentNumF + currentNumM) / gK)) * (currentNumFW * (currentNumM / (gAlleleEffect + currentNumM))) - gDelta * currentNumMW + currentWMControl);
                
                //check to see if any of the next values are out of bounds
                //todo: fix this! This is leading to incorrect boundary condititions
                if ((nextR > gRegularMaxNodes * gNumMosquitoesPerNode) || (nextFW > gFemaleWMaxNodes * gNumMosquitoesPerNode) || (nextMW > gMaleWMaxNodes * gNumMosquitoesPerNode)){
                    currentValueFunction = gInfinity;
                    optimalAWM = gControlInfty;
                    optimalAWF = gControlInfty;
                }//end check on bounds
                
                //put in clause to handle the zero-faces
                
                else if ((nextR < 0) || (nextFW < 0) || (nextMW < 0)){
                    currentValueFunction = gInfinity;
                    optimalAWM = gControlInfty;
                    optimalAWF = gControlInfty;
                }
                
                else{
                    //Step one: Convert the next values back to grid coordinates:
                    double nextRGrid = nextR / gNumMosquitoesPerNode;
                    double nextFWGrid = nextFW / gNumMosquitoesPerNode;
                    double nextMWGrid = nextMW / gNumMosquitoesPerNode;
                    
                    //Determine the bounding vertices in each direction
                    //Mallory Note: keep in mind rounding - be careful of rounding outside of current cell
                    int lowerR = int(floor(nextRGrid));
                    int upperR = int(ceil(nextRGrid));
                    int lowerFW = int(floor(nextFWGrid));
                    int upperFW = int(ceil(nextFWGrid));
                    int lowerMW = int(floor(nextMWGrid));
                    int upperMW = int(ceil(nextMWGrid));
                    
                    //Interpolation
                    //Upper Cube Face in FW - MW plane
                    
                    //UTL = upper top left, UBL = upper bottom left, UTR = upper top right, UBR = upper bottom right
                    
                    //pull these from arrays - be careful, may want / need to scale back to grid coords
                    double UTL = (*aValueFunction)[lowerFW][upperMW][upperR][nextTimeIndex];
                    double UBL = (*aValueFunction)[lowerFW][lowerMW][upperR][nextTimeIndex];
                    double UTR = (*aValueFunction)[upperFW][upperMW][upperR][nextTimeIndex];
                    double UBR = (*aValueFunction)[upperFW][lowerMW][upperR][nextTimeIndex];
                    
                    double betaUpper = nextMWGrid - lowerMW;
                    double gammaUpper = nextFWGrid - lowerFW;
                    
                    double Q2Upper = (1 - betaUpper) * UBL + betaUpper * UTL;
                    double Q1Upper = (1 - betaUpper) * UBR + betaUpper * UTR;
                    
                    double interpValueUpper = gammaUpper * Q1Upper + (1 - gammaUpper) * Q2Upper;
                    
                    //Bilinear interpolation on lower face
                    //pull these from arrays
                    double LTL = (*aValueFunction)[lowerFW][upperMW][lowerR][nextTimeIndex];
                    double LBL = (*aValueFunction)[lowerFW][lowerMW][lowerR][nextTimeIndex];
                    double LTR = (*aValueFunction)[upperFW][upperMW][lowerR][nextTimeIndex];
                    double LBR = (*aValueFunction)[upperFW][lowerMW][lowerR][nextTimeIndex];
                    
                    double betaLower = nextMWGrid - lowerMW;
                    double gammaLower = nextFWGrid - lowerFW;
                    
                    double Q2Lower = (1 - betaLower) * LBL + betaLower * LTL;
                    double Q1Lower = (1 - betaLower) * LBR + betaLower * LTR;
                    
                    double interpValueLower = gammaLower * Q1Lower + (1 - gammaLower) * Q2Lower;
                    
                    //Interpolation between lower and upper slices in R direction
                    double betaR = (nextRGrid - lowerR);
                    assert((betaR + (1-betaR)) == 1);
                    
                    double interpValue = betaR * interpValueUpper + (1 - betaR) * interpValueLower;
                    
                    //double convertedInterpValue = interpValue * gNumMosquitoesPerNode;
                    
                    assert(interpValue >= 0);
                    
                    currentValueFunction = interpValue + runningCost(currentWFControl, currentWMControl, currentNumR, gDT); //Note: assuming tau = gDT, hence why no interp in time now
                    //currentValueFunction = convertedInterpValue + runningCost(currentWFControl, currentWMControl, currentNumR, gDT); //Note: assuming tau = gDT, hence why no interp in time now
                    
                    
                }//end interpolation
                
                //check to see if this value is less than the current saved optimal.
                if (currentValueFunction < optimalValueFunction){
                    optimalValueFunction = currentValueFunction;
                    optimalAWF = currentWFControl;
                    optimalAWM = currentWMControl;
                }
                
            }//end inner control loop
        }//end outer control loop
    }//end else
    
    if (optimalValueFunction > 0){
        
    }
    
    //Save optimal control values and value function
    (*aValueFunction)[aGridNumFW][aGridNumMW][aGridNumR][aTimeIndex] = optimalValueFunction;
    (*aFemaleWControls)[aGridNumFW][aGridNumMW][aGridNumR][aTimeIndex] = optimalAWF;
    (*aMaleWControls)[aGridNumFW][aGridNumMW][aGridNumR][aTimeIndex] = optimalAWM;
    
}//end function

//dynamic programming loop
void dynamicProgLoop(multiarray *aValueFunction, multiarray *aFemaleWControls, multiarray *aMaleWControls){
    /*
     * Function: dynamicProgLoop
     *
     * Purpose: Solves the HJB equation backwards in time starting from the second
     * to last time index (Note: arrays aValueFunction and aOptimalControls are pre-populated
     * at the last timestep according to the boundary conditions before passing them to this function).
     */
    
    for(int t = gNt - 1; t >= 0; t--) {//time loop
        for(short int k = 1; k < gRegularMaxNodes + 1; k++) {//regular type loop
            for(short int i = 0; i < gFemaleWMaxNodes + 1; i++) {// fw loop
                for(short int j = 0; j < gMaleWMaxNodes + 1; j++) {// mw loop
                    optimalValue(aValueFunction, aFemaleWControls, aMaleWControls, k, i, j, t); //compute u(d,v,t)
              }//end mw loop
            } //end fw loop
        } //end r loop
        
        //Write results for this timeslice to file if desired
        if (t == 0){
            writeToFile(aValueFunction, aFemaleWControls, aMaleWControls, 0);
        }
        
    }//end time loop
}


void initializeArray (multiarray *aArray, const double aIllegalValue){
    /*
     * Function: initialzeArray()
     *
     * Purpose: Populates an array with values while imposing the proper
     * boundary conditions in the context of the problem.
     *
     */
    for(int t = gNt; t >= 0; t--) {//time loop
        for(short int k = 0; k < gRegularMaxNodes + 1; k++) {
            for(short int i = 0; i < gFemaleWMaxNodes + 1; i++) {// fw loop
                for(short int j = 0; j < gMaleWMaxNodes + 1; j++) {// mw loop
                    if (t == gNt){
                        (*aArray)[i][j][k][t] = 0; //check on this
                    }
                    //along the line R = 0, value function should be 0
                    if (k == 0){
                        (*aArray)[i][j][k][t] = 0;
                    }
                    
                    if ((k == 0) && (i ==0) && (j == 0)){
                        (*aArray)[i][j][k][t] = 0; //0 at the origin - no mosquitoes, no problem!
                    }
                    
              }//end mw loop
            } //end fw loop
        } //end r loop
    }//end time loop
}


/*==============================================================================
 * Write results to a file
 *============================================================================*/
void writeToFile(multiarray *aValueFunction, multiarray *aFemaleWControls, multiarray *aMaleWControls, int aSelectedTimeslice){

    ofstream ufile;
    ofstream femaleWControlsFile;
    ofstream maleWControlsFile;

    ufile.open ("value_function_output.txt");

    femaleWControlsFile.open("female_controls_output.txt"); //optimal control values
        
    maleWControlsFile.open("male_controls_output.txt"); //error values
        
   //ufile << gN << " " << gDMax << " " << gVMax << " " << gLightLength << " "<< gDNum << " " << gVNum <<endl;
   for(short int k = 0; k < gRegularMaxNodes + 1; k++) {           // time
        for(short int i = 0; i < gFemaleWMaxNodes + 1; i++) {      // distance
            for(short int j = 0; j < gMaleWMaxNodes + 1; j++) {  // velocity

                double valueFunction = (*aValueFunction)[i][j][k][aSelectedTimeslice];
                ufile << valueFunction;
                    
                double femaleWControl = (*aFemaleWControls)[i][j][k][aSelectedTimeslice];
                femaleWControlsFile << femaleWControl;
                
                double maleWControl = (*aMaleWControls)[i][j][k][aSelectedTimeslice];
                maleWControlsFile << maleWControl;
                
                if (j < gMaleWMaxNodes) {
                    ufile << " ";
                    femaleWControlsFile << " ";
                    maleWControlsFile << " ";
                }
            }
            ufile << endl;
            femaleWControlsFile << endl;
            maleWControlsFile << endl;
        }
    }
    
    ufile.close();
    femaleWControlsFile.close();
    maleWControlsFile.close();
}

//Solver function
void zikaHJBSolver(){
    //init
    multiarray u(boost::extents[gFemaleWMaxNodes+1][gMaleWMaxNodes+1][gRegularMaxNodes+1][gNt + 1]);
    multiarray controlsFW(boost::extents[gFemaleWMaxNodes+1][gMaleWMaxNodes+1][gRegularMaxNodes+1][gNt +1]);
    multiarray controlsMW(boost::extents[gFemaleWMaxNodes+1][gMaleWMaxNodes+1][gRegularMaxNodes+1][gNt + 1]);
    
    //initialize arrays at last step. Note that terminal cost is zero.
    initializeArray(&u, gInfinity);
    initializeArray(&controlsFW, gInfinity);
    initializeArray(&controlsMW, gInfinity);
    
    dynamicProgLoop(&u, &controlsFW, &controlsMW);
    
}
