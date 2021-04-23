//
//  SolverFunctions.cpp
//  cs6241 Zika HJB Solver
//
//  Created by Mallory Gaspard on 4/22/21.
//

#include "SolverFunctions.hpp"
#include "helper.h"
#include "boost/multi_array.hpp"

typedef boost::multi_array<double, 4> multiarray;

//Running cost
double runningCost (double aAWF, double aAWM, double aR, double aTau){
    double rCostVal;
    
    //compute integral approximation
    rCostVal = (aAWF + aAWM + aR * gAlpha) * aTau; 
    
    return rCostVal;
}

//SL interpolation scheme 
void optimalValue (multiarray *aValueFunction, multiarray *aFemaleWControls, multiarray *aMaleWControls, int aGridNumR, int aGridNumFW, int aGridNumMW, int aTimeIndex){
    //takes in number of R, FW, and MW type mosquitoes, carries out trilinear interpolation with four corners of the grid cell to determine value function at
    
    //Plan: Bilinear on upper face and lower face of interp cube, then linear interpolation between the two faces to get value at the point in the cube.
    
    //Mallory Note: For now, I'd prefer to NOT do interpolation in time as well. With that said, we need to think about the grid / timescales to make sure they make sense. Also, be careful with typecasting!
    
    //Current number of mosquitoes of each type at the current  timestep.
    int currentNumR = gNumMosquitoesPerNode * aGridNumR;
    int currentNumFW = gNumMosquitoesPerNode * aGridNumFW;
    int currentNumMW = gNumMosquitoesPerNode * aGridNumMW;
    int currentNumF = currentNumR * gAlpha + currentNumFW;
    int currentNumM = currentNumR * (1 - gAlpha) + currentNumMW;
    
    //Initialize placeholders for optimal controls and value function values
    int optimalAWF;
    int optimalAWM;
    double optimalValueFunction;

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
            int nextR = currentNumR + gDT * (gRate * (1 - ((currentNumF + currentNumM) / gK)) * ((gAlpha * currentNumR * (1 - gAlpha) * currentNumR) / (gAlleleEffect + currentNumM)) - gDelta * currentNumR);
            int nextFW = currentNumFW + gDT * (gRate * gAlpha *  (1 - ((currentNumF + currentNumM) / gK)) * (currentNumFW * (currentNumM / (gAlleleEffect + currentNumM))) - gDelta * currentNumFW + currentWFControl);
            int nextMW = currentNumMW + gDT * (gRate * (1 - gAlpha) * (1 - ((currentNumF + currentNumM) / gK)) * (currentNumFW * (currentNumM / (gAlleleEffect + currentNumM))) - gDelta * currentNumMW + currentWMControl);
            
            //check to see if any of the next values are out of bounds
            if ((nextR > gRegularMaxNodes * gNumMosquitoesPerNode) || (nextFW > gFemaleWMaxNodes * gNumMosquitoesPerNode) || (nextMW > gMaleWMaxNodes * gNumMosquitoesPerNode)){
                currentValueFunction = gInfinity;
            }//end check on bounds
            
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
                
                currentValueFunction = interpValue + runningCost(currentWFControl, currentWMControl, currentNumR, gDT); //Note: assuming tau = gDT, hence why no interp in time now
                
            }//end interpolation
            
            //check to see if this value is less than the current saved optimal.
            if (currentValueFunction < optimalValueFunction){
                optimalValueFunction = currentValueFunction;
                optimalAWF = currentWFControl;
                optimalAWM = currentWMControl;
            }
            
        }//end inner control loop
    }//end outer control loop
    
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
        for(short int k = 0; k < gRegularMaxNodes + 1; k++) {
            for(short int i = 0; i < gFemaleWMaxNodes + 1; i++) {// fw loop
                for(short int j = 0; j < gMaleWMaxNodes + 1; j++) {// mw loop
                    optimalValue(aValueFunction, aFemaleWControls, aMaleWControls, k, i, j, t); //compute u(d,v,t)
              }//end mw loop
            } //end fw loop
        } //end r loop
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
    for(int t = gNt - 1; t >= 0; t--) {//time loop
        for(short int k = 0; k < gRegularMaxNodes + 1; k++) {
            for(short int i = 0; i < gFemaleWMaxNodes + 1; i++) {// fw loop
                for(short int j = 0; j < gMaleWMaxNodes + 1; j++) {// mw loop
                    
                    (*aArray)[i][j][k][t] = 0; //check on this
                    
              }//end mw loop
            } //end fw loop
        } //end r loop
    }//end time loop
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
