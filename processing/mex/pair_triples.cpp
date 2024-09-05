#include "mex.h"
#include <vector>
#include <array>
#include <iostream>
#include <cmath>

#define N_ROWS 9
#define ID 0
#define X 1
#define Y 2
#define Z 3
#define TIME 4
#define EN 5
#define TX 6
#define AX 7
#define CP 8

#if USE_SECREJ
    #define BASIC_ARGS 8
#else //other policy no need for secondary window
    #define BASIC_ARGS 7
#endif

struct event {
    double eID, x, y, z, time, energy, tx, ax, cp;
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


    // Check for the right number of input and output arguments
    bool isGeoRej;
    if (nrhs == BASIC_ARGS) {
        isGeoRej = false;
    } else if (nrhs == BASIC_ARGS + 8) {
        isGeoRej = true;
    } else {
        mexErrMsgTxt("Invalid number of input arguments.");
    }

    if (nlhs != 1) {
        mexErrMsgTxt("Invalid number of output arguments.");
    }

    // Check input
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 2 || mxGetM(prhs[0]) != N_ROWS) {
        mexErrMsgTxt("Input must be a real matrix.");
    }  

    if (!mxIsDouble(prhs[1]) || !mxIsScalar(prhs[1])) {
        mexErrMsgTxt("coincWidth must be a double scalar.");
    } 

    if (!mxIsDouble(prhs[2]) || !mxIsScalar(prhs[2])) {
        mexErrMsgTxt("coincStartTime must be a double scalar.");
    } 
        
    if (!mxIsDouble(prhs[3]) || mxGetM(prhs[3]) != 1 || mxGetN(prhs[3]) != 2) {
        mexErrMsgTxt("PGWin must be a 1x2 double.");
    }

    if (!mxIsDouble(prhs[4]) || !mxIsScalar(prhs[4])) {
        mexErrMsgTxt("pGammaMiddleTime must be a double scalar.");
    }
        
    if (!mxIsDouble(prhs[5]) || mxGetM(prhs[5]) != 1 || mxGetN(prhs[5]) != 2) {
        mexErrMsgTxt("annihEn must be a 1x2 array.");
    }
        
    if (!mxIsDouble(prhs[6]) || mxGetM(prhs[6]) != 1 || mxGetN(prhs[6]) != 2) {
        mexErrMsgTxt("pGammaEn must be a 1x2 array.");
    }

    // parameters initialization
    double coincWidth = mxGetDoubles(prhs[1])[0];
    double coincStartTime = mxGetDoubles(prhs[2])[0];
    double pGammaWidth = mxGetDoubles(prhs[3])[1] - mxGetDoubles(prhs[3])[0];
    double pGammaMiddleTime = mxGetDoubles(prhs[4])[0]; // offset of PG window
    double pGammaStartTime = mxGetDoubles(prhs[3])[0] + pGammaMiddleTime;
    double LowAn = mxGetDoubles(prhs[5])[0];
    double HighAn = mxGetDoubles(prhs[5])[1];
    double LowPG = mxGetDoubles(prhs[6])[0];
    double HighPG = mxGetDoubles(prhs[6])[1];
    // to be filled later
    int minTxDiff_p_a = 0;
    int maxTxDiff_p_a = 0;
    int minAxDiff_p_a = 0;
    int minTxDiff_a_a = 0;
    int maxTxDiff_a_a = 0;
        
#if USE_SECREJ
    if (!mxIsDouble(prhs[7]) || mxGetM(prhs[7]) != 1 || mxGetN(prhs[7]) != 2) {
        mexErrMsgTxt("secWin must be a 1x2 double.");
    }
    double secWinStart = mxGetDoubles(prhs[7])[0];
    double secWinEnd = mxGetDoubles(prhs[7])[1];
#endif
        
    if (isGeoRej) {
        if (!mxIsScalar(prhs[BASIC_ARGS]) || !mxIsDouble(prhs[BASIC_ARGS])) {
            mexErrMsgTxt("minTxDiff_p_a must be a double scalar.");
        } else {
            minTxDiff_p_a = (int) mxGetDoubles(prhs[BASIC_ARGS])[0];
        }
            
        if (!mxIsScalar(prhs[BASIC_ARGS+1]) || !mxIsDouble(prhs[BASIC_ARGS+1])) {
            mexErrMsgTxt("maxTxDiff_p_a must be a double scalar.");
        } else {
            maxTxDiff_p_a = (int) mxGetDoubles(prhs[BASIC_ARGS+1])[0];
        }

        if (!mxIsScalar(prhs[BASIC_ARGS+2]) || !mxIsDouble(prhs[BASIC_ARGS+2])) {
            mexErrMsgTxt("minAxDiff_p_a must be a double scalar.");
        } else {
            minAxDiff_p_a = (int) mxGetDoubles(prhs[BASIC_ARGS+2])[0];
        }
            
        if (!mxIsScalar(prhs[BASIC_ARGS+3]) || !mxIsDouble(prhs[BASIC_ARGS+3])) {
            mexErrMsgTxt("minTxDiff_a_a must be a double scalar.");
        } else {
            minTxDiff_a_a = (int) mxGetDoubles(prhs[BASIC_ARGS+3])[0];
        }
            
        if (!mxIsScalar(prhs[BASIC_ARGS+4]) || !mxIsDouble(prhs[BASIC_ARGS+4])) {
            mexErrMsgTxt("maxTxDiff_a_a must be a double scalar.");
        } else {
            maxTxDiff_a_a = (int) mxGetDoubles(prhs[BASIC_ARGS+4])[0];
        }
            
        if (!mxIsDouble(prhs[BASIC_ARGS+5]) || mxIsComplex(prhs[BASIC_ARGS+5]) || mxGetNumberOfDimensions(prhs[BASIC_ARGS+5]) != 2 || mxGetM(prhs[BASIC_ARGS+5]) != 2) {
            mexErrMsgTxt("m_xtal_pos_xy must be a 2xn real array.");
        }
            
        if (!mxIsDouble(prhs[BASIC_ARGS+6]) || mxIsComplex(prhs[BASIC_ARGS+6]) || mxGetNumberOfDimensions(prhs[BASIC_ARGS+6]) != 2 || mxGetM(prhs[BASIC_ARGS+6]) != 1)
            mexErrMsgTxt("m_ring_z must be a row vector.");

        if (!mxIsDouble(prhs[BASIC_ARGS+7]) || mxIsComplex(prhs[BASIC_ARGS+7]) || mxGetNumberOfDimensions(prhs[BASIC_ARGS+7]) != 2 || mxGetM(prhs[BASIC_ARGS+7]) != 2 || mxGetN(prhs[BASIC_ARGS+7]) != 3)
            mexErrMsgTxt("fov must be a 2x3 real matrix.");
    }
        
#if DEBUG        
    mexPrintf("coincStartTime: %f\n", coincStartTime*1000000000.);
    mexPrintf("coincWidth: %f\n", coincWidth*1000000000.);
    mexPrintf("pGammaStartTime: %f\n", pGammaStartTime*1000000000.);
    mexPrintf("pGammaWidth: %f\n", pGammaWidth*1000000000.);
    mexPrintf("pGammaMiddleTime: %f\n", pGammaMiddleTime*1000000000.);
#if USE_SECREJ
    mexPrintf("secWinStart: %f\n", secWinStart*1000000000.);
    mexPrintf("secWinEnd: %f\n", secWinEnd*1000000000.);
#endif 
    mexPrintf("isGeoRej: %s\n", isGeoRej ? "true" : "false");
    mexPrintf("LowAn: %f\n", LowAn);
    mexPrintf("HignAn: %f\n", HighAn);
    mexPrintf("LowPG: %f\n", LowPG);
    mexPrintf("HighPG: %f\n", HighPG);
    mexPrintf("minTxDiff_p_a: %d\n", minTxDiff_p_a);
    mexPrintf("maxTxDiff_p_a: %d\n", maxTxDiff_p_a);
    mexPrintf("minAxDiff_p_a: %d\n", minAxDiff_p_a);
    mexPrintf("minTxDiff_a_a: %d\n", minTxDiff_a_a);
    mexPrintf("maxTxDiff_a_a: %d\n", maxTxDiff_a_a);
#endif


    // Get the input matrix dimensions
    mwSize n = mxGetN(prhs[0]);

    // Indexes
    unsigned long long coincEndIndex = 0;
    unsigned long long coincStartIndex = 0;
    unsigned long long pGammaEndIndex = 0;
    unsigned long long pGammaStartIndex = 0;

    // Temporary variables
    double dt; // temp lifetime measurement
    std::array<double, N_ROWS> columnTemp1;
    std::array<double, N_ROWS> columnTemp2;
    std::array<double, N_ROWS> columnTemp3;
    // Use a vector to store the columns with the first element as zero
    std::vector<std::array<double, N_ROWS>> outputColumns;
    

    // Iterate through columns
    for (unsigned long long refIndex = 0; refIndex < n; refIndex++)
    {
        if (mxGetDoubles(prhs[0])[refIndex*N_ROWS + EN] > HighAn || mxGetDoubles(prhs[0])[refIndex*N_ROWS + EN] < LowAn) {
            continue;
        }

        // Prepare start and end index for qualified second 511 and the PG. 
        while (mxGetDoubles(prhs[0])[pGammaEndIndex * N_ROWS + TIME] - mxGetDoubles(prhs[0])[refIndex * N_ROWS + TIME] 
            <= pGammaStartTime + pGammaWidth + coincWidth) {
                pGammaEndIndex++;
                if (pGammaEndIndex == n) {
                    break;
                }
        } 

        if (pGammaEndIndex == n) {
            break;
        }
        
        while (mxGetDoubles(prhs[0])[pGammaStartIndex * N_ROWS + TIME] - mxGetDoubles(prhs[0])[refIndex * N_ROWS + TIME]
                < pGammaStartTime) {
            pGammaStartIndex++;
        }
            
        while (mxGetDoubles(prhs[0])[coincEndIndex * N_ROWS + TIME] - mxGetDoubles(prhs[0])[refIndex * N_ROWS + TIME]
                <= coincStartTime + coincWidth) {
            coincEndIndex++;
            if (coincEndIndex == n) {
                break;
            }
        }

        if (coincEndIndex == n) {
            break;
        }

        while (mxGetDoubles(prhs[0])[coincStartIndex * N_ROWS + TIME] - mxGetDoubles(prhs[0])[refIndex * N_ROWS + TIME]
                < coincStartTime || refIndex == coincStartIndex) {
            if (coincStartIndex == coincEndIndex) {
                break; //avoid surpassing when coinc is delayed and coincEndIndex==refIndex
            }
            coincStartIndex++;
        }

        if (coincEndIndex <= coincStartIndex || pGammaEndIndex <= pGammaStartIndex) {
            continue;
        }

        // Search eligible 511 and PG singles
        for (unsigned long long ia = coincStartIndex; ia < coincEndIndex; ia++) {

            if (mxGetDoubles(prhs[0])[ia*N_ROWS + EN] > HighAn || mxGetDoubles(prhs[0])[ia*N_ROWS + EN] < LowAn) {
                continue;
            }
#if USE_TAKEFIRST
            bool isFirstFromLeft = true;
#endif
            for (unsigned long long ip = pGammaStartIndex; ip < pGammaEndIndex; ip++) {

                if (mxGetDoubles(prhs[0])[ip*N_ROWS + EN] < LowPG || mxGetDoubles(prhs[0])[ip*N_ROWS + EN] > HighPG || ip == refIndex) {
                    continue;
                }

                if (isGeoRej) {
                    int tx1 = (int) mxGetDoubles(prhs[0])[refIndex*N_ROWS + TX];
                    int ax1 = (int) mxGetDoubles(prhs[0])[refIndex*N_ROWS + AX];
                    int tx2 = (int) mxGetDoubles(prhs[0])[ia*N_ROWS + TX];
                    int ax2 = (int) mxGetDoubles(prhs[0])[ia*N_ROWS + AX];
                    int tx3 = (int) mxGetDoubles(prhs[0])[ip*N_ROWS + TX];
                    int ax3 = (int) mxGetDoubles(prhs[0])[ip*N_ROWS + AX];
                    double x1 = mxGetDoubles(prhs[BASIC_ARGS+5])[2 * tx1];
                    double y1 = mxGetDoubles(prhs[BASIC_ARGS+5])[2 * tx1 + 1];
                    double z1 = mxGetDoubles(prhs[BASIC_ARGS+6])[ax1];
                    double x2 = mxGetDoubles(prhs[BASIC_ARGS+5])[2 * tx2];
                    double y2 = mxGetDoubles(prhs[BASIC_ARGS+5])[2 * tx2 + 1];
                    double z2 = mxGetDoubles(prhs[BASIC_ARGS+6])[ax2];

                    if ((std::abs(tx1 - tx3) < minTxDiff_p_a || std::abs(tx1 - tx3) > maxTxDiff_p_a) && std::abs(ax1 - ax3) < minAxDiff_p_a) {
                        continue;
                    } else if ((std::abs(tx2 - tx3) < minTxDiff_p_a || std::abs(tx2 - tx3) > maxTxDiff_p_a) && std::abs(ax2 - ax3) < minAxDiff_p_a){
                        continue;
                    } else if (std::abs(tx1 - tx2) < minTxDiff_a_a || std::abs(tx1-tx2) > maxTxDiff_a_a) {
                        continue;
                    } else { 

                        double fov_x1 = mxGetDoubles(prhs[BASIC_ARGS+7])[0];
                        double fov_x2 = mxGetDoubles(prhs[BASIC_ARGS+7])[1];
                        double fov_y1 = mxGetDoubles(prhs[BASIC_ARGS+7])[2];
                        double fov_y2 = mxGetDoubles(prhs[BASIC_ARGS+7])[3];
                        double fov_z1 = mxGetDoubles(prhs[BASIC_ARGS+7])[4];
                        double fov_z2 = mxGetDoubles(prhs[BASIC_ARGS+7])[5];

                        if ((((y2-y1)*(fov_z1-z1)/(z2-z1+1e-7)+y1) < fov_y1 || ((y2-y1)*(fov_z1-z1)/(z2-z1+1e-7)+y1) > fov_y2 || ((x2-x1)*(fov_z1-z1)/(z2-z1+1e-7)+x1) < fov_x1 || ((x2-x1)*(fov_z1-z1)/(z2-z1+1e-7)+x1) > fov_x2) && 
                            (((y2-y1)*(fov_z2-z1)/(z2-z1+1e-7)+y1) < fov_y1 || ((y2-y1)*(fov_z2-z1)/(z2-z1+1e-7)+y1) > fov_y2 || ((x2-x1)*(fov_z2-z1)/(z2-z1+1e-7)+x1) < fov_x1 || ((x2-x1)*(fov_z2-z1)/(z2-z1+1e-7)+x1) > fov_x2) &&
                            (((z2-z1)*(fov_y1-y1)/(y2-y1+1e-7)+z1) < fov_z1 || ((z2-z1)*(fov_y1-y1)/(y2-y1+1e-7)+z1) > fov_z2 || ((x2-x1)*(fov_y1-y1)/(y2-y1+1e-7)+x1) < fov_x1 || ((x2-x1)*(fov_y1-y1)/(y2-y1+1e-7)+x1) > fov_x2) &&
                            (((z2-z1)*(fov_y2-y1)/(y2-y1+1e-7)+z1) < fov_z1 || ((z2-z1)*(fov_y2-y1)/(y2-y1+1e-7)+z1) > fov_z2 || ((x2-x1)*(fov_y2-y1)/(y2-y1+1e-7)+x1) < fov_x1 || ((x2-x1)*(fov_y2-y1)/(y2-y1+1e-7)+x1) > fov_x2) &&
                            (((z2-z1)*(fov_x1-x1)/(x2-x1+1e-7)+z1) < fov_z1 || ((z2-z1)*(fov_x1-x1)/(x2-x1+1e-7)+z1) > fov_z2 || ((y2-y1)*(fov_x1-x1)/(x2-x1+1e-7)+y1) < fov_y1 || ((y2-y1)*(fov_x1-x1)/(x2-x1+1e-7)+y1) > fov_y2) &&
                            (((z2-z1)*(fov_x2-x1)/(x2-x1+1e-7)+z1) < fov_z1 || ((z2-z1)*(fov_x2-x1)/(x2-x1+1e-7)+z1) > fov_z2 || ((y2-y1)*(fov_x2-x1)/(x2-x1+1e-7)+y1) < fov_y1 || ((y2-y1)*(fov_x2-x1)/(x2-x1+1e-7)+y1) > fov_y2)) {
                            continue;
                        }

                    }
                }

#if USE_SECREJ
                // if additional 511 in the window is >= 1, reject this triple
                unsigned long long ip1 = ip;
                unsigned long long ip2 = ip;
                // rewind to the lowest time tag inside the secondary window
                double tSecPast = mxGetDoubles(prhs[0])[ip*N_ROWS + TIME] + secWinStart;
                double tSecFuture = mxGetDoubles(prhs[0])[ip*N_ROWS + TIME] + secWinEnd;
                while (mxGetDoubles(prhs[0])[ip1*N_ROWS + TIME] >= tSecPast && ip1 > 0) {
                    ip1--;
                }
                if (mxGetDoubles(prhs[0])[ip1*N_ROWS + TIME] < tSecPast) {
                    ip1++; //
                }
                // move to the first time tag outside to the right side of the window
                while (mxGetDoubles(prhs[0])[ip2*N_ROWS + TIME] <= tSecFuture && ip2 < n-1) {
                    ip2++;
                }
                if (mxGetDoubles(prhs[0])[ip2*N_ROWS + TIME] <= tSecFuture) {
                    ip2++; // make sure it counts the last one
                }

                // count 511s in the window, should not include the refIndex, ia, and possibly ip
                unsigned long countSec = 0;
                for (unsigned long long ip3 = ip1; ip3 < ip2; ip3++) {
                    if ((mxGetDoubles(prhs[0])[ip3*N_ROWS + EN] < HighAn && mxGetDoubles(prhs[0])[ip3*N_ROWS + EN] > LowAn) && (ip3 != refIndex) && (ip3 != ia) && (ip3 != ip)) {
                        countSec++;
                    }
                }
                // decision
                if (countSec >= 1) {
                    continue;
                }
#endif
                
                dt = 0.5 * (mxGetDoubles(prhs[0])[refIndex*N_ROWS + TIME] + mxGetDoubles(prhs[0])[ia*N_ROWS + TIME] - coincStartTime) - (mxGetDoubles(prhs[0])[ip*N_ROWS + TIME] - pGammaMiddleTime);
                if (dt > (pGammaMiddleTime - pGammaStartTime) || -dt > (pGammaStartTime + pGammaWidth - pGammaMiddleTime)) {
                    continue;
                }

                // Copy the columnTemp to the temp vector
                for (unsigned int row = 0; row < N_ROWS; row++)
                {
                    columnTemp1[row] = mxGetDoubles(prhs[0])[refIndex*N_ROWS + row];
                    columnTemp2[row] = mxGetDoubles(prhs[0])[ia*N_ROWS + row];
                    columnTemp3[row] = mxGetDoubles(prhs[0])[ip*N_ROWS + row];
                }

#if USE_TAKEFIRST // if a triple has been identified in the ip loop previously and dt > 0, replace the previous one
                if (isFirstFromLeft) {
                    isFirstFromLeft = false;
                } else if (!isFirstFromLeft && dt > 0) {
                    outputColumns.pop_back();
                    outputColumns.pop_back();
                    outputColumns.pop_back();
                }    
#endif 

                outputColumns.push_back(columnTemp1);
                outputColumns.push_back(columnTemp2);
                outputColumns.push_back(columnTemp3);  
            }
        }
        
     
    }
    

    // Create an output matrix with dynamically allocated size
    mwSize outputCols = outputColumns.size();
    plhs[0] = mxCreateDoubleMatrix(N_ROWS, outputCols, mxREAL);
    double *outputData = mxGetDoubles(plhs[0]);

    // Copy columns with the first element as zero to the output matrix
    for (mwIndex colIndex = 0; colIndex < outputCols; colIndex++)
    {
        for (mwIndex row = 0; row < N_ROWS; row++)
        {
            outputData[colIndex * N_ROWS + row] = outputColumns[colIndex][row];
        }
    }
}
