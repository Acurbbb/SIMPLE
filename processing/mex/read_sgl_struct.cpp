//mex -R2018a read_sgl_struct.cpp
#include "mex.h"
#include <cstdio>

#define N_ROWS 9

struct MyData {
	// eventID to be casted from int (G4int)
	// edep is the energy deposited in each hit (in MeV)
	uint32_T eventID;
	float sourcePosX, sourcePosY, sourcePosZ;
	double time;
	float energy;
	int16_T txID, axID;
	uint8_T comptonPhantom;
} __attribute__((packed)); // avoid struct padding https://www.geeksforgeeks.org/how-to-avoid-structure-padding-in-c/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check input arguments
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "One input required.");
    }

    // Check if the input is a string
    if (!mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notString", "Input must be a filename string.");
    }

    // Get the filename from the input
    char *filename;
    filename = mxArrayToString(prhs[0]);

    // Open the binary file for reading in binary mode
    FILE* inputFile = fopen(filename, "rb");

    // Check if the file is opened successfully
    if (!inputFile) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:fileOpenError", "Error opening file: %s", filename);
    }

    // Determine the dataSize based on the file size
    fseek(inputFile, 0, SEEK_END);
    long fileSize = ftell(inputFile);
    rewind(inputFile);
    const int dataSize = fileSize / sizeof(MyData);

    // Declare an array of doubles to store your data
    MyData* dataArray = new MyData[dataSize]; // Each entry has three elements

    // Read binary data into the array
    fread(dataArray, sizeof(MyData), dataSize, inputFile);

    // Close the file
    fclose(inputFile);

    // Create MATLAB output array
    plhs[0] = mxCreateDoubleMatrix(N_ROWS, dataSize, mxREAL);
    double* outputArray = mxGetPr(plhs[0]);

    // Copy data to the output array
    for (int i = 0; i < dataSize ; ++i) {
        outputArray[i*N_ROWS] = (double) dataArray[i].eventID;
        outputArray[i*N_ROWS + 1] = (double) dataArray[i].sourcePosX;
        outputArray[i*N_ROWS + 2] = (double) dataArray[i].sourcePosY;
        outputArray[i*N_ROWS + 3] = (double) dataArray[i].sourcePosZ;
        outputArray[i*N_ROWS + 4] = (double) dataArray[i].time;
        outputArray[i*N_ROWS + 5] = (double) dataArray[i].energy;
        outputArray[i*N_ROWS + 6] = (double) dataArray[i].txID;
        outputArray[i*N_ROWS + 7] = (double) dataArray[i].axID;
        outputArray[i*N_ROWS + 8] = (double) dataArray[i].comptonPhantom;
    }

    // Clean up allocated memory
    delete[] dataArray;
}
