// Author: Akash Patel (apatel435@gatech.edu)

// Methods that deal with file operations

// Includes
#include <dart/dart.hpp>
#include <iostream>
#include <fstream>

// Namespaces
using namespace std;
using namespace dart::math;

// Defines
#define MAXBUFSIZE ((int) 1e6)

// Function Prototypes
// // Read file as matrix
Eigen::MatrixXd readInputFileAsMatrix(string inputPosesFilename);

// // Extract filename
string extractFilename(string filename);
