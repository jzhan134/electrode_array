#pragma once
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <thread>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cmath>
#include <memory>
#include <random>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <queue>
#include "omp.h"
#include <cassert>
#include <map>

using namespace std;


struct BD_Particle{
    double x;
    double y;
    double Ex;
    double Ey;
    double D;
    double Fx;
    double Fy;
    
};

// const int np = 1000;      // particle number
// size of field look-up table
