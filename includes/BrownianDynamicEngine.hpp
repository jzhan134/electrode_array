#ifndef BROWNIANDYNAMICENGINE_HPP
#define BROWNIANDYNAMICENGINE_HPP
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
#include <set>

using namespace std;
extern const int total_time;
extern ofstream summary;

struct BD_Particle{
    std::pair<double, double> coord;
    std::pair<double, double> F;
    std::pair<double, double> EP;
    std::pair<double, double> DEP;
    double D;
    int group;
    BD_Particle( std::pair<double, double> coord_) : coord(coord_){}
};
class BrownianDynamicEngine{
public:
    BrownianDynamicEngine(){
        prevCheckSum = 0;
        ReadDiffusivity();
    }
    ~BrownianDynamicEngine(){
        if (trajOs.is_open()){ trajOs.close();}
    }
    
    void initialization
    (
        std::string startingConfig
    );
    
    void DEP_quench
    (
        string fieldFileName,
        double fieldStrength_, 
        double execTime
    );

    void EP_move(
        int targetNum
    );
        
    void GenerateRandomConfig();

    void ReadSavedConfig
    (
        std::string
    );
    
    void CalDss(); 

    void barycentricInterpolation(
        bool isEP
    );
        
    void ReadDiffusivity();
    
    void ReadE
    (
        std::string filename, 
        bool isEP
    );
        
    void OutputTrajectory
    (
        ostream& osFileName
    );

    inline int countNum
    (
        int gridX, 
        int gridY
    ){
        int num = 0;
        for (int i = 0; i < p.size(); i++){
            if (fabs(p[i].coord.first - (gridX-1)*60000) <= 30000 && 
                fabs(p[i].coord.second- (gridY-1)*60000) <= 30000){
                num++;
            }
        }
        return num;
    }

    inline double Distance2D
    (
        std::vector<BD_Particle>::iterator p1, 
        std::vector<BD_Particle>::iterator p2
    ) const {
            double x1 = p1->coord.first;
            double x2 = p2->coord.first;
            double y1 = p1->coord.second;
            double y2 = p2->coord.second; 
        return sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
    }

protected:
    // constants
    const double pi = 3.14159265;
    const int nstep = 10000;        // steps to simulate 1s
    const double dt = 0.1;          // (ms) simulation time of each step
    const double a = 800;        // particle radius (nm)
    const double tempr = 20.0;      // temperature (C)
    const double fcm = -0.4667;     //Clausius–Mossotti Factor
    const double kb = 1.380658e-23; // Boltzmann constant
    const double kappa = a/10;      // Debye length
    const double pfpp = 2.2975*a;   // electrostatic repulsion
    const double fac1 = 5.9582e7/a;
    const double fac2 = 40.5622*sqrt((273 + tempr) / a) / sqrt(dt);
    const double Fhw = 0.417; // hard-sphere repulsion 
    const double rcut = 5.0*a; // range of neighbor for force calculation
    const double ppPreFactor = 1e9 * (273 + tempr) * kb * 21.86 * 6.25 * 8 * pow(a,3);
    const double pfPreFactor = 1e9 * 2 * pi * 80 * 8.85e-12 * pow(a*1e-9, 3) * fcm;
    const double epPreFactor = 1.5 * 80 * 8.85e-12 * 1e-4 / (a*1e-9);
    
    
    static const int rgdssbin = 25;
    static const int distdssbin = 50;
    int dsscount[rgdssbin][distdssbin];
    double dssarray[rgdssbin][distdssbin]; 
    
    
    // path of particle starting configurations
    std::string configPath = "./library/configurations/";

    //BD_Particle p[np]; // coordinate and force of each particle
    std::vector<BD_Particle> p;
    const int defaultParticleNumber = 400;
    const double periodicWindow = 30000;
    
    // data structure for lookup tables
    static const int E_tot = 801; 
    const double d_idx = 250; // (nm) resolution of field loop-up table
    double EP_Table[E_tot][E_tot][2];
    double DEP_Table[E_tot][E_tot][2];
    
    double fieldStrength; // applied field strength with respect to the table
    int elapsedTime; // elapsed time in a complete simulation cycle
    static int repCycle;
    int prevCheckSum;
    ofstream trajOs; // output file names

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> rand_normal{0.0,1.0};
    std::uniform_real_distribution<> rand_uniform{0.0,1.0};
};

#endif /* BROWNIANDYNAMICS_HPP */