#include "BrownianDynamicEngine.hpp"
int BrownianDynamicEngine::repCycle(0);
/***************************** DYNAMICS SIMULATION ****************************/
void BrownianDynamicEngine::initialization(
        std::string startingConfig)
{
    // reset timer for the current cycle
    elapsedTime = 0;
    
    // update finished number of cycles
    stringstream repCycleOS;
    std::string cnt = std::to_string(repCycle++);
    
    // read initial configuration
    p.clear();
    if (startingConfig.compare("random") == 0) {
        GenerateRandomConfig();
    } else {
        ReadSavedConfig(configPath + startingConfig);
    }
    
    // create new trajectory output file
    if (trajOs.is_open()) trajOs.close();
    trajOs.open("./traj/xyz" + cnt + ".dat");  
    // OutputTrajectory(trajOs);
}

// time step: 1e-4s
// total time: 0.1s
// output freq: 0.1s
void BrownianDynamicEngine::EP_move(
    int targetNum)  
{
    double rijsep,xij, yij;
    double dx, dy;
    double randx(0), randy(0);
    double EPStrength;
    double pfFactor = pfPreFactor * pow(0,2);
    int iniNum = countNum(1, 1);

    ReadE("./library/fields/array/3x3_EP.txt", true);
    ReadE("./library/fields/array/test1.txt", false);
    
    while (countNum(1, 1) > targetNum){
        
        // update local particle diffusivity and field magnitudes per second of simulation 
        CalDss();
        barycentricInterpolation(true);
        
        // simulation is discretized into 10000 steps of 0.1ms
        for (int step = 0; step < nstep; step++){ 
            EPStrength = (countNum(1, 1) == iniNum) ? 0.4 : 0.05;
            std::for_each(p.begin(), p.end(), [](BD_Particle& p) {p.F = std::make_pair(0.0,0.0);});

            #pragma omp parallel for
            for (auto p1 = p.begin(); p1 != p.end(); ++p1) {
                for (auto p2 = std::next(p1); p2 != p.end(); ++p2){
                    rijsep = Distance2D(p1,p2);
                    if (rijsep < rcut){
                        xij = p2->coord.first - p1->coord.first;
                        yij = p2->coord.second - p1->coord.second;
                        if (rijsep < 2*a + 40){
                            dx = xij * (2*a + 40 - rijsep)  /sqrt(xij * xij + yij * yij);
                            dy = yij * (2*a + 40 - rijsep)  /sqrt(xij * xij + yij * yij);
                            p1->coord.first -= dx;
                            p1->coord.second -= dy;
                            p2->coord.first += dx;
                            p2->coord.second += dy;
                        } else {
                            double Fpp = 1e18*kb*(tempr+273)*kappa*pfpp*exp(-kappa*(rijsep-2.0*a)/a)/a;
                            p1->F.first -= Fpp*xij/rijsep;
                            p1->F.second -= Fpp*yij/rijsep;
                            p2->F.first += Fpp*xij/rijsep;
                            p2->F.second += Fpp*yij/rijsep;
                        }
                    }
                }
                // EP particle-field interaction
                p1->F.first += epPreFactor*EPStrength*p1->EP.first;
                p1->F.second += epPreFactor*EPStrength*p1->EP.second;

                // DEP particle-field interaction
                p1->F.first += pfFactor*p1->DEP.first;
                p1->F.second += pfFactor*p1->DEP.second;

                // random displacement
                randx =  rand_normal(gen) * sqrt(1.0 / p1->D);
                randy =  rand_normal(gen) * sqrt(1.0 / p1->D);

                // final displacement update is the summation of deterministic and random
                p1->coord.first += p1->D * (p1->F.first*fac1 + randx*fac2) * dt;
                p1->coord.second += p1->D * (p1->F.second*fac1 + randy*fac2) * dt;

                // periodic boundary
                if (p1->coord.first > periodicWindow) p1->coord.first = (p1->coord.first - 2*periodicWindow);
                if (p1->coord.first  < -periodicWindow) p1->coord.first  = (p1->coord.first + 2*periodicWindow);
                if (p1->coord.second > periodicWindow) p1->coord.second = (p1->coord.second - 2*periodicWindow);
                if (p1->coord.second < -periodicWindow) p1->coord.second = (p1->coord.second + 2*periodicWindow);
            }
            if ((step+1) % 1000 == 0){
                std::cout << "move " << (double) (elapsedTime++)/10 << endl;
                OutputTrajectory(trajOs);
            }

            // exit as soon as target number is obtained
            if (step % 500 == 0 && countNum(1, 1) <= targetNum) return;
        }
    }
}

// time step: 1e-4s
// total time: 0.1s
// output freq: 0.1s
void BrownianDynamicEngine::DEP_quench(
    string fieldFileName,
    double fieldStrength_, 
    double execTime)  
{
    double rijsep,xij, yij;
    double dx, dy;
    double randx(0), randy(0);
    double pfFactor = pfPreFactor * pow(fieldStrength_,2);
    int iniNum = countNum(1, 1);
    
    ReadE(fieldFileName, false);
    
    for (int t = 0; t < execTime*10; t++){
        // update local particle diffusivity and field magnitudes per second of simulation 
        CalDss();
        barycentricInterpolation(false);

        // update particle group, the group is used to reduce amount of pair-wise 
        // comparison if two particles belong to different groups
        for (auto it = p.begin(); it != p.end(); ++it) {
            it->group = round(it->coord.first/60000)*3 + round(it->coord.second/60000);
        }
        
        // simulation is discretized into 10000 steps of 0.1ms
        for (int step = 0; step < nstep/10; step++){ 
            std::for_each(p.begin(), p.end(), [](BD_Particle& p) {p.F = std::make_pair(0.0,0.0);});
            #pragma omp parallel for
            for (auto p1 = p.begin(); p1 != p.end(); ++p1) {
                for (auto p2 = std::next(p1); p2 != p.end(); ++p2){
                    if (p1->group == p2->group){
                        rijsep = Distance2D(p1,p2);
                        if (rijsep < rcut){
                            xij = p2->coord.first - p1->coord.first;
                            yij = p2->coord.second - p1->coord.second;
                            if (rijsep < 2*a + 40){
                                dx = xij * (2*a + 40 - rijsep)  /sqrt(xij * xij + yij * yij);
                                dy = yij * (2*a + 40 - rijsep)  /sqrt(xij * xij + yij * yij);
                                p1->coord.first -= dx;
                                p1->coord.second -= dy;
                                p2->coord.first += dx;
                                p2->coord.second += dy;
                            } else {
                                double Fpp = 1e18*kb*(tempr+273)*kappa*pfpp*exp(-kappa*(rijsep-2.0*a)/a)/a;
                                p1->F.first -= Fpp*xij/rijsep;
                                p1->F.second -= Fpp*yij/rijsep;
                                p2->F.first += Fpp*xij/rijsep;
                                p2->F.second += Fpp*yij/rijsep;
                            }
                        }
                    }
                }
                // DEP particle-field interaction
                p1->F.first += pfFactor*p1->DEP.first;
                p1->F.second += pfFactor*p1->DEP.second;

                // random displacement
                randx =  rand_normal(gen) * sqrt(1.0 / p1->D);
                randy =  rand_normal(gen) * sqrt(1.0 / p1->D);

                // final displacement update is the summation of deterministic and random
                p1->coord.first += p1->D * (p1->F.first*fac1 + randx*fac2) * dt;
                p1->coord.second += p1->D * (p1->F.second*fac1 + randy*fac2) * dt;

                // periodic boundary
                if (p1->coord.first > periodicWindow) p1->coord.first = (p1->coord.first - 2*periodicWindow);
                if (p1->coord.first  < -periodicWindow) p1->coord.first  = (p1->coord.first + 2*periodicWindow);
                if (p1->coord.second > periodicWindow) p1->coord.second = (p1->coord.second - 2*periodicWindow);
                if (p1->coord.second < -periodicWindow) p1->coord.second = (p1->coord.second + 2*periodicWindow);
            }
            if ((step+1) % 1000 == 0){
                std::cout << (double) (elapsedTime++)/10 << endl;
                OutputTrajectory(trajOs);
            }
        }
    }
}

// at the beginning of every simulation second, the field applied on
// every particle is approximated based on its instantaneous position
// with barycentric interpolation
void BrownianDynamicEngine::barycentricInterpolation(bool isEP)
{
    int x[3],y[3];
    double dx,dy;
    int idx_x, idx_y;
    double det, lambda1, lambda2, lambda3;

    for (auto p1 = p.begin(); p1 != p.end(); ++p1){
        dx = p1->coord.first/d_idx + (E_tot-1)/2;
        dy = p1->coord.second/d_idx + (E_tot-1)/2;
        idx_x = floor (dx);
        idx_y = floor (dy);
        
        // pivot grid points 
        if (dx-idx_x >= 0.5) {
            x[0] = (double)idx_x + 1;
            x[1] = (double)idx_x;
            x[2] = (double)idx_x + 1;
        } else {
            x[0] = (double)idx_x;
            x[1] = (double)idx_x + 1;
            x[2] = (double)idx_x;
        }
        
        if (dy-idx_y >= (0.5)) {
            y[0] = (double)idx_y + 1;
            y[1] = (double)idx_y + 1;
            y[2] = (double)idx_y;
        } else{
            y[0] = (double)idx_y;
            y[1] = (double)idx_y;
            y[2] = (double)idx_y + 1;
        }
        
        // barycentric weights
        det = (y[1] - y[2])*(x[0] - x[2]) + (x[2] - x[1])*(y[0] - y[2]);
        lambda1 = ((y[1] - y[2]) * (dx - x[2]) + (x[2] - x[1]) * (dy - y[2]))/det;
        lambda2 = ((y[2] - y[0]) * (dx - x[2]) + (x[0] - x[2]) * (dy - y[2]))/det;
        lambda3 = 1.0 - lambda1 - lambda2;
        
        // find weighted local field magnitudes
        if (isEP){
            p1->EP.first =  lambda1*EP_Table[x[0]][y[0]][0] + lambda2*EP_Table[x[1]][y[1]][0] + lambda3*EP_Table[x[2]][y[2]][0];
            p1->EP.second = lambda1*EP_Table[x[0]][y[0]][1] + lambda2*EP_Table[x[1]][y[1]][1] + lambda3*EP_Table[x[2]][y[2]][1];
        }
        p1->DEP.first =  lambda1*DEP_Table[x[0]][y[0]][0] + lambda2*DEP_Table[x[1]][y[1]][0] + lambda3*DEP_Table[x[2]][y[2]][0];
        p1->DEP.second = lambda1*DEP_Table[x[0]][y[0]][1] + lambda2*DEP_Table[x[1]][y[1]][1] + lambda3*DEP_Table[x[2]][y[2]][1];
        
    }
}

/******************************* OTHER OPERATIONS *******************************/

void BrownianDynamicEngine::OutputTrajectory(ostream& os) {
    os << fixed;
    os << setprecision(2);
    os << setw(7);
    for (int i = 0; i < p.size(); i++) {
        os << i << "\t";
        os << p[i].coord.first/a << "\t";
        os << p[i].coord.second/a << endl;
    }
}

void BrownianDynamicEngine::ReadE(const std::string filename, bool EPFlag) {
    // if the same lookup table is used in the previous run, skip the reading process
    int chksum = 0;
    for (int i = 0; i < filename.length(); i++) chksum += (int) filename[i];
    if (prevCheckSum == chksum) return;
    else prevCheckSum = chksum;
    ifstream is;
    is.open(filename.c_str());
    assert(is.is_open());
    string line;
    double dum;
    for (int i = 0; i < E_tot; i++) {
        for (int j = 0; j < E_tot; j++) {
            getline(is, line);
            stringstream linestream(line);
            linestream >> dum; // x
            linestream >> dum; // y
            if (EPFlag){
                linestream >> EP_Table[j][i][0];
                linestream >> EP_Table[j][i][1];
            } else{
                linestream >> dum; // Ex
                linestream >> dum; // Ey
                linestream >> DEP_Table[j][i][0];
                linestream >> DEP_Table[j][i][1];
            }
        }
    }
    is.close();
}

void BrownianDynamicEngine::ReadSavedConfig(string filename) 
{    
    ifstream is;
    string line;
    double dum, x, y;

    is.open(filename.c_str());
    assert(is.is_open());
    while(std::getline(is,line)) {
        stringstream  linestream(line);
        linestream >> dum;
        linestream >> x;
        linestream >> y;
        p.push_back(BD_Particle(std::make_pair(x*a,y*a)));
    }
    is.close();
}

void BrownianDynamicEngine::ReadDiffusivity() 
{
    const int rgbin(25), rbin(50); // rg bin step size
    double dum;
    string line, fileName("./library/Diffusivity.txt");
    ifstream is;
    
    is.open(fileName.c_str());
    assert(is.is_open());
    for (int i = 0; i < rbin * rgbin; i++) {
        getline(is, line);
        std::stringstream linestream(line);
        linestream >> dum;
        linestream >> dum;
        linestream >> dum;
        linestream >> dum;
        linestream >> *(&(dssarray[0][0])+i);
        linestream >> *(&(dsscount[0][0])+i);
    }
    is.close();
}

void BrownianDynamicEngine::GenerateRandomConfig()
{ 
    std::set<int> hist;
    int randMeshGrid = 178;
    int temp;
    for (int i = 0; i < defaultParticleNumber; ++i)
    {
        do {
            temp = (int) floor(rand_uniform(gen) * randMeshGrid * randMeshGrid);
        } while (hist.find(temp) != hist.end());
        hist.insert(temp);
        double detX = (temp/randMeshGrid-randMeshGrid/2+rand_uniform(gen)*2-1)*1000;
        double detY = (temp%randMeshGrid-randMeshGrid/2+rand_uniform(gen)*2-1)*1000;
        p.push_back(BD_Particle(std::make_pair(detX,detY)));
    }
    for (int i = 0; i < 100; i++){
        for (auto p1 = p.begin(); p1 != p.end(); ++p1) {
            for (auto p2 = std::next(p1); p2 != p.end(); ++p2){
                double rijsep = Distance2D(p1,p2);
                if (rijsep < 2*a + 40){
                    double xij = p2->coord.first - p1->coord.first;
                    double yij = p2->coord.second - p1->coord.second;
                    double dx = xij * (2*a + 40 - rijsep)  /sqrt(xij * xij + yij * yij);
                    double dy = yij * (2*a + 40 - rijsep)  /sqrt(xij * xij + yij * yij);
                    p1->coord.first -= dx;
                    p1->coord.second -= dy;
                    p2->coord.first += dx;
                    p2->coord.second += dy;
                }
            }
        }
    }
}

void BrownianDynamicEngine::CalDss() 
{
    const double rgdsmin = 22250;
    const double delrgdsmin = -250;
    const double dssmin = 0.15;
    const double dssmax = 0.50;
    const double distmin = 0.0;
    const double deldist = 1400;
    // calculate rg
    double xmean(0), ymean(0), rgmean(0), rg;
    for (int i = 0; i < p.size(); i++) {
        xmean += p[i].coord.first/p.size();
        ymean += p[i].coord.second/p.size();
    }

    for (int i = 0; i < p.size(); i++) {
        rgmean += (p[i].coord.first - xmean) * (p[i].coord.first - xmean);
        rgmean += (p[i].coord.second - ymean) * (p[i].coord.second - ymean);
    }
    rg = sqrt(rgmean/p.size());
    int rgbinindex = (int) ((rg - rgdsmin) / delrgdsmin) + 1;
    if (rgbinindex <= 0) {rgbinindex = 1;}
    
    // approximate diffusivity by rg and partilce position
    for (int i = 0; i < p.size(); i++) {
        double disttemp = sqrt(pow(p[i].coord.first - xmean, 2) + pow(p[i].coord.second - ymean, 2));
        int distbinindex = (int) ((disttemp - distmin) / deldist) + 1;
        if (rgbinindex >= 1 && rgbinindex <= rgdssbin) {
            if (distbinindex >= 1 && distbinindex <= distdssbin) {
                if (dsscount[rgbinindex-1][distbinindex-1] >= 1) {
                    p[i].D  = dssarray[rgbinindex - 1][distbinindex - 1];
                } else {
                    p[i].D  = dssmax;
                }
            } else {
                p[i].D = dssmax;
            } 
        } else if (rgbinindex >= rgdssbin) {
            p[i].D = dssmin;
        }
    }
}