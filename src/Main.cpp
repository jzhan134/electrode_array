/*
 * starting configuration:
 *      file_name - use a specific file
 *      library - use one of the configurations from library folder
 *      random - randomize particles into dilute phase
 */
#include "BrownianDynamicEngine.hpp"

using namespace std;

void controlWithElectrodeArray(BrownianDynamicEngine* model);
void equalizeClusterSize(BrownianDynamicEngine* model);

const int num_thread(1);
int totExecCycle(1), execCycle(0);
std::string fileName;
int targetNum; 

int main(int argc,char* argv[]) {
    if (argc == 3){
        fileName = argv[1];
        targetNum = atoi(argv[2]);
    }
    thread *threads = new thread[num_thread];
    for (int idx = 0; idx < min(num_thread, totExecCycle); idx++){
        BrownianDynamicEngine* model = new BrownianDynamicEngine();
        // threads[idx] = thread(controlWithElectrodeArray, model);
       threads[idx] = thread(equalizeClusterSize, model);
    }
    
    for (int idx = 0; idx < num_thread && threads[idx].joinable(); idx++){
        threads[idx].join();
    }
    delete[] threads;
    return 0;
}

void controlWithElectrodeArray(BrownianDynamicEngine* model){
    model->initialization("fluid.txt");
    model->coreBD("./library/fields/array/3x3_frame.txt", 0, 2);
    model->DEP_quench("./library/fields/array/3x3_frame.txt", 1.75, 20);
    }


void equalizeClusterSize(BrownianDynamicEngine* model){
    model->initialization(fileName+".txt");
    while (model->countNum(1, 1) > targetNum){
        model->EP_move(targetNum);
        model->DEP_quench("./library/fields/array/3x3_frame.txt", 1.75, 0.3);
    }
    model->DEP_quench("./library/fields/array/3x3_frame.txt", 2.5, 2);
    model->DEP_quench("./library/fields/array/3x3_frame.txt", 4, 2);
    model->DEP_quench("./library/fields/array/3x3_frame.txt", 1.75, 30);
}