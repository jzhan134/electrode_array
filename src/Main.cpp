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

int main() {
    thread *threads = new thread[num_thread];
    for (int idx = 0; idx < min(num_thread, totExecCycle); idx++){
        BrownianDynamicEngine* model = new BrownianDynamicEngine();
        threads[idx] = thread(controlWithElectrodeArray, model);
    //    threads[idx] = thread(equalizeClusterSize, model);
    }
    
    for (int idx = 0; idx < num_thread && threads[idx].joinable(); idx++){
        threads[idx].join();
    }
    delete[] threads;
    return 0;
}

void controlWithElectrodeArray(BrownianDynamicEngine* model){
    model->initialization("test.txt");
    model->coreBD("./library/fields/array/3x3_frame.txt", 0, 5, false);
    // model->coreBD("./library/fields/array/3x3_frame.txt", 0.5, 2, false);
    // model->coreBD("./library/fields/array/3x3_frame.txt", 1, 2, false);
    // model->coreBD("./library/fields/array/3x3_frame.txt", 1.5, 2, false);
    // model->coreBD("./library/fields/array/3x3_frame.txt", 2, 2, false);
    model->coreBD("./library/fields/array/3x3_frame.txt", 5, 5, false);
    model->coreBD("./library/fields/array/3x3_EP.txt", 0.5, 10, true, 100);
}


void equalizeClusterSize(BrownianDynamicEngine* model){
    model->initialization("test.txt");
    model->coreBD("./library/fields/array/3x3_quad.txt", 0.5, 10, false);
    model->coreBD("./library/fields/array/3x3_EP.txt", 0.5, 10, true, 400);
    model->coreBD("./library/fields/array/3x3_quad.txt", 5, 10, false);
}