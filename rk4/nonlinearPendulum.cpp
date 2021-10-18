#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;
#define NUM_THREADS 24

typedef long double dtype; //determines floating pt precision
typedef vector<dtype (*)(vector<dtype>)> dfns; //vector of func ptrs

#include "rk4.cpp"

#define PI 3.141592653589793
#define g 9.8
#define L 1.0

//return vector containing every nth element of V
template<class T>
vector<T> nthValues(vector<T> V, int n){
    vector<T> result;
    for(unsigned i = 0; i < V.size(); i++){
        if(i%n==0){result.push_back(V[i]);}
    }
    return result;
}

//print results
template<class T>
void printDV(vector<T> v){printf("time:% 5.5f   vel: % 5.5f   ang: % 5.5f\n",
                                (float)v[0], (float)v[1], (float)v[2]);}

//derivative functions
//V is vector <time, y', y>, at some timestep;
//d2y returns 2nd derivative of y, dy returns derivative of y
dtype d2y(vector<dtype> V){return -1*(g/L)*sin(V[2]);}
dtype dy(vector<dtype> V){return V[1];}

struct thread_params{
    int tid;
    dfns F;
    vector<vector<dtype>> DV;
    dtype step;
    dtype tmax;
};

void *rk4Worker(void *t_arg){
    struct thread_params *p = (struct thread_params*)t_arg;
    vector<vector<dtype>> resultDV = 
        rk4<dtype>((*p).F, (*p).DV, (*p).step, (*p).tmax);

    ofstream resultfile;
    resultfile.open("data/pendulumData"+std::to_string((*p).tid)+".csv");
    for(int j = 0; j < resultDV.size(); j++){
        resultfile << std::to_string(resultDV[j][0]) << ",";
        resultfile << std::to_string(resultDV[j][1]) << ",";
        resultfile << std::to_string(resultDV[j][2]) << "\n";
    }
    resultfile.close();
    pthread_exit(NULL);
}

int main(){
    dtype initStep = 0.01;
    vector<dtype> h = {initStep}; 
    int smallestSplit = 25;
    for(int i = 2; i < smallestSplit; i++){
        h.push_back(h[0] / i);
    }
    dtype tmax = 1000.0;
    
    //DV is vector <time, vel, ang>, at each timestep;
    //DV[0] = initial conditions, angle in radians
    vector<vector<dtype>> inits = {{0.0, 0.0, PI-0.0001}};
    vector<vector<dtype>> DV = inits;

    //F is vector of pointers to derivative functions
    dfns F = {&d2y, &dy};

    //prepare thread parameters
    pthread_t threads[NUM_THREADS];
    thread_params args[NUM_THREADS];
    for(int i = 0; i < NUM_THREADS; i++){
        args[i].tid = i;
        args[i].F = F;
        args[i].DV = DV;
        args[i].step = h[i];
        args[i].tmax = tmax;
    }
    
    //start rk4 threads
    unsigned counter = 0;
    for(int i = 0; i < NUM_THREADS; i++){
        cout << "starting thread " << 1+counter++ << "/" << smallestSplit-1;
        cout <<  " for stepsize: " << h[i] << endl;
        pthread_create(&threads[i], NULL, rk4Worker, (void *)&args[i]);
    }

    for(int i = 0; i < NUM_THREADS; i++){
        pthread_join(threads[i], NULL);
    }
    return 0;
}
