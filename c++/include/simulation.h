#ifndef SIMULATION
#define SIMULATION

#include<vector>
#include<memory>
#include<algorithm>
#include<iostream>

#include <boost/heap/fibonacci_heap.hpp>

#include "dist.h"
#include "node.h"
#include "network.h"
#include "algorithm.h"


//Network type: (1) erdos-renyi (2) square lattice
double average(std::vector<double> values);
double standarddeviation(std::vector<double> values);
double average(std::vector<int> values);
double standarddeviation(std::vector<int> values);
bool not_in(int k,std::vector<int> list);


void scaling(char* filename, int network_type, int N, double alpha_min, double alpha_max, double alpha_step, double fmin, double pmin, int nb);
void simulation_path(char* filename, char* filenamedata, int network_type, int n_min, int n_max, int n_step, int nb, double fmin, double pmin);
void simulation_star(char* filename, char* filenamedata, int network_type, int n_min, int n_max, int n_step, int nb, int nt, double fmin, double pmin);


#endif