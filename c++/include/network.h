#ifndef NETWORK
#define NETWORK

#include<vector>
#include<memory>
#include<algorithm>
#include<iostream>

#include <boost/heap/fibonacci_heap.hpp>


#include "dist.h"
#include "node.h"

struct pair_hash {
  template <class T1, class T2>
    std::size_t operator()(std::pair<T1, T2> const &pair) const
  {
    std::size_t h1 = std::hash<T1>()(pair.first);
    std::size_t h2 = std::hash<T2>()(pair.second);

    return h1 ^ h2;
  }
};

void erdos_renyi_model_2(int N, double p, double fmin, double p_gen_min, double p_swap_min, std::vector<std::shared_ptr<node> >& nodes);

void random_geometric(int N, double r, double fmin, double p_gen_min, double p_swap_min, std::vector<std::shared_ptr<node> >& nodes, double tmax=100.);

void squared_lattice_cyclic(int Nx, int Ny, double fmin, double p_gen_min, double p_swap_min, std::vector<std::shared_ptr<node> >& nodes);

void create_grid_from_file(std::string filename, std::vector<std::shared_ptr<node> >& nodes,double fmin, double p_gen_min, double p_swap_min, int n_nodes);

std::vector<std::shared_ptr<node> > copy_network(std::vector<std::shared_ptr<node> >& nodes);

#endif