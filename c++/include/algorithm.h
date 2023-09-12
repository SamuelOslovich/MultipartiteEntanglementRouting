
#ifndef ALGORITHM1
#define ALGORITHM1

#include<vector>
#include<memory>
#include<algorithm>
#include<iostream>
#include <memory>
#include <chrono>
#include <random>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <boost/heap/fibonacci_heap.hpp>

#include "dist.h"
#include "node.h"
#include "network.h"
#include "tree.h"


//Auxiliary
std::shared_ptr<node> convert_ptr(std::shared_ptr<node> n1,std::vector<std::shared_ptr<node> >& net1,std::vector<std::shared_ptr<node> >& net2);
std::vector<std::shared_ptr<node> > convert_vec_ptr(std::vector<std::shared_ptr<node> > n1,std::vector<std::shared_ptr<node> >& net1,std::vector<std::shared_ptr<node> >& net2);
std::vector<std::shared_ptr<node> > reconstruct_path(std::shared_ptr<node> n1, int np,std::vector<std::shared_ptr<node> >& net);
dist lower_bound(std::vector<tree> trees);
dist upper_bound(std::vector<tree> trees);
bool allnetshavepaths(const std::vector<std::vector<std::shared_ptr<node> > >& nets,int k);
void remove_paths(std::vector<std::vector<std::shared_ptr<node> > >& nets,dist trunc);
std::vector<std::vector<int> > combinations(std::vector<int> np); 
bool dominates(std::vector<tree>& list, tree t);


//Algorithms for shortest multi-objective paths and stars
void shortest_path(std::vector<std::shared_ptr<node> >& nodes, std::shared_ptr<node> origin, dist trunc,bool flag_keep_track);
void shortest_path_3(std::vector<std::shared_ptr<node> >& nodes, std::shared_ptr<node> origin, dist trunc,bool flag_keep_track);
std::vector<tree> shortest_star_simple(std::vector<std::shared_ptr<node> >& net, std::vector<int> terminal, dist trunc_path, dist trunc_star);
std::vector<tree> shortest_star(std::vector<std::shared_ptr<node> >& net, std::vector<int> terminal, dist trunc_path, dist trunc_star, std::ostream & pathdata);
std::vector<tree> shortest_star(std::vector<std::shared_ptr<node> >& net, std::vector<int> terminal, dist trunc_path, dist trunc_star);



#endif

