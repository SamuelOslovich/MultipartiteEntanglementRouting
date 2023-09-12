#include<vector>
#include<memory>
#include<algorithm>
#include<iostream>

#include <boost/heap/fibonacci_heap.hpp>

#include "tree.h"
#include <math.h>   


tree::tree() : tree_struct{}, tree_dist{dist()}, connected{0} {}

tree::tree(std::vector<std::vector<std::shared_ptr<node> > > s,std::vector<dist> pd,dist d, int c) : tree_struct{s}, path_dist{pd}, tree_dist{d}, connected{c} {}

tree::tree(std::vector<std::vector<std::shared_ptr<node> > > s,std::vector<dist> pd) {
	tree_struct = s;
	path_dist = pd;
	tree_dist = maketree(pd);
	connected = pd.size();
}


tree::tree(std::vector<std::shared_ptr<node> > s,dist d) {
	tree_struct.push_back(s);
	path_dist.push_back(d);
	tree_dist = maketree(path_dist);
	connected = 0; 
}

tree::tree(dist d) {
	tree_dist = d;
	connected = 0; 
}

tree::tree(const tree& d){
	tree_struct = d.tree_struct;
	path_dist = d.path_dist;
	tree_dist = d.tree_dist;
	connected = d.connected;
}


void tree::operator=(const tree& d) {
 	tree_struct = d.tree_struct;
 	path_dist = d.path_dist;
	tree_dist = d.tree_dist;
	connected = d.connected;
}

tree tree::operator+(const tree& d) {
	std::vector<std::vector<std::shared_ptr<node> > > ts = tree_struct;
  	std::vector<dist> pd = path_dist;
  	dist td = tree_dist;
  	int c = connected;

	int n_branches = tree_struct.size();

	if(d.connected == 0){
		for( int i=0 ; i<d.tree_struct[0].size() ; ++i)
			ts[n_branches-1].push_back(d.tree_struct[0][i]);
		pd[n_branches-1] = pd[n_branches-1] + d.path_dist[0];
		td = maketree(pd);
	}
	else if(d.connected == 1){
		for( int i=0 ; i<d.tree_struct.size() ; ++i){
			ts.push_back(d.tree_struct[i]);
			pd.push_back(d.path_dist[i]);
		}
		td = maketree(pd);
	}
	return tree(ts,pd,td,c+d.connected);
}

bool tree::operator<(const tree& d) {
	if(tree_dist.pathfidelity < d.tree_dist.pathfidelity && tree_dist.p_end < d.tree_dist.p_end)
		return true;
	else
		return false;
}

bool tree::D(const tree& d) {
	if( (tree_dist.pathfidelity > d.tree_dist.pathfidelity && tree_dist.p_end >= d.tree_dist.p_end) 
		|| (tree_dist.pathfidelity >= d.tree_dist.pathfidelity && tree_dist.p_end > d.tree_dist.p_end))
		return true;
	else if( tree_dist.pathfidelity == d.tree_dist.pathfidelity && tree_dist.p_end == d.tree_dist.p_end )
		return true;
	else return false;
}

void tree::print(bool flag) const {
	std::cout << "Tree connected to " << connected << " with ";
	tree_dist.print();
	if(!flag){
		for( int i=0 ; i<tree_struct.size() ; ++i){
			for( int j=0 ; j<tree_struct[i].size() ; ++j){
				if(tree_struct[i][j])
					std::cout << tree_struct[i][j]->get_label() << " << " ;
			}
			std::cout << std::endl;
		}	
	}
	
}

double ftreefidelity(std::vector<double> fi){

	if(fi.size() ==2){
		double eta = (fi[0]+1.)/2.;
		double gamma = fi[0];

		for( int i=1 ; i<fi.size() ; ++i){
			eta = eta * ((fi[i]+1.)/2.);
			gamma = gamma * fi[i];
		}
		return (eta+gamma)/2;
	}	
	
	if(fi.size() > 2){
		double eta = (fi[0]+1.)/2.;
		double sigma = (1.-fi[0])/2.;
		double gamma = fi[0];

		for( int i=1 ; i<fi.size() ; ++i){
			eta = eta * ((fi[i]+1.)/2.);
			sigma = sigma * ((1.-fi[i])/2.);
			gamma = gamma * fi[i];
		}
		return (eta+sigma+gamma)/2;

	}
}

dist maketree(std::vector<dist> paths){
	std::vector<double> fids;
	fids.push_back(paths[0].pathfidelity);
	double pend = paths[0].p_end;
	double tmax = paths[0].t_wait;
	for( int i=1 ; i<paths.size() ; ++i){
		pend = pend * paths[i].p_end;
		if(paths[i].t_wait>tmax)
			tmax = paths[i].t_wait;
		double dec_factor = exp(-paths[i].t_wait/paths[i].t_mem);
		fids.push_back(paths[i].pathfidelity*dec_factor);
	}

	double pfid = ftreefidelity(fids);
	return dist(pfid,1,1,pend/tmax*1000000,0,HUGE_VAL);	
}



