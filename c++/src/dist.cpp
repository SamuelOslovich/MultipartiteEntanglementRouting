#include<vector>
#include<memory>
#include<algorithm>
#include<iostream>
#include <stdio.h>
#include <math.h>

#include <boost/heap/fibonacci_heap.hpp>

#include "dist.h"


dist::dist() : pathfidelity{1}, p_gen{1}, p_swap{1}, p_end{1}, t_wait{0.}, t_mem{HUGE_VAL} {}

dist::dist(double f, double p1, double p2,double p3,double t1, double t2) :pathfidelity{f}, p_gen{p1}, p_swap{p2}, p_end{p3}, t_wait{t1}, t_mem{t2} {}

dist::dist(const dist& d){
	pathfidelity = d.pathfidelity;
	p_gen = d.p_gen;
	p_swap = d.p_swap;
	p_end = d.p_end;
	t_wait = d.t_wait;
	t_mem = d.t_mem;
}


void dist::operator=(const dist& d) {
	pathfidelity = d.pathfidelity;
	p_gen = d.p_gen;
	p_swap = d.p_swap;
	p_end = d.p_end;
	t_wait = d.t_wait;
	t_mem = d.t_mem;
}

dist dist::operator+(const dist& d) {
	double p_end_aux = 1;
	if(p_end == 1)
		p_end_aux = p_end * d.p_gen;
	else
		p_end_aux = p_end * d.p_gen * d.p_swap;
	
	double t_mem_aux = 1/(1/t_mem+1/d.t_mem);
	
	return dist(d.pathfidelity*pathfidelity,p_end_aux,1,p_end_aux,t_wait + d.t_wait,t_mem_aux);
}

/*dist dist::operator-(const dist& d) {
	return dist(pathfidelity/d.pathfidelity,hops-d.hops);
}*/

bool dist::operator>(const dist& d) {
	if(pathfidelity < d.pathfidelity)
		return true;
	else if(pathfidelity == d.pathfidelity && p_end < d.p_end)
		return true;
	else if(pathfidelity == d.pathfidelity && p_end == d.p_end && t_wait > d.t_wait)
		return true;
	else if(pathfidelity == d.pathfidelity && p_end == d.p_end && t_wait == d.t_wait && t_mem < d.t_mem)
		return true;
	else
		return false;
}

bool dist::operator<(const dist& d) {
	if(pathfidelity >= d.pathfidelity && p_end >= d.p_end && t_wait <= d.t_wait && t_mem >= d.t_mem)
		return true;
	else
		return false;
}

bool dist::D(const dist& d) {
	if( ( pathfidelity > d.pathfidelity && p_end >= d.p_end && t_wait <= d.t_wait && t_mem >= d.t_mem ) ||
		( pathfidelity >= d.pathfidelity && p_end > d.p_end && t_wait <= d.t_wait && t_mem >= d.t_mem ) ||
		( pathfidelity >= d.pathfidelity && p_end >= d.p_end && t_wait < d.t_wait && t_mem >= d.t_mem ) ||
		( pathfidelity >= d.pathfidelity && p_end >= d.p_end && t_wait <= d.t_wait && t_mem > d.t_mem ))
		return true;
	else if(pathfidelity == d.pathfidelity && p_end == d.p_end && t_wait == d.t_wait && t_mem == d.t_mem )
		return true;
	else return false;
}


void dist::print() const{
	std::cout << "Dist: " << pathfidelity << " || P_END: " << p_end << " || P_SWAP: " << p_swap << " || T_WAIT: " << t_wait << " || T_MEM: " << t_mem << std::endl;
}

