#ifndef DIST
#define DIST

#include<vector>
#include<memory>
#include<algorithm>
#include<iostream>

#include <boost/heap/fibonacci_heap.hpp>

class dist{
public:
	dist();
	dist(double,double,double,double,double,double);
	dist(const dist&);
	~dist() {}

	void operator=(const dist&);
	dist operator+(const dist&);
	//dist operator-(const dist&);
	bool operator>(const dist&);  //lexicographic order
	bool operator<(const dist&); //verify if path is possible (if x < trunc, x is possible)
	bool D(const dist&); //dominance relation x.D(y)
	
	void set_dist_p_swap(double p) {p_swap = p;}

	void print() const;  

	//weights
	double pathfidelity;
	double p_gen;
	double p_swap;
	double p_end;
	double t_wait;
	double t_mem;
};



#endif