#ifndef TREE
#define TREE

#include<vector>
#include<memory>
#include<algorithm>
#include<iostream>

#include <boost/heap/fibonacci_heap.hpp>

#include "node.h"
#include "dist.h"

class node;

class tree{
public:
	tree();
	tree(std::vector<std::vector<std::shared_ptr<node> > >,std::vector<dist>,dist,int);
	tree(std::vector<std::vector<std::shared_ptr<node> > >,std::vector<dist>);
	tree(std::vector<std::shared_ptr<node> >,dist);
	tree(dist);
	tree(const tree&);
	~tree() {}

	void operator=(const tree&);
	tree operator+(const tree&);
	bool operator<(const tree&);
	bool D(const tree&);

	void print(bool=false) const;

	//Elements
  	std::vector<std::vector<std::shared_ptr<node> > > tree_struct;
  	std::vector<dist> path_dist;
  	dist tree_dist;
  	int connected;
};

dist maketree(std::vector<dist> paths);
double ftreefidelity(std::vector<double> fi);


#endif