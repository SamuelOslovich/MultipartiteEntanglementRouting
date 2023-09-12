#ifndef NODE
#define NODE

#include<vector>
#include<memory>
#include<algorithm>
#include<iostream>

#include <boost/heap/fibonacci_heap.hpp>

#include "dist.h"


template <class T> struct p_greater {
  bool operator() (const T& x, const T& y) const { return *x>*y; }
};

class node{
public:
	node();
	node(const node&);
	~node();

  //Constructing the network
  void set_label(int l) { label = l; };
  int get_label() { return label; };
  void add(std::shared_ptr<node>, dist);
  void print(bool=false);

  //Using in the algorithm
  bool operator>(const node&);
  operator bool() { return visited; }

  //Sets
  void set_priority_and_prev(dist, std::shared_ptr<node>);
  void set_visited(bool flag=true) { visited = flag; }
  void set_p_swap(double);

  //Gets
  std::vector<std::vector<std::shared_ptr<node> > > get_prev() {return prev;};
  std::vector<std::shared_ptr<node> > get_prev(int i) { return prev[i]; }
  std::shared_ptr<node> get_prev(int i,int j) { return prev[i][j]; }
  dist get_priority(int i) const { return priority[i]; }
  std::vector<dist> get_paths() {return priority;}
  bool get_visited() { return visited; }
  int get_n_paths() const { return npaths; }
  int get_n_neighbours() {return neighbours.size(); }
  dist get_weight(int i) {return weights[i]; }
  std::shared_ptr<node> get_neighbour(int i) {return neighbours[i].lock(); }
  dist dist_between_neighbour(std::shared_ptr<node>); //if are neighbour return the dist between

  
  //Algorithm
  void add_paths(std::shared_ptr<node>,int,dist,bool&,bool);
  void remove_path(int);
  void update_path(int);
  void update_list_of_paths();

  boost::heap::fibonacci_heap<std::shared_ptr<node>,boost::heap::compare<p_greater<std::shared_ptr<node> > > >::handle_type handle_heap;
 
private:
	unsigned int label;
	std::vector<std::weak_ptr<node> > neighbours; //# number of neighbours
	std::vector<dist> weights; //# number of neighbours
  bool visited;

  int npaths;
	std::vector<dist> priority; //# number of paths
	std::vector<std::vector<std::shared_ptr<node> > > prev; //# number of paths

};

#endif