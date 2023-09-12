#include<vector>
#include<memory>
#include<algorithm>
#include<iostream>
#include <stdio.h>
#include <math.h>

#include <boost/heap/fibonacci_heap.hpp>

#include "node.h"

using HEAP = boost::heap::fibonacci_heap<std::shared_ptr<node>,boost::heap::compare<p_greater<std::shared_ptr<node> > > >;

node::node() : label{0} , neighbours{}, weights{}, visited{false}, npaths{0}, priority{}, prev{} {}

node::node(const node& n) : label{n.label} , neighbours{n.neighbours}, weights{n.weights}, visited{n.visited}, npaths{n.npaths}, priority{n.priority}, prev{n.prev} {}

node::~node() {}

void node::add(std::shared_ptr<node> prev, dist w) {
  neighbours.push_back(prev);
  weights.push_back(w);
}

bool node::operator>(const node& n1) {
  // Priority queue ordering
  if(npaths != 0 && n1.get_n_paths() != 0) {
    for( int i=0 ; i< n1.get_n_paths() ; ++i)
    	for( int j=0 ; j< npaths ; ++j)
    		if( n1.get_priority(i).D(priority[j])) //n1 is better is it has some path that dominates some path in *this
    			return true;
  }
  else if(npaths == 0) return true;
  else return false;
}

void node::print(bool flagfast) {
	std::cout << "#####################################" << std::endl;

	std::cout << " Label: " << label << std::endl;

	std::cout << "-------------------------------------" << std::endl;

	for( int i=0 ; i<neighbours.size() ; ++i){
		auto weight = weights[i];
		auto nlabel = neighbours[i].lock();
		std::cout << "Neighbour number: " << nlabel->label << " -  with fidelity: " ;
		weight.print();
	}

	std::cout <<  "Number of paths: " << npaths << "  " << std::endl;
	if(npaths!=0){
		if(flagfast){
			std::cout << "Shortest path has dist: ";
			priority[0].print();

			std::cout << label << " > " ;
			
			auto ptr = prev[0][0];

	    	while(ptr) {
	      		std::cout << ptr->get_label();
	      		auto aux = ptr->get_prev(0);
	      		ptr = aux[0];
	      		if(ptr) std::cout << " > ";
	     	 	else std::cout << std::endl;
			}
			std::cout << std::endl;
		}
		else{
			for(int i=0 ; i<npaths ; ++i){
				std::cout << "Path " << i << " :" ;
				priority[i].print();
				for( int j=0 ; j< prev[i].size() ; ++j)
					if(prev[i][j] != nullptr)
						std::cout << prev[i][j]->get_label() << " << ";
				std::cout << label << std::endl;
			}
		}
	}
	else 
	   std::cout << "Impossible to establish the connection." << std::endl;
}

void node::set_priority_and_prev(dist p, std::shared_ptr<node> n) {
	if(npaths != 0) {
		priority.clear();
		prev.clear();
	}
	priority.push_back(p);
	std::vector<std::shared_ptr<node> > aux;
	aux.push_back(n);
	prev.push_back(aux);
	npaths = 1;
}

void node::set_p_swap(double p) {
	for( int i=0 ; i< weights.size() ; ++i)
		weights[i].set_dist_p_swap(p);
}

dist node::dist_between_neighbour(std::shared_ptr<node> n2){
	for( int i=0 ; i<neighbours.size() ; ++i)
		if( neighbours[i].lock() == n2)
			return weights[i];

	std::cout << "Not neighbour!!" << std::endl;
	return dist();
}

void node::add_paths(std::shared_ptr<node> origin, int j, dist trunc,bool& flag,bool flag_keep_track){
	auto possible_paths = origin->get_paths();
	std::vector<dist> selected_paths;
	std::vector<int> selected_pos;
	
	for( int i=0 ; i<possible_paths.size() ; ++i){
		dist faux;
		faux = possible_paths[i] + origin->get_weight(j);
		if( !(faux < trunc))
			continue;

		if(npaths == 0){
			selected_paths.push_back(faux);
			if(flag_keep_track)
				selected_pos.push_back(i);
		}

		else if(npaths != 0){
			bool flagdominate;
			for(int it=npaths-1 ; it >= 0 ; --it){
				//some path dominates the new one, immediately disregard
				if(priority[it].D(faux)){
					flagdominate = false;
					break;
				}
				//new path dominates one or more of the paths, erase them and add the new one
				else if( faux.D(priority[it])){
					priority.erase(priority.begin()+it);
					prev.erase(prev.begin()+it);
					npaths = npaths-1;
					flagdominate = true;
				}
				//new path doesn't dominate nor is dominated, add to the list
				else {
					flagdominate = true;
				}	
			}
			if(flagdominate){
				selected_paths.push_back(faux);
				if(flag_keep_track)
					selected_pos.push_back(i);
			}
		}
	}

	//no paths, add all, dominance inheritance
	if(npaths == 0 && selected_paths.size() != 0){
		npaths = selected_paths.size();
		flag = true;
		priority = selected_paths;
		for( int i=0 ; i<selected_paths.size() ; ++i){
			if(!flag_keep_track){
				std::vector<std::shared_ptr<node> > aux;
				aux.push_back(origin);
				prev.push_back(aux);
			}
			else{
				auto prevaux = origin->get_prev(selected_pos[i]);
				prevaux.push_back(origin);
				prev.push_back(prevaux);
			}
		}
	}

	//some paths, add set of non-dominated ones and remove if necessary
	else if( npaths != 0 && selected_paths.size() != 0){
		npaths = npaths + selected_paths.size();
		flag = true;
		for( int i=0 ; i<selected_paths.size() ; ++i){
			priority.push_back(selected_paths[i]);
			if(!flag_keep_track){
				std::vector<std::shared_ptr<node> > aux;
				aux.push_back(origin);
				prev.push_back(aux);
			}
			else{
				auto prevaux = origin->get_prev(selected_pos[i]);
				prevaux.push_back(origin);
				prev.push_back(prevaux);
			}
			
		}
		
	}	
}

void node::remove_path(int i){
	priority.erase(priority.begin() + i);
	prev.erase(prev.begin()+i);
	npaths = npaths-1;
}

void node::update_path(int i){
	auto d = priority[i];
	double fidnew = d.pathfidelity * exp(-d.t_wait/d.t_mem);
	if(fidnew < 0.33333333){
		this->remove_path(i);
	}
	else
		priority[i] = dist(fidnew,0,0,d.p_end,d.t_wait,HUGE_VAL);
	
}

void node::update_list_of_paths(){
	auto aux_priority = priority;
	auto aux_prev = prev;
	
	priority.clear();
	prev.clear();

	priority.push_back(aux_priority[0]);
	prev.push_back(aux_prev[0]);
	
	npaths=1;

	
	for(int j=npaths-1 ; j >= 0 ; --j){
		for( int i=1 ; i<aux_priority.size(); ++i){
			//some path dominates the new one, immediately disregard
			if(priority[j].D(aux_priority[i]))
				break;
			//new path dominates one or more of the paths, erase them and add the new one
			else if( aux_priority[i].D(priority[j])){
				this->remove_path(j);
				priority.push_back(aux_priority[i]);
				prev.push_back(aux_prev[i]);
				npaths = npaths+1;
			}
			//new path doesn't dominate nor is dominated, add to the list
			else{
				priority.push_back(aux_priority[i]);
				prev.push_back(aux_prev[i]);
				npaths=npaths+1;
			}
		}
	}
	npaths = priority.size();
			
}




