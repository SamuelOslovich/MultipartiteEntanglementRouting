#include<vector>
#include<memory>
#include<algorithm>
#include<iostream>
#include <chrono>
#include <iterator>
#include <math.h>
#include <cstdio>



#include <boost/heap/fibonacci_heap.hpp>

#include "dist.h"
#include "node.h"
#include "network.h"
#include "algorithm.h"
#include "tree.h"


using HEAP = boost::heap::fibonacci_heap<std::shared_ptr<node>,boost::heap::compare<p_greater<std::shared_ptr<node> > > >;

//--------------------------------------------- AUXILIARY FUNCTIONS -------------------------------------------------------------------------------------

std::shared_ptr<node> convert_ptr(std::shared_ptr<node> n1,std::vector<std::shared_ptr<node> >& net1,std::vector<std::shared_ptr<node> >& net2){
  //std::cout << "This convert_ptr only works for copies of the same network!" << std::endl;
  bool belongsto = false;
  int iaux;
  for( int i=0 ; i<net1.size() ; ++i){
    if( n1 == net1[i] ){
      belongsto = true;
      iaux = i
      ;
    }
  }

  if(belongsto)
    return net2[iaux];
  else if(n1 == nullptr)
    return nullptr;
  else {
    std::cout << "The node doesn't belong to the first network! Returning nullptr! " << std::endl;
    return nullptr;
  }
}

std::vector<std::shared_ptr<node> > convert_vec_ptr(std::vector<std::shared_ptr<node> > n1,std::vector<std::shared_ptr<node> >& net1,std::vector<std::shared_ptr<node> >& net2){
  std::vector<std::shared_ptr<node> > res;
  for( int i=0 ; i<n1.size() ; ++i)
    res.push_back(convert_ptr(n1[i],net1,net2));
  return res;
}

std::vector<std::shared_ptr<node> > reconstruct_path(std::shared_ptr<node> n1, int np,std::vector<std::shared_ptr<node> >& net){
  auto ptr = n1->get_prev(np,0);
  std::vector<std::shared_ptr<node> > pathfinal;

  pathfinal.push_back(n1);

  while(ptr) {
    pathfinal.push_back(ptr);
    auto aux = ptr->get_prev(0);
    ptr = aux[0];
  }

  return pathfinal;  
}

dist lower_bound(std::vector<tree> trees){
  std::vector<dist> truncs;
  for( int i=0; i<trees.size() ; ++i)
    truncs.push_back(trees[i].tree_dist);

  double fid_aux = truncs[0].pathfidelity;
  double p_end_aux = truncs[0].p_end;
  double t_wait_aux = truncs[0].t_wait;
  double t_mem_aux = truncs[0].t_mem;

  for( int i=1 ; i< truncs.size() ; ++i){
    if( truncs[i].pathfidelity < fid_aux)
      fid_aux = truncs[i].pathfidelity;
    if( truncs[i].p_end < p_end_aux)
      p_end_aux = truncs[i].p_end;
    if( truncs[i].t_wait > t_wait_aux)
      t_wait_aux = truncs[i].t_wait;
    if( truncs[i].t_mem < t_mem_aux)
      t_mem_aux = truncs[i].t_mem;
  }

  return dist(fid_aux,1,1,p_end_aux,t_wait_aux,t_mem_aux);
}

dist upper_bound(std::vector<tree> trees){
  std::vector<dist> truncs;
  for( int i=0; i<trees.size() ; ++i)
    truncs.push_back(trees[i].tree_dist);

  double fid_aux = truncs[0].pathfidelity;
  double p_end_aux = truncs[0].p_end;
  double t_wait_aux = truncs[0].t_wait;
  double t_mem_aux = truncs[0].t_mem;

  for( int i=1 ; i< truncs.size() ; ++i){
    if( truncs[i].pathfidelity > fid_aux)
      fid_aux = truncs[i].pathfidelity;
    if( truncs[i].p_end > p_end_aux)
      p_end_aux = truncs[i].p_end;
    if( truncs[i].t_wait < t_wait_aux)
      t_wait_aux = truncs[i].t_wait;
    if( truncs[i].t_mem > t_mem_aux)
      t_mem_aux = truncs[i].t_mem;
  }

  return dist(fid_aux,1,1,p_end_aux,t_wait_aux,t_mem_aux);
}

bool allnetshavepaths(const std::vector<std::vector<std::shared_ptr<node> > >& nets,int k){
  for(int i=nets.size()-1 ; i>=0 ; --i)
    if(nets[i][k]->get_n_paths() == 0)
      return false;
    
  return true;
}

void remove_paths(std::shared_ptr<node> n,dist trunc){
  for( int i=n->get_n_paths()-1 ; i>=0 ; --i){
    auto path = n->get_priority(i);
    auto star = dist(path.pathfidelity,0,0,path.t_wait/path.p_end,0,HUGE_VAL);
    if(!(star < trunc))
      n->remove_path(i);
  }
}

std::vector<std::vector<int> > combinations(std::vector<int> np) { 
    // number of arrays 

    std::vector<std::vector<int> > arr;
    

    for( int i=0 ; i<np.size() ; ++i){
      std::vector<int> aux;
      for( int j=0 ; j<np[i] ; ++j)
        aux.push_back(j);
      arr.push_back(aux);
    }


    int n = arr.size(); 

    std::vector<std::vector<int> > sol;
  
    // to keep track of next element in each of 
    // the n arrays 
    int* indices = new int[n]; 
  
    // initialize with first element's index 
    for (int i = 0; i < n; i++) 
        indices[i] = 0;  

    while (1) { 

        std::vector<int> aux;
        // print current combination 
        for (int i = 0; i < n; i++) 
            aux.push_back(arr[i][indices[i]]);
        sol.push_back(aux);
  
        // find the rightmost array that has more 
        // elements left after the current element  
        // in that array 
        int next = n - 1; 
        while (next >= 0 &&  
              (indices[next] + 1 >= arr[next].size())) 
            next--; 
  
        // no such array is found so no more  
        // combinations left 
        if (next < 0) 
            break; 
  
        // if found move to next element in that  
        // array 
        indices[next]++; 
  
        // for all arrays to the right of this  
        // array current index again points to  
        // first element 
        for (int i = next + 1; i < n; i++) 
            indices[i] = 0; 
    } 
    return sol;
} 

bool dominates(std::vector<tree>& list, tree t){
  bool flagdominate;
  for(int it=list.size()-1 ; it >= 0 ; --it){
    //some path dominates the new one, immediately disregard
    if(list[it].D(t)){
      flagdominate = false;
      break;
    }
    //new path dominates one or more of the paths, erase them and add the new one
    else if( t.D(list[it])){
      list.erase(list.begin()+it);
      flagdominate = true;
    }
    //new path doesn't dominate nor is dominated, add to the list
    else {
      flagdominate = true;
    } 
  }
  return flagdominate;
}

void add_solution(std::vector<tree>& sol, dist trunc_star, tree possible_solution){
  if(!(possible_solution.tree_dist < trunc_star))
    return;
  else if(sol.size() == 0)
    sol.push_back(possible_solution);
  else {
    if(dominates(sol,possible_solution))
      sol.push_back(possible_solution); 
  }
}

void find_potential_solutions(std::vector<std::vector<std::shared_ptr<node> > >& nets, int n, std::vector<int> terminal, std::vector<tree>& sol, dist trunc_star ){
  
  int t_max = terminal.size();

  std::vector<dist>  v_dist[t_max];
  std::vector<std::vector<std::shared_ptr<node> > >  v_prev[t_max];
  std::vector<dist>::iterator  it_dist[t_max];
  std::vector<std::vector<std::shared_ptr<node> > >::iterator  it_prev[t_max];

  for( int j=0 ; j<t_max ; ++j){
    v_dist[j] = nets[j][n]->get_paths();
    v_prev[j] = nets[j][n]->get_prev();
    it_dist[j] = v_dist[j].begin();
    it_prev[j] = v_prev[j].begin();
  }

  while (it_dist[0] != v_dist[0].end()) {

    // process the pointed-to elements

    std::vector<dist> aux_dist;
    std::vector<std::vector<std::shared_ptr<node> > > aux_prev;

    for( int j=0 ; j<terminal.size() ; ++j){
      aux_dist.push_back(*it_dist[j]);
      aux_prev.push_back(*it_prev[j]);
    }

    tree possible_solution = tree(aux_prev,aux_dist);
    //possible_solution.print(false);

    add_solution(sol,trunc_star,possible_solution);

    // the following increments the "odometer" by 1
    ++it_dist[terminal.size()-1];
    ++it_prev[terminal.size()-1];
    for (int i1 = terminal.size()-1; (i1 > 0) && (it_dist[i1] == v_dist[i1].end()); --i1) {
      it_dist[i1] = v_dist[i1].begin();
      it_prev[i1] = v_prev[i1].begin();
      ++it_dist[i1-1];
      ++it_prev[i1-1];
    }
  }
}

void find_potential_solutions(std::vector<std::shared_ptr<node> >& net, int n, std::vector<int> terminal, std::vector<tree>& sol, dist trunc_star){
  int t_max = terminal.size();

  std::vector<dist>  v_dist[t_max];
  std::vector<std::vector<std::shared_ptr<node> > >  v_prev[t_max];
  std::vector<dist>::iterator  it_dist[t_max];
  std::vector<std::vector<std::shared_ptr<node> > >::iterator  it_prev[t_max];

  for( int j=0 ; j<t_max ; ++j){
    v_dist[j] = net[j]->get_paths();
    v_prev[j] = net[j]->get_prev();
    it_dist[j] = v_dist[j].begin();
    it_prev[j] = v_prev[j].begin();
  }

  while (it_dist[0] != v_dist[0].end()) {

    // process the pointed-to elements

    std::vector<dist> aux_dist;
    std::vector<std::vector<std::shared_ptr<node> > > aux_prev;

    for( int j=0 ; j<terminal.size() ; ++j){
      aux_dist.push_back(*it_dist[j]);
      aux_prev.push_back(*it_prev[j]);
    }

    tree possible_solution = tree(aux_prev,aux_dist);
    //possible_solution.print(false);

    add_solution(sol,trunc_star,possible_solution);

    // the following increments the "odometer" by 1
    ++it_dist[terminal.size()-1];
    ++it_prev[terminal.size()-1];
    for (int i1 = terminal.size()-1; (i1 > 0) && (it_dist[i1] == v_dist[i1].end()); --i1) {
      it_dist[i1] = v_dist[i1].begin();
      it_prev[i1] = v_prev[i1].begin();
      ++it_dist[i1-1];
      ++it_prev[i1-1];
    }
  }
}

void simplify_paths(std::shared_ptr<node> n){
  //n->print();
  auto nmax = n->get_n_paths();
  for( int i=nmax-1 ; i>=0 ; --i)
    n->update_path(i);
  //n->print();
  if(n->get_n_paths() == 0) return;
  n->update_list_of_paths();
  //n->print();
}

//-------------------------------- MAIN ALGORITHMS FOR SHORTEST STAR ------------------------------------------------------------------------------------------------------------

void shortest_path(std::vector<std::shared_ptr<node> >& nodes, std::shared_ptr<node> origin, dist trunc, bool flag_keep_track){
  origin->set_priority_and_prev(dist(),nullptr);

  HEAP Q;
  origin->handle_heap = Q.push(origin);

  while(!Q.empty()){
    auto selected_node = Q.top();

    //Remove visited nodes
    selected_node->set_visited();
    Q.pop();

    //Go to each neighbour and write the path in each node
    for(int i=0 ; i < selected_node->get_n_neighbours() ; ++i){
      
      if(selected_node->get_n_paths()==0 || selected_node->get_n_neighbours() == 1) continue ;
      auto neighbour = selected_node->get_neighbour(i);
      
      bool flagadd = false;

      //If neighbour wasn't checked and not visited (paths=0) and there are paths to there, add to the queue
      if(neighbour->get_n_paths()==0){
        neighbour->add_paths(selected_node,i,trunc,flagadd,flag_keep_track);
        if(!flagadd) continue;
        neighbour->handle_heap = Q.push(neighbour);
      }
      //If neighbour has been checked or visited, add paths and update him on queue
      else if(neighbour->get_n_paths()!=0){
        if(*neighbour){
          neighbour->add_paths(selected_node,i,trunc,flagadd,flag_keep_track);
          if(flagadd){
            neighbour->set_visited(false);
            neighbour->handle_heap = Q.push(neighbour);
          }
        }
        else{
          neighbour->add_paths(selected_node,i,trunc,flagadd,flag_keep_track);
          Q.update(neighbour->handle_heap);
        }
      }
    }
  }
}

void shortest_path_3(std::vector<std::shared_ptr<node> >& nodes, std::shared_ptr<node> origin, dist trunc, bool flag_keep_track){
  origin->set_priority_and_prev(dist(),nullptr);

  HEAP Q;
  origin->handle_heap = Q.push(origin);

  while(!Q.empty()){
    auto selected_node = Q.top();

    //std::cout << "Visiting node: " << selected_node->get_label() << std::endl;

    //Remove visited nodes
    selected_node->set_visited();
    Q.pop();

    //Go to each neighbour and write the path in each node
    for(int i=0 ; i < selected_node->get_n_neighbours() ; ++i){
      
      if(selected_node->get_n_paths()==0) continue ;
      auto neighbour = selected_node->get_neighbour(i);
      
      bool flagadd = false;

      //If neighbour wasn't checked and not visited (paths=0) and there are paths to there, add to the queue
      if(neighbour->get_n_paths()==0){
        neighbour->add_paths(selected_node,i,trunc,flagadd,flag_keep_track);
        if(!flagadd) continue;
        neighbour->handle_heap = Q.push(neighbour);
      }
      //If neighbour has been checked or visited, add paths and update him on queue
      else if(neighbour->get_n_paths()!=0){
        if(*neighbour){
          neighbour->add_paths(selected_node,i,trunc,flagadd,flag_keep_track);
          if(flagadd){
            neighbour->set_visited(false);
            neighbour->handle_heap = Q.push(neighbour);
          }
        }
        else{
          neighbour->add_paths(selected_node,i,trunc,flagadd,flag_keep_track);
          Q.update(neighbour->handle_heap);
        }
      }
    }
  }

  //std::cout << "eliminating extra paths" << std::endl;

  for( int i=0 ; i<nodes.size() ; ++i){
    //nodes[i]->print();
    if(nodes[i]->get_n_paths() == 0)
      continue;
    simplify_paths(nodes[i]);
    //nodes[i]->print();
  }
}

std::vector<tree> shortest_star(std::vector<std::shared_ptr<node> >& net, std::vector<int> terminal, dist trunc_path, dist trunc_star){
  std::vector<std::vector<std::shared_ptr<node> > > nets;

  dist trunc = trunc_path;
  std::vector<tree> sol;

  for( int i=0 ; i< terminal.size() ; ++i){
    //Solving the shortest path algorithm for every terminal node
    nets.push_back(copy_network(net));
    shortest_path_3(nets[i],nets[i][terminal[i]],trunc_path,true);

    for( int j=0 ; j<terminal.size() ; ++j){
      if(j==i) //don't search in the terminal of each round
        continue;
      else if(nets[i][terminal[j]]->get_n_paths() ==0){
        std::cout << "No solutions, no path connecting two terminals" << std::endl;
        return sol;
      }
    }
  }

  //Finding all possible choices for the center node, excluding terminals already analised in the begining and already part of the solution
  std::vector<int> nodeslist;

  for(int i=0 ; i<net.size() ; ++i)
    if( allnetshavepaths(nets,i) && net[i]->get_n_neighbours() >= terminal.size())
      nodeslist.push_back(i);

  if(nodeslist.size() == 0){
    return sol;
  }
  else{
    for( int i=0 ; i<nodeslist.size() ; ++i){

      int t_max = terminal.size();

      std::vector<dist>  v_dist[t_max];
      std::vector<std::vector<std::shared_ptr<node> > >  v_prev[t_max];
      std::vector<dist>::iterator  it_dist[t_max];
      std::vector<std::vector<std::shared_ptr<node> > >::iterator  it_prev[t_max];

      for( int j=0 ; j<t_max ; ++j){
        v_dist[j] = nets[j][nodeslist[i]]->get_paths();
        v_prev[j] = nets[j][nodeslist[i]]->get_prev();
        it_dist[j] = v_dist[j].begin();
        it_prev[j] = v_prev[j].begin();
      }

      while (it_dist[0] != v_dist[0].end()) {

        // process the pointed-to elements

        std::vector<dist> aux_dist;
        std::vector<std::vector<std::shared_ptr<node> > > aux_prev;

        for( int j=0 ; j<terminal.size() ; ++j){
          aux_dist.push_back(*it_dist[j]);
          aux_prev.push_back(*it_prev[j]);
        }

        tree possible_solution = tree(aux_prev,aux_dist);
        //possible_solution.print(false);

        add_solution(sol,trunc_star,possible_solution);

        // the following increments the "odometer" by 1
        ++it_dist[terminal.size()-1];
        ++it_prev[terminal.size()-1];
        for (int i1 = terminal.size()-1; (i1 > 0) && (it_dist[i1] == v_dist[i1].end()); --i1) {
          it_dist[i1] = v_dist[i1].begin();
          it_prev[i1] = v_prev[i1].begin();
          ++it_dist[i1-1];
          ++it_prev[i1-1];
        }
      }
    }
  }
  return sol;
}

std::vector<tree> shortest_star(std::vector<std::shared_ptr<node> >& net, std::vector<int> terminal, dist trunc_path, dist trunc_star,std::ostream & pathdata){
  std::vector<std::vector<std::shared_ptr<node> > > nets;

  dist trunc = trunc_path;
  std::vector<tree> sol;


  for( int i=0 ; i< terminal.size() ; ++i){
    //Solving the shortest path algorithm for every terminal node ----------
    nets.push_back(copy_network(net));
    clock_t start, end;
    start = clock();
    std::cout << "Running shortest path!" << std::endl;
    shortest_path_3(nets[i],nets[i][terminal[i]],trunc_path,true);
    end = clock();

    //Write the paths in auxiliary file! -----------------------------------
    int count =0;
    int npathsall = 0;
    int N_connected = 0;
    int N = nets[i].size();

    for( int j=0 ; j<nets[i].size() ; ++j){
        if( nets[i][j]->get_n_paths() != 0){
            ++count;
            npathsall += nets[i][j]->get_n_paths();
        }
        if( nets[i][j]->get_n_neighbours() != 0)
            ++N_connected;
        //auxnet[i]->print(false);
    }
    std::cout << N << " " << (double) end-start << " " << (double) count/(N_connected) << " " << npathsall / (double) (count) << std::endl;
    pathdata << N << " " << end-start << " " << (double) count/(N_connected) << " " << npathsall / (double) (count) << std::endl;

    //Verify necessaary condiition if there might exist a solution ----------
    for( int j=0 ; j<terminal.size() ; ++j){
      if(j==i) //don't search in the terminal of each round
        continue;
      else if(nets[i][terminal[j]]->get_n_paths() ==0){
        std::cout << "No solutions, no path connecting two terminals" << std::endl;
        return sol;
      }
    }
  }

  //Finding all possible choices for the center node, excluding terminals already analised in the begining and already part of the solution
  std::vector<int> nodeslist;

  for(int i=0 ; i<net.size() ; ++i)
    if( allnetshavepaths(nets,i) && net[i]->get_n_neighbours() >= terminal.size())
      nodeslist.push_back(i);

  if(nodeslist.size() == 0){
    return sol;
  }
  else{
    for( int i=0 ; i<nodeslist.size() ; ++i){

      int t_max = terminal.size();

      std::vector<dist>  v_dist[t_max];
      std::vector<std::vector<std::shared_ptr<node> > >  v_prev[t_max];
      std::vector<dist>::iterator  it_dist[t_max];
      std::vector<std::vector<std::shared_ptr<node> > >::iterator  it_prev[t_max];

      for( int j=0 ; j<t_max ; ++j){
        v_dist[j] = nets[j][nodeslist[i]]->get_paths();
        v_prev[j] = nets[j][nodeslist[i]]->get_prev();
        it_dist[j] = v_dist[j].begin();
        it_prev[j] = v_prev[j].begin();
      }

      while (it_dist[0] != v_dist[0].end()) {

        // process the pointed-to elements

        std::vector<dist> aux_dist;
        std::vector<std::vector<std::shared_ptr<node> > > aux_prev;

        for( int j=0 ; j<terminal.size() ; ++j){
          aux_dist.push_back(*it_dist[j]);
          aux_prev.push_back(*it_prev[j]);
        }

        tree possible_solution = tree(aux_prev,aux_dist);
        //possible_solution.print(false);

        add_solution(sol,trunc_star,possible_solution);

        // the following increments the "odometer" by 1
        ++it_dist[terminal.size()-1];
        ++it_prev[terminal.size()-1];
        for (int i1 = terminal.size()-1; (i1 > 0) && (it_dist[i1] == v_dist[i1].end()); --i1) {
          it_dist[i1] = v_dist[i1].begin();
          it_prev[i1] = v_prev[i1].begin();
          ++it_dist[i1-1];
          ++it_prev[i1-1];
        }
      }
    }
  }
  return sol;
}

std::vector<tree> shortest_star_simple(std::vector<std::shared_ptr<node> >& net, std::vector<int> terminal, dist trunc_path, dist trunc_star){
  std::vector<std::vector<std::shared_ptr<node> > > nets;

  std::vector<tree> sol;

 for( int i=0 ; i< terminal.size() ; ++i){
    //Solving the shortest path algorithm for every terminal node
    nets.push_back(copy_network(net));
    shortest_path_3(nets[i],nets[i][terminal[i]],trunc_path,true);

    for( int j=0 ; j<terminal.size() ; ++j){
      if(j==i) //don't search in the terminal of each round
        continue;
      else if(nets[i][terminal[j]]->get_n_paths() ==0){
        std::cout << "No solutions, no path connecting two terminals" << std::endl;
        return sol;
      }
    }
    find_potential_solutions(nets[i],terminal[i],terminal,sol,trunc_star);
  }

  return sol;
  //Finding all possible choices for the center node, excluding terminals already analised in the begining and already part of the solution
  std::vector<int> nodeslist;
  std::vector<int> aux = terminal;

  for(int i=0 ; i<net.size() ; ++i){
    if( aux[0] == i){
      aux.erase(aux.begin());
      continue;
    }
    if( allnetshavepaths(nets,i) && net[i]->get_n_neighbours() >= terminal.size())
      nodeslist.push_back(i);
  }

  if(nodeslist.size() == 0){
    return sol;
  }
  else{
    for( int i=0 ; i<nodeslist.size() ; ++i){

      int t_max = terminal.size();

      std::vector<dist>  v_dist[t_max];
      std::vector<std::vector<std::shared_ptr<node> > >  v_prev[t_max];
      std::vector<dist>::iterator  it_dist[t_max];
      std::vector<std::vector<std::shared_ptr<node> > >::iterator  it_prev[t_max];

      for( int j=0 ; j<t_max ; ++j){
        v_dist[j] = nets[j][nodeslist[i]]->get_paths();
        v_prev[j] = nets[j][nodeslist[i]]->get_prev();
        it_dist[j] = v_dist[j].begin();
        it_prev[j] = v_prev[j].begin();
      }

      while (it_dist[0] != v_dist[0].end()) {

        // process the pointed-to elements

        std::vector<dist> aux_dist;
        std::vector<std::vector<std::shared_ptr<node> > > aux_prev;

        for( int j=0 ; j<terminal.size() ; ++j){
          aux_dist.push_back(*it_dist[j]);
          aux_prev.push_back(*it_prev[j]);
        }

        tree possible_solution = tree(aux_prev,aux_dist);
        //possible_solution.print(false);

        add_solution(sol,trunc_star,possible_solution);

        // the following increments the "odometer" by 1
        ++it_dist[terminal.size()-1];
        ++it_prev[terminal.size()-1];
        for (int i1 = terminal.size()-1; (i1 > 0) && (it_dist[i1] == v_dist[i1].end()); --i1) {
          it_dist[i1] = v_dist[i1].begin();
          it_prev[i1] = v_prev[i1].begin();
          ++it_dist[i1-1];
          ++it_prev[i1-1];
        }
      }
    }
  }
  
  return sol;      
}




