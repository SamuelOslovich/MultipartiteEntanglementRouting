#include <cmath>
#include <vector>
#include <memory>
#include <random>
#include <chrono>
#include <utility>
#include <unordered_set>
#include <fstream>

#include "network.h"
#include "dist.h"
#include "node.h"

double distance(std::vector<double> pos1, std::vector<double> pos2){
  if(pos1.size() != 2 || pos2.size() !=2 )
    return HUGE_VAL;
  return sqrt(pow(pos1[0]-pos2[0],2)+pow(pos1[1]-pos2[1],2));
}

void erdos_renyi_model_2(int N, double p, double fmin, double p_gen_min, double p_swap_min, std::vector<std::shared_ptr<node> >& nodes) {

  if(nodes.size() != 0) nodes.clear();

  for(int i=0; i<N; i++) {
    nodes.push_back(std::make_shared<node>());
    nodes[i]->set_label(i+1);
  }

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> uniform(0.0,1.0);
  std::uniform_real_distribution<> fid(fmin, 1);
  std::uniform_real_distribution<> pgen(p_gen_min, 1);
  std::uniform_real_distribution<> pswap(p_swap_min, 1);
  std::uniform_real_distribution<> twait(1, 100);
  std::uniform_real_distribution<> tmem(5000, 10000);

  double a;
  double f1;
  double p1;
  double p2;
  double t1;
  double t2;
  for(int i=0; i<N; i++) {
    for(int j=i+1; j<N; j++) {
      a = uniform(generator);
      f1 = fid(generator);
      p1 = pgen(generator);
      t1 = twait(generator);
      t2 = tmem(generator);
      if(a <= p) {
        nodes[i]->add(nodes[j],dist(f1,p1,1,p1,t1,t2));
        nodes[j]->add(nodes[i],dist(f1,p1,1,p1,t1,t2));
      }
    }
  }

  for(int i=0; i<N; i++) {
    p2 = pswap(generator);
    nodes[i]->set_p_swap(p2);
  }
}

void random_geometric(int N, double r, double fmin, double p_gen_min, double p_swap_min, std::vector<std::shared_ptr<node> >& nodes,double tmax) {

  if(nodes.size() != 0) nodes.clear();

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> uniformx(0.0,1.0);
  std::uniform_real_distribution<double> uniformy(0.0,1.0);


  std::vector<std::vector<double> > positions;

  for(int i=0; i<N; i++) {
    nodes.push_back(std::make_shared<node>());
    std::vector<double> posxy;
    posxy.push_back(uniformx(generator));
    posxy.push_back(uniformy(generator));
    positions.push_back(posxy);
    nodes[i]->set_label(i+1);
  }

  std::uniform_real_distribution<double> uniform(0.0,1.0);
  std::uniform_real_distribution<> fid(fmin, 1);
  std::uniform_real_distribution<> pgen(p_gen_min, 1);
  std::uniform_real_distribution<> pswap(p_swap_min, 1);
  std::uniform_real_distribution<> twait(1,tmax);
  std::uniform_real_distribution<> tmem(tmax*sqrt(N)/log(N)*100,tmax*sqrt(N)/log(N)*1000);

  double a;
  double f1;
  double p1;
  double p2;
  double t1;
  double t2;
  for(int i=0; i<N; i++) {
    for(int j=i+1; j<N; j++) {
      a = uniform(generator);
      f1 = fid(generator);
      p1 = pgen(generator);
      t2 = tmem(generator);
      t1 = twait(generator);
      double d=distance(positions[i],positions[j]);
      if(d <= r) {
        nodes[i]->add(nodes[j],dist(f1,p1,1,p1,t1,t2));
        nodes[j]->add(nodes[i],dist(f1,p1,1,p1,t1,t2));
      }
    }
  }

  for(int i=0; i<N; i++) {
    p2 = pswap(generator);
    nodes[i]->set_p_swap(p2);
  }
}

void squared_lattice_cyclic(int Nx, int Ny, double fmin,double p_gen_min, double p_swap_min, std::vector<std::shared_ptr<node> >& nodes) {

  if(nodes.size() != 0) nodes.clear();

  for(int i=0; i<Nx*Ny; i++) {
    nodes.push_back(std::make_shared<node>());
    nodes[i]->set_label(i+1);
  }

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<> fid(fmin, 1);
  std::uniform_real_distribution<> pgen(p_gen_min, 1);
  std::uniform_real_distribution<> pswap(p_swap_min, 1);
  std::uniform_real_distribution<> twait(1, 100);
  std::uniform_real_distribution<> tmem(1000, 10000);

  double f1, f2, f3, f4, p1, p2, p3, p4, t1a, t1b, t1c, t1d, t2a, t2b, t2c, t2d;
  for( int i=0 ; i<Nx*Ny ; ++i ){
    f1 = fid(generator);
    f2 = fid(generator);
    f3 = fid(generator);
    f4 = fid(generator);
    p1 = pgen(generator);
    p2 = pgen(generator);
    p3 = pgen(generator);
    p4 = pgen(generator);
    t1a = twait(generator);
    t1b = twait(generator);
    t1c = twait(generator);
    t1d = twait(generator);
    t2a = tmem(generator);
    t2b = tmem(generator);
    t2c = tmem(generator);
    t2d = tmem(generator);

    //std::cout << p1 << p2 << p3 << p4 << std::endl;
    nodes[i]->add(nodes[(i+Nx)%(Nx*Ny)],dist(f1,p1,1,p1,t1a,t2a));  //all upwards
    nodes[i]->add(nodes[(i-Nx + Nx*Ny)%(Nx*Ny)],dist(f2,p2,1,p2,t1b,t2b)); //all downwards
    nodes[i]->add(nodes[i/Nx*Nx + (i-1 + Ny)%Ny ],dist(f3,p3,1,p3,t1c,t2c)); //all left ones
    nodes[i]->add(nodes[i/Nx*Nx + (i+1 + Ny)%Ny ],dist(f4,p4,1,p4,t1d,t2d)); //all right ones
  }
  for(int i=0; i<Nx*Ny; i++) {
    p2 = pswap(generator);
    nodes[i]->set_p_swap(p2);
  }
}

void create_grid_from_file(std::string filename, std::vector<std::shared_ptr<node> >& nodes, double fmin, double p_gen_min, double p_swap_min, int n_nodes) {
  std::ifstream ifs;
  ifs.open(filename);

  nodes.clear();

  for(int i=0; i<n_nodes; i++) {
    nodes.push_back(std::make_shared<node>());
    nodes[i]->set_label(i+1);
  }

  unsigned int x, y;

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<> fid(fmin, 1);
  std::uniform_real_distribution<> pgen(p_gen_min, 1);
  std::uniform_real_distribution<> pswap(p_swap_min, 1);
  std::uniform_real_distribution<> twait(1, 100);
  std::uniform_real_distribution<> tmem(1000, 10000);

  double f1;
  double p1;
  double p2;
  double t1;
  double t2;

  while(!ifs.eof()) {
    ifs >> x >> y;
    f1 = fid(generator);
    p1 = pgen(generator);
    t1 = twait(generator);
    t2 = tmem(generator);
    //std::cout << x << " " << y << " " << f1 << std::endl;
    nodes[x]->add(nodes[y],dist(f1,p1,1,p1,t1,t2));
    nodes[y]->add(nodes[x],dist(f1,p1,1,p1,t1,t2));
  }

  for(int i=0; i<n_nodes; i++) {
    p2 = pswap(generator);
    nodes[i]->set_p_swap(p2);
  }

}

std::vector<std::shared_ptr<node> > copy_network(std::vector<std::shared_ptr<node> >& nodes){
  std::vector<std::shared_ptr<node> > copy;

  //std::cout << "WATCHOUT! This copy_network only works for labeled 1, 2, 3,.... nodes" << std::endl;
  for( int i=0 ; i< nodes.size() ; ++i){
    copy.push_back(std::make_shared<node>());
    copy[i]->set_label(nodes[i]->get_label());
  }
  for( int i=0 ; i< nodes.size() ; ++i)
    for( int j=0 ; j< nodes[i]->get_n_neighbours() ; ++j)
      copy[i]->add(copy[nodes[i]->get_neighbour(j)->get_label()-1], nodes[i]->get_weight(j));
    
  return copy;
}

