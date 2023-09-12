#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <chrono>
#include <random>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <boost/heap/fibonacci_heap.hpp>

#include "dist.h"
#include "node.h"
#include "network.h"
#include "algorithm.h"
#include "simulation.h"

double average(std::vector<double> values){
	double res = 0;
	for( int i=0 ; i<values.size() ; ++i)
		res = res + values[i];
	return res/values.size();
}
double average(std::vector<int> values){
	double res = 0;
	for( int i=0 ; i<values.size() ; ++i)
		res = res + (double) values[i];
	return res/values.size();
}
double standarddeviation(std::vector<double> values){
	double mean = average(values);
	double res = 0;
	for( int i=0 ; i<values.size() ; ++i)
		res = res + (values[i]-mean)*(values[i]-mean) ;
	return sqrt(res/((double)values.size()-1));
}
double standarddeviation(std::vector<int> values){
	double mean = average(values);
	double res = 0;
	for( int i=0 ; i<values.size() ; ++i)
		res = res + ((double) values[i]-mean)*((double) values[i]-mean);
	return sqrt(res/((double)values.size()-1));
}
bool not_in(int k,std::vector<int> list){
	for( int i=0 ; i<list.size() ; ++i)
		if( k == list[i])
			return false;
	return true;	
}

// Simulations to find the value of alpha such that the network is connected
// alpha in [alpha_min,alpha_max] in alpha_step steps
// N size of the network
// fmin - minimum fidelity; pmin - minimum probability of success; nb - average degree (number of neighbours)
void scaling(char* filename, int network_type, int N, double alpha_min, double alpha_max, double alpha_step, double fmin, double pmin, int nb){
	std::ofstream outfile (filename);
	
	int NNETWORKS = 10;
	int NSIMUL = 10;
	int NNODES = N*pow(N,network_type-1);
	
	
	double alpha = alpha_min;
	while( alpha < alpha_max){
		std::vector<double> time_elapsed_aux;
		std::vector<double> nodes_reached_aux;
		std::vector<double> n_paths_aux;
		for( int iNET=0 ; iNET < NNETWORKS ; ++iNET){
			std::vector<std::shared_ptr<node> > net;
			if( network_type == 1 ){
				double average_neighbours_aux = 0;
				while(average_neighbours_aux < 2){
					average_neighbours_aux = 0;
					erdos_renyi_model_2(N,(double) nb / N,(4*pow(fmin,log(nb)*alpha/log(N))-1)/3, pmin,pmin,net);
					for( int j=0 ; j<net.size() ; ++j){
				        average_neighbours_aux = average_neighbours_aux + (double) net[j]->get_n_neighbours();
				    }
				    average_neighbours_aux = average_neighbours_aux/N;
				}
			}
			else if( network_type == 2 )
				squared_lattice_cyclic( N,N,(4*pow(fmin,alpha/N)-1)/3,pmin,pmin,net);
			else {
				std::cout << "Invalid network type!" << std::endl;
				break;
			}
			for (int iSIM=0 ; iSIM<NSIMUL; ++iSIM){
				std::vector<std::shared_ptr<node> > netcopy;
				netcopy = copy_network(net);
				dist trunc = dist(((4*fmin-1)/3),0,0,0,1000000,0);
				clock_t start, end;
				double elapsed_time;
				//Algorithm		        			
				start = clock();
				shortest_path(netcopy,netcopy[iSIM],trunc,false);
				end = clock();

				elapsed_time = (double)(end - start);
				time_elapsed_aux.push_back(elapsed_time);
				int ActualNNODES = 0;
				int count =0;
        		int npathsall = 0;
        		for( int i=0 ; i<net.size() ; ++i){
        			if( net[i]->get_n_neighbours() > 0)
        				++ActualNNODES;
		            if( netcopy[i]->get_n_paths() != 0){
		                ++count;
		                npathsall += netcopy[i]->get_n_paths();
		            }
		        }
		        nodes_reached_aux.push_back( (double) count / ActualNNODES);
		        n_paths_aux.push_back(npathsall/count);
			}
		}
		std::cout << alpha << "   " << average(time_elapsed_aux) << "   "  << standarddeviation(time_elapsed_aux) << "  " << average(nodes_reached_aux) << "  " << standarddeviation(nodes_reached_aux) << "  " << average(n_paths_aux) << "  " << standarddeviation(n_paths_aux) << std::endl << std::flush;
		outfile << alpha << "   " << average(time_elapsed_aux) << "   "  << standarddeviation(time_elapsed_aux) << "  " << average(nodes_reached_aux) << "  " << standarddeviation(nodes_reached_aux) << "  " << average(n_paths_aux) << "  " << standarddeviation(n_paths_aux) << std::endl;
		alpha = alpha + alpha_step;
	}
	
	outfile.close();	
}

// Simulations of MOSP algorithm 
// filename/filenamedata - containers for storing results
// n_nodes - size of the network in [n_min,n_max] with steps of size n_step
// fmin - minimum fidelity; pmin - minimum probability of success; nb - average degree (number of neighbours)
void simulation_path(char* filename, char* filenamedata, int network_type, int n_min, int n_max, int n_step, int nb, double fmin, double pmin){
	std::ofstream outfile (filename);
	std::ofstream data (filenamedata);

	int aux_n_nodes = (n_max-n_min)/n_step;
	std::vector<int> n_nodes;
	std::vector<double> time_elapsed_MEAN;
	std::vector<double> time_elapsed_STD;
	std::vector<double> nodes_reached_MEAN;
	std::vector<double> nodes_reached_STD;
	std::vector<double> n_neighbours_MEAN;
	std::vector<double> n_neighbours_STD;
	std::vector<double> n_paths_MEAN;
	std::vector<double> n_paths_STD;

	int NNETWORKS = 10;
	int NSIMUL = 20;

	
	for(int i=0 ; i<aux_n_nodes ; ++i)
		n_nodes.push_back(n_min+ i*n_step);

	for( int i=0 ; i<aux_n_nodes ; ++i){
		std::vector<double> time_elapsed_aux;
		std::vector<double> nodes_reached_aux;
		std::vector<double> n_neighbours_aux;
		std::vector<double> n_paths_aux;
		std::cout << n_nodes[i] << std::endl;	
		for( int iNET=0 ; iNET < NNETWORKS ; ++iNET){
			std::vector<std::shared_ptr<node> > net;

			// Erdos-Renyi networks
			if( network_type == 1 ){
				double average_neighbours_aux = 0;
				while(average_neighbours_aux < 2){
					average_neighbours_aux = 0;
					// note minimum fidelity is calculated such that the network is functionally connected
					erdos_renyi_model_2(n_nodes[i],(double) nb / n_nodes[i],(4*pow(fmin,nb*2/3/log(n_nodes[i]))-1)/3, pmin, pmin,net); 
					for( int j=0 ; j<net.size() ; ++j){
				        average_neighbours_aux = average_neighbours_aux + (double) net[j]->get_n_neighbours();
				    }
				    average_neighbours_aux = average_neighbours_aux/n_nodes[i];
				}
				n_neighbours_aux.push_back(average_neighbours_aux);
			}

			// Random-geometric_networks
			else if( network_type == 2 ){
				squared_lattice_cyclic( n_nodes[i],n_nodes[i],(4*pow(fmin,4./n_nodes[i])-1)/3,pmin,pmin,net);
				n_neighbours_aux.push_back(4);
			}
			else {
				std::cout << "Invalid network type!" << std::endl;
				break;
			}
			
			for (int iSIM=0 ; iSIM<NSIMUL; ++iSIM){
				data << n_nodes[i] << "  ";
				dist trunc = dist(((4*fmin-1)/3),0,0,0,1000000,0);
				std::vector<std::shared_ptr<node> > netcopy;
				netcopy = copy_network(net);
				clock_t start, end;
				double elapsed_time;
				start = clock();
				shortest_path(netcopy,netcopy[iSIM],trunc,false);
				end = clock();
				elapsed_time = (double)(end - start) ;
				time_elapsed_aux.push_back(elapsed_time);
				int count = 0;
				int npathstotal = 0;
				for( int k=0 ; k<netcopy.size() ; ++k){
					int npathsaux = netcopy[k]->get_n_paths();
					npathstotal = npathstotal + netcopy[k]->get_n_paths();
					if( npathsaux != 0){
						n_paths_aux.push_back(npathsaux);
						++count;
					}
				}
				nodes_reached_aux.push_back((double)count/net.size());
				std::cout << "|" << std::flush;
				data << elapsed_time << "  " << (double)count/net.size() << "  " << (double)npathstotal/count << std::endl;
			}
		}
		time_elapsed_MEAN.push_back(average(time_elapsed_aux));
		time_elapsed_STD.push_back(standarddeviation(time_elapsed_aux));
		nodes_reached_MEAN.push_back(average(nodes_reached_aux));
		nodes_reached_STD.push_back(standarddeviation(nodes_reached_aux));
		n_neighbours_MEAN.push_back(average(n_neighbours_aux));
		n_neighbours_STD.push_back(standarddeviation(n_neighbours_aux));
		n_paths_MEAN.push_back(average(n_paths_aux));
		n_paths_STD.push_back(standarddeviation(n_paths_aux));
		std::cout << std::endl;
		std::cout << "Time elapsed: " << time_elapsed_MEAN[i] << " +/- " << time_elapsed_STD[i] <<	std::endl;
		std::cout << "Nodes reached: " << nodes_reached_MEAN[i] << " +/- " << nodes_reached_STD[i] <<	std::endl;
		std::cout << "Neighbours: " << n_neighbours_MEAN[i] << " +/- " << n_neighbours_STD[i] <<	std::endl;
		std::cout << "N paths per node: " << n_paths_MEAN[i] << " +/- " << n_paths_STD[i] <<	std::endl;
		outfile << n_nodes[i]*pow(n_nodes[i],network_type-1) << "   " << time_elapsed_MEAN[i] << "   "  << time_elapsed_STD[i] << "  " << nodes_reached_MEAN[i] << "  " << nodes_reached_STD[i] << "  " << n_neighbours_MEAN[i] << "  " << n_neighbours_STD[i] << "  " << n_paths_MEAN[i] << "  " << n_paths_STD[i] << std::endl;
	}
	
	outfile.close();
	data.close();
}


// Simulations of Multi-obective Shortest-Star algorithm 
// filename/filenamedata - containers for storing results
// n_nodes - size of the network in [n_min,n_max] with steps of size n_step
// fmin - minimum fidelity; pmin - minimum probability of success; nb - average degree (number of neighbours); nt - number of terminals
void simulation_star(char* filename, char* filenamedata, int network_type, int n_min, int n_max, int n_step, int nb, int nt, double fmin, double pmin){
	std::ofstream outfile (filename);
	std::ofstream data (filenamedata);

	int aux_n_nodes = (n_max-n_min)/n_step;
	std::vector<int> n_nodes;
	std::vector<double> time_elapsed_MEAN;
	std::vector<double> trees_created_MEAN;
	std::vector<double> time_elapsed_STD;
	std::vector<double> trees_created_STD;
	std::vector<double> n_neighbours_MEAN;
	std::vector<double> n_neighbours_STD;

	int NNETWORKS = 10;
	int NSIMUL = 10;

	dist trunc_star = dist(0.5,0,0,0,HUGE_VAL,0);
	
	for(int i=0 ; i<aux_n_nodes ; ++i)
		n_nodes.push_back(n_min+ i*n_step);

	for( int i=0 ; i<aux_n_nodes ; ++i){
		std::vector<double> time_elapsed_aux;
		std::vector<double> trees_created_aux;
		std::vector<double> n_neighbours_aux;
		std::cout << n_nodes[i] << std::endl;	
		for( int iNET=0 ; iNET < NNETWORKS ; ++iNET){
			std::vector<std::shared_ptr<node> > net;
			if( network_type == 1 ){
				double average_neighbours_aux = 0;
				while(average_neighbours_aux < 2){
					average_neighbours_aux = 0;
					// note minimum fidelity is calculated such that the network is functionally connected
					erdos_renyi_model_2(n_nodes[i],(double) nb / n_nodes[i],(4*pow(fmin,log((double)nb)/log(n_nodes[i]))-1)/3, pmin,pmin,net);
					for( int j=0 ; j<net.size() ; ++j){
				        average_neighbours_aux = average_neighbours_aux + (double) net[j]->get_n_neighbours();
				    }
				    average_neighbours_aux = average_neighbours_aux/n_nodes[i];
				}
				n_neighbours_aux.push_back(average_neighbours_aux);
			}
			else if( network_type == 2 )
				// note minimum fidelity is calculated such that the network is functionally connected
				squared_lattice_cyclic( n_nodes[i],n_nodes[i],(4*pow(fmin,4./n_nodes[i])-1)/3,pmin,pmin,net);
			else {
				std::cout << "Invalid network type!" << std::endl;
				break;
			}
			int breakflag = 0;
			for (int iSIM=0 ; iSIM<NSIMUL; ++iSIM){
				data << n_nodes[i] << "  ";
				std::vector<std::shared_ptr<node> > netcopy;
				netcopy = copy_network(net);
				dist trunc = dist(((4*fmin-1)/3),0,0,0,1000000,0);
				clock_t start, end;
				double elapsed_time;
				//Generating random terminals
				std::vector<int> terminal;
    			unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    			std::default_random_engine generator(seed);
    			std::uniform_int_distribution<int> uniform(0,n_nodes[i]*pow(n_nodes[i],network_type-1)-1);

    			int k = uniform(generator);
    			terminal.push_back(k);

    			while(terminal.size() != nt){
			        k = uniform(generator);
			        if( not_in(k,terminal))
			            terminal.push_back(k);
			    }

    			std::sort(terminal.begin(),terminal.end());

  				//Algorithm
				start = clock();
				std::vector<tree> trees;
				trees = shortest_star(netcopy,terminal,trunc,trunc);
				end = clock();
				elapsed_time = (double)(end - start);
				if(trees.size() != 0){
					time_elapsed_aux.push_back(elapsed_time);
					trees_created_aux.push_back(trees.size());
					std::cout << "|" << std::flush;
					data << elapsed_time << "  " << trees.size() << std::endl;
				}
				else if( breakflag > 20 ){
					//std::cout << "maximum number of attempts reached, breaking to next graph" << std::endl;
					break;
				}
				else {
					--iSIM;
					++breakflag;
				}

			}
		}
		time_elapsed_MEAN.push_back(average(time_elapsed_aux));
		time_elapsed_STD.push_back(standarddeviation(time_elapsed_aux));
		trees_created_MEAN.push_back(average(trees_created_aux));
		trees_created_STD.push_back(standarddeviation(trees_created_aux));
		n_neighbours_MEAN.push_back(average(n_neighbours_aux));
		n_neighbours_STD.push_back(standarddeviation(n_neighbours_aux));
		std::cout << std::endl;
		std::cout << "Time elapsed: " << time_elapsed_MEAN[i] << " +/- " << time_elapsed_STD[i] <<	std::endl;
		std::cout << "Trees Created: " << trees_created_MEAN[i] << " +/- " << trees_created_STD[i] <<	std::endl;
		outfile << n_nodes[i]*pow(n_nodes[i],network_type-1) << "   " << time_elapsed_MEAN[i] << "   "  << time_elapsed_STD[i] << "  " << n_neighbours_MEAN[i] << "  " << n_neighbours_STD[i] << "  " << trees_created_MEAN[i] << "  " << trees_created_STD[i] << std::endl;

	}

	std::cout << " ------------------------------------------------------------------------------------ " << std::endl;

	outfile.close();	
	data.close();
}
 