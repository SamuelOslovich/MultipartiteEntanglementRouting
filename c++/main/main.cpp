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
#include <math.h>


#include "node.h"
#include "dist.h"
#include "algorithm.h"
#include "network.h"
#include "simulation.h"


int main(void)
{
    
    std::vector<std::shared_ptr<node> > net;

    //SCALLING
    char  filename1[] = "../bin/ER_1000_scaling_3.txt";
    scaling(filename1,1,1000,1.,10.,0.5,0.9,0.9,3);
    char  filename2[] = "../bin/ER_5000_scaling_3.txt";
    scaling(filename2,1,5000,1.,10.,0.5,0.9,0.9,3);
    
    char  filename3[] = "../bin/SQ_10_scaling.txt";
    scaling(filename3,2,10,1.,10.,0.5,0.9,0.9); 
    char  filename3[] = "../bin/SQ_50_scaling.txt";
    scaling(filename3,2,50,3.,10.,0.5,0.9,0.9);
    
    
    //PATHS
    char  filenameERpath[] = "../bin/ER_path_3.txt";
    char  filenameERpath_d[] = "../bin/ER_path_3data.txt";
    simulation_path(filenameERpath,filenameERpath_d,1,100,10100,100,3,0.9,0.5);

    char  filenameGEOpath[] = "../bin/Geometric_path .txt";
    char  filenameGEOpathd[] = "../bin/Geometric_path_data.txt";
    simulation_path(filenameGEOpath,filenameGEOpathd,2,100,1250,50,3,0.9,0.5);
    
    //TREES
    //3 terminais
    char  filenameERtree3[] = "../bin/ER_tree_3.txt";
    char  filenameERtree3_d[] = "../bin/ER_tree_3data.txt";
    simulation_star(filenameERtree3,filenameERtree3_d,1,100,5100,100,3,3,0.9,0.5);

    char  filenameERtree4[] = "../bin/ER_tree_4.txt";
    char  filenameERtree4_d[] = "../bin/ER_tree_4data.txt";
    simulation_star(filenameERtree4,filenameERtree4_d,1,100,5100,100,3,4,0.9,0.5);

    //TREES
    //4 terminais
    char  filenameGEOtree3[] = "../bin/Geometric_tree_3.txt";
    char  filenameGEOtree3d[] = "../bin/Geometric_tree_3_data.txt";
    simulation_star(filenameGEOtree3,filenameGEOtree3d,1,100,1250,50,3,3,0.9,0.5);

    char  filenameGEOtree4[] = "../bin/Geometric_tree_4.txt";
    char  filenameGEOtree4d[] = "../bin/Geometric_tree_4data.txt";
    simulation_star(filenameGEOtree4,filenameGEOtree4d,1,100,1250,50,3,4,0.9,0.5);


    /*
    //USING ALGORITHM TO FIND THE STARS ------------------------------------------

    N=1000

    std::cout << "Constructing network ... building ...." << std::endl;
    double ftrunc = 0.5;
    double gammatrunc = (4*ftrunc-1)/3;
    double pmin = 0.5;

    clock_t start, end;

    dist trunc = dist(gammatrunc,0,0,0,HUGE_VAL,0);
    dist truncstar = dist(ftrunc,0,0,0,HUGE_VAL,0);
    int NSIMUL = 100;


    for( int j= 0 ; j<NSIMUL ; ++j){
        double gammamin = pow(gammatrunc,log(N)/sqrt(N));
        create_grid_from_file("../networks/Internet/", std::vector<std::shared_ptr<node> >& nodes,double fmin, double p_gen_min, double p_swap_min, int n_nodes);
        int k=0;
        while( k< 10 ){
            auto auxnet = copy_network(net);
            
            std::vector<int> terminal;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator(seed);
            std::uniform_int_distribution<int> uniform(0,N-1);

            int k1 = uniform(generator);
            terminal.push_back(k1);

            while(terminal.size() != 3){
                k1 = uniform(generator);
                if( not_in(k1,terminal))
                    terminal.push_back(k1);
            }

            std::sort(terminal.begin(),terminal.end());

            for( int i=0 ; i<terminal.size() ; ++i)
                std::cout << terminal[i] << std::endl;

            std::vector<tree> sol;
            start = clock();
            sol = shortest_star(net,terminal,trunc,truncstar,datapath);
            end = clock();

            double elapsed_time = (double)(end - start);

            datatree << N << " " << elapsed_time << " " << sol.size() << std::endl;
            std::cout << N << " " << elapsed_time << " " << sol.size() << std::endl;
        
            ++k;

        }
    }

    datapath.close();
    datatree.close();

    //------------------------------------------------------------------------------
    //USING ALGORITHM TO FIND STARS
    /*
    std::cout << "Constructing network ... building ...." << std::endl;
    int N=1000;
    int Nterminais = 3;
    double ftrunc = 0.9;
    double gammatrunc = (4*ftrunc-1)/3;
    double gammamin = pow((4*ftrunc-1)/3,log(3.)/log(N));
    double pmin = 0.9;

    clock_t start, end;
    start = clock();
    erdos_renyi_model_2(N, 3./ N, gammamin, pmin, pmin, net);
    end = clock();

    std::cout << "Elapsed time (clock cycles): " << (double) end-start << std::endl;

    dist trunc1 = dist(gammatrunc,0,0,0,HUGE_VAL,0);
    dist trunc2 = dist(0.5,0,0,0,HUGE_VAL,0);

    std::vector<int> terminal;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> uniform(0,N-1);

    int k = uniform(generator);
    terminal.push_back(k);

    while(terminal.size() != Nterminais){
        k = uniform(generator);
        if( not_in(k,terminal))
            terminal.push_back(k);
    }

    std::sort(terminal.begin(),terminal.end());

    for( int i=0 ; i<terminal.size() ; ++i)
        std::cout << terminal[i] << std::endl;

    std::vector<tree> sol;
    start = clock();
    sol = shortest_star(net,terminal,trunc1,trunc2);
    end = clock();


    for( int i=0 ; i<sol.size() ; ++i){
        sol[i].print();
    }

    std::cout << "Number of trees solution: " << sol.size() << std::endl;
    std::cout << "Elapsed time (clock cycles): " << (double) end-start << std::endl;
    */
    //--------------------------------------------------------------------------------
    //USING ALGORITHM TO FIND THE PATHS SQUARE NETWORK ------------------------------------------
    /*
    std::cout << "Constructing network ... building ...." << std::endl;
    int N=22000;
    double fmin = 0.9;
    double pmin = 0.9;

    clock_t start, end;
    start = clock();
    //squared_lattice_cyclic(N , N, (4*pow(fmin,4./N)-1)/3, pmin, pmin, net);
    erdos_renyi_model_2(N,3./N,(4*pow(fmin,1./log(N))-1)/3, pmin, pmin, net);
    end = clock();

    std::cout << "Elapsed time (clock cycles): " << (double) end-start << std::endl;

    dist trunc = dist((4*fmin-1)/3,0,0,0,1000000,0);

    start = clock();
    shortest_path(net,net[0],trunc,true);
    end = clock();

    int count =0;
    int npathsall = 0;

    for( int i=0 ; i<net.size() ; ++i){
        if( net[i]->get_n_paths() != 0){
            ++count;
            npathsall += net[i]->get_n_paths();
        }
        //net[i]->print(false);
    }
    

    

    std::cout << "Number of nodes reached: " << count << std::endl;
    std::cout << "Percentage of nodes reached: " << (double) count/(N) << std::endl;
    std::cout << "Number of paths per node: " << (double) npathsall / (double) (count) << std::endl;

    std::cout << "Elapsed time (clock cycles): " << (double) end-start << std::endl;
    */
    //------------------------------------------------------------------------------
    //USING ALGORITHM TO FIND STARS IN SQUARE NETWORK ------------------------------------------
    /*
    std::cout << "Constructing network ... building ...." << std::endl;
    int N=10;
    int Nterminais = 3;
    double fmin = 0.9;
    double pmin = 0.9;

    clock_t start, end;
    start = clock();
    squared_lattice_cyclic(N , N, (4*pow(fmin,1./N)-1)/3, pow(pmin,2.5/N), pow(pmin,2./N), net);
    end = clock();

    std::cout << "Elapsed time (clock cycles): " << (double) end-start << std::endl;

    dist trunc = dist(N,4*(fmin-1)/3,0,0,pmin);

    std::vector<int> terminal;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> uniform(0,N-1);

    int k = uniform(generator);
    terminal.push_back(k);

    while(terminal.size() != Nterminais){
        k = uniform(generator);
        if( not_in(k,terminal))
            terminal.push_back(k);
    }

    std::sort(terminal.begin(),terminal.end());

    for( int i=0 ; i<terminal.size() ; ++i)
        std::cout << terminal[i] << std::endl;

    std::vector<tree> sol;
    start = clock();
    sol = shortest_star(net,terminal,trunc);
    end = clock();

    for( int i=0 ; i<sol.size() ; ++i)
        sol[i].print();

    std::cout << "Number of trees solution: " << sol.size() << std::endl;

    std::cout << "Elapsed time (clock cycles): " << (double) end-start << std::endl;
    */
    //------------------------------------------------------------------------------------------------
    /*
    //USING ALGORITHM TO FIND THE PATHS ------------------------------------------

    std::cout << "Constructing network ... building ...." << std::endl;

    double fmin = 0.5;
    int N =10;

    clock_t start, end;
    start = clock();
    create_grid_from_file("grid10.txt", net, N);
    end = clock();

    std::cout << "Elapsed time (clock cycles): " << (double) end-start << std::endl;

    dist trunc = dist(4*(fmin-1)/3,0,0,0);

    start = clock();
    shortest_path(net,net[0],trunc,true);
    end = clock();

    std::cout << "oal" << std::endl;
    int count =0;
    int npathsall = 0;

    for( int i=0 ; i<net.size() ; ++i){
        if( net[i]->get_n_paths() != 0){
            ++count;
            npathsall += net[i]->get_n_paths();
        }
        net[i]->print(false);
        //netcopy[i]->print(); 
    }

    std::cout << "Number of nodes reached: " << count << std::endl;
    std::cout << "Percentage of nodes reached: " << (double) count/(N) << std::endl;
    std::cout << "Number of paths per node: " << (double) npathsall / (double) (count) << std::endl;

    std::cout << "Elapsed time (clock cycles): " << (double) end-start << std::endl;
    */
    //USING ALGORITHM TO FIND STARS -----------------------------------------------------------------
    /*
    std::cout << "Constructing network ... building ...." << std::endl;
    
    double fmin = 0.5;
    int N =10;

    clock_t start, end;

    start = clock();
    create_grid_from_file("grid10.txt", net, N);
    end = clock();

    std::cout << "Elapsed time (clock cycles): " << (double) end-start << std::endl;

    dist trunc = dist(4*(fmin-1)/3,0,0,0);
    std::vector<int> terminal;

    terminal.push_back(0);
    terminal.push_back(4);
    terminal.push_back(9);

    std::sort(terminal.begin(),terminal.end());

    std::vector<tree> sol;
    start = clock();
    sol = shortest_star(net,terminal,trunc);
    end = clock();

    for( int i=0 ; i<sol.size() ; ++i)
        sol[i].print(false);

    std::cout << "Number of trees solution: " << sol.size() << std::endl;

    std::cout << "Elapsed time (clock cycles): " << (double) end-start << std::endl;
    */
    return 0;
}
