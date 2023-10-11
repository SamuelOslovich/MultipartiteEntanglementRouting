import sys
import random
import math
import numpy as np
import networkx as nx
from graphs import *
from weights import *
from algorithms import *
from network import *
from scipy import stats
import pickle
import time
timestr = time.strftime("%Y%m%d-%H")
print(timestr)

## Sam Comment
# Note that anytime fidelity is mentioned, it is gamma (from the paper)
# Is the rate calculation incorrect? Shouldn't it be prob / (2*max time)

# sys.argv[1] - type of simulation
# 	. test - debugging on small graphs
#	. convergence - reproduce convergence results on internet-like networks
#	. ER - reproduce optimality results on Erdos-Renyi networks
# 	. internet - reproduce optimality results on internet-like networks

# ----------------------------- TESTING -----------------------------

if sys.argv[1] == "test":
	N_nodes = 100
	N_terminals = 3
	g1,g2,g3 = create_grid_from_file("Dummy/0",0.9,0.5,N_nodes,dist_t="uniform")

	g1.print()
	print(g1.average_degree())

	terminal = random.sample(range(N_nodes), N_terminals)
	
	## Sam comment
	# used to test that the optimal paths algorithm is working
	# optimal paths finds all of the paths in g1 between a source (g1[terminal[0]]) and all other nodes.
	# paths that don't satisfy the weight_trunc are not included
	optimal_paths(g1,g1[terminal[0]],weight_trunc=Weight(fid=0.72, prob=0, time=math.inf, sigma=0))
	g1.print_paths(g1[terminal[0]])
	g1.clear()
	g1.print_paths(g1[terminal[0]])
	
	## Sam comment
	# verifies that shortest_path is working
	sol = nx.shortest_path(g2,source=terminal[0])
	print(sol)
	
	## Sam comment
	# verifies that shortest_paths_general is working
	solpaths = shortest_paths_general(graph=g3,source=terminal[0],special_signature=Capacity(0,math.inf),neutral_signature=Capacity())
	print(solpaths)

	## Sam comment
	# checks that no paths are repeated
	for key,value in solpaths.items():
		seen = []
		for number in value[1]:
			if number in seen:
				print("Number repeated!")
			else:
				seen.append(number)

	stars = optimal_star(g1,terminal, qkd_flag=True)
	star2 = shortest_star(g2,terminal)
	star3 = shortest_star_general(graph=g3,terminal=terminal,special_signature=Capacity(0,math.inf),neutral_signature=Capacity())
	print(star3)

	for star in stars:
		print(star)

	if len(stars) > 0:
		print(compare_stars(g1,stars,star2,star3))

# ----------------------------- SIMULATING - CONVERGENCE OF METHODS ---------------------------

if sys.argv[1] == "convergence":

	N_nodes = 1000
	N_tries = 100
	N_simul = 40
	p_min = 0.5
	delta_p = (1.-p_min)/N_simul

	N_terminals = 3
	values_p = [p_min + i*delta_p for i in range(N_simul+1)]
	values_p.reverse()
	print(values_p)

	convergence_points = []
	filename = 'Convergence_points_' + timestr 
	fileConv=open( filename + '.p', 'wb')

	for i in range(len(values_p)-1,-1,-1):
		for j in range(N_tries):	
			g1,g2,g3 = create_grid_from_file_simple("Internet/"+str(j),values_p[0:i+1],N_nodes)
			print(values_p[0:i+1])
			print(" ------ TRY " + str(j) + " - " + str(values_p[i]) + " -----------")

			terminal = random.sample(range(N_nodes), N_terminals)

			stars = optimal_star(g1,terminal)
			star2 = shortest_star(g2,terminal)
			star3 = shortest_star_general(g3,terminal,special_signature=Capacity(0,math.inf),neutral_signature=Capacity())

			dict_aux = compare_stars(g1,stars,star2,star3)
			dict_aux["graph"] = "Internet/"+str(j)
			dict_aux["dist_type"] = "uniform"
			dict_aux["p"] = values_p[i]

			print(dict_aux)

			pickle.dump(dict_aux.copy(),fileConv)
			convergence_points.append(dict_aux)

			g1.clear()

			del g1, g2, g3

	np.save(filename + ".npy",convergence_points)
	fileConv.close()


# ----------------------------- SIMULATING - UNIFORM DISTRIBUTION -----------------------------

if sys.argv[1] == "ER":

	N_nodes = 1000
	N_simul = 200
	N_tries = 20
	N_terminals = 3

	ER_points = []
	filename = 'ER_points_' + timestr 
	fileER=open( filename + '.p', 'wb')


	for i in range(0,N_simul):	
		g1,g2,g3 = create_grid_from_file("ER/"+str(i),0.9,0.5,N_nodes,dist_t="uniform")
		
		for j in range(0,N_tries):
			print(" ------ TRY " + str(i) + " - " + str(j) + " -----------")
			terminal = random.sample(range(N_nodes), N_terminals)

			stars = optimal_star(g1,terminal)
			star2 = shortest_star(g2,terminal)
			star3 = shortest_star_general(g3,terminal,special_signature=Capacity(0,math.inf),neutral_signature=Capacity())

			dict_aux = compare_stars(g1,stars,star2,star3)
			dict_aux["graph"] = "ER/"+str(i)
			dict_aux["dist_type"] = "uniform"

			print(dict_aux)

			pickle.dump(dict_aux.copy(),fileER)
			ER_points.append(dict_aux)

		g1.clear()

		del g1, g2, g3

	np.save(filename + ".npy",ER_points)
	fileER.close()

if sys.argv[1] == "internet":
	
	N_nodes = 1000
	N_simul = 200
	N_tries = 20
	N_terminals = 3

	Internet_points = []
	filename = 'Internet_points_' + timestr
	fileInternet=open(filename + '.p', 'wb')

	for i in range(0,N_simul):	
		g1,g2,g3 = create_grid_from_file("Internet/"+str(i),0.9,0.5,N_nodes,dist_t="uniform")
			
		for j in range(0,N_tries):
			print(" ------ TRY " + str(i) + " - " + str(j) + " -----------")

			terminal = random.sample(range(N_nodes), N_terminals)

			stars = optimal_star(g1,terminal)
			star2 = shortest_star(g2,terminal)
			star3 = shortest_star_general(g3,terminal,special_signature=Capacity(0,math.inf),neutral_signature=Capacity())

			dict_aux = compare_stars(g1,stars,star2,star3)
			dict_aux["graph"] = "Internet/"+str(i)
			dict_aux["dist_type"] = "uniform"

			print(dict_aux)

			pickle.dump(dict_aux.copy(),fileInternet)
			Internet_points.append(dict_aux)

		g1.clear()

		del g1, g2, g3

	np.save(filename + ".npy",Internet_points)
	fileInternet.close()








