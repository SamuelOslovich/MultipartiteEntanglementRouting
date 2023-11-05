import sys
import random
import math
import numpy as np
import networkx as nx
from graphs import *
from weights import *
from algorithms import *
from scipy.stats import powerlaw


def get_power_law(alpha,min_value,max_value):
	x = powerlaw.rvs(alpha, size=1)
	return max_value - x[0]*(max_value-min_value)

## Sam Comment
# Calculates the binary entropy of p
def H2(p):
	if p == 0: 
		return 0
	else:
	#print(- p*math.log2(p) - (1-p)*math.log2(1-p))
		return - p*math.log2(p) - (1-p)*math.log2(1-p)

def k(p):
	return 1-H2(3*p/4)

## Sam Comment
# filename: a file that contains the edges for the graph
# fid_min: the minimum fidelity
# prob_min: the min successful prob of entang
# n_nodes: number of nodes
# dist_t: distribution to use for the probability of successful entang
def create_grid_from_file(filename,fid_min,prob_min,n_nodes,dist_t="uniform"):
	t_mem_min = 100000
	t_mem_max = 1000000

	time_min = 1
	time_max = 100

	fid_min = fid_min
	fid_max = 1.0

	prob_min = prob_min
	prob_max = 1.0
	dist_type = random.uniform



	file = open("../networks/"+filename+".txt","r")
	edges = file.readlines()

	## Sam Comment
	# G1 uses probability and fidelity for each link
	# G2 is just like from the graph paper (no prob, no fidelity, everything is perfect)
	# G3 just uses bounds-based: "where we perform the routing based on the maximum bound
	#                             on what a protocol could achieve on distributing
	#                             GHZ states. To each edge, one attributes a capacity that 
	#                             translates the loss in the channel, and the probabilistic 
	# 							  nature of entanglement generation. Then the capacity of 
	#                             distribution is bounded by the minimum capacity along the
	#                             links that make a given path."
	# nx.Graph is a default graph from a library
	# Graph is a custom implementation in Graphs.py that implements a graph with "quantum channels" as the links

	G1 = Graph()
	G2 = nx.Graph()
	G3 = nx.Graph()
    
	## Sam Comment
	# Create all of the nodes
	# Only G1 cares about the probability of successful entag swapping or the quantum memory
	for i in range(0,n_nodes):
		t_mem = dist_type(t_mem_min,t_mem_max)
		prob_swap = dist_type(prob_min,prob_max)
		name = i
		G1.add_node(Node(name=name,neighbours={},prob_swap=prob_swap,t_mem=t_mem))
		G2.add_node(name)
		G3.add_node(name)
		#print(G1.nodes[name])

	## Sam Comment
	# Add all of the edges (read from the passed in file)
	# Add (and generate) the details for the links for G1. (fidelity, probability of entang generation, comm time)
	# For G3 add the max capacity (loss of the channel) as the weight
	for edge in edges: 
		e = edge.strip().split(' ')
		s = int(e[0])
		t = int(e[1])

		fid = None
		if dist_t == "uniform":
			fid = dist_type(fid_min,fid_max)
		elif dist_t == "power_law":
			fid = get_power_law(2,fid_min,fid_max)
		time = dist_type(time_min,time_max)
		prob_gen = dist_type(prob_min,prob_max)
		linka = Link(source=G1.nodes[s], target=G1.nodes[t], fid=fid, prob_gen=prob_gen, time=time)
		#print(linka)
		G1.add_link(linka)
		#G1.print()
		G2.add_edge(s,t)
		G3.add_edge(s,t,weight=Capacity(cap=prob_gen*k(1-fid),length=1))

	return G1, G2, G3

## Sam Comment
# New function I created to create a graph
# nodeFilename: a space separated file that contains all of the nodes (in order), 
#               the probability of successful entanglement swapping for the node,
#               and the memory decoherence time for the node  
# edgeFilename: a file that contains the edges, the probability of generating 
# 				entanglement for each edge, the fidelity for each edge, and the 
# 				communication time for each edge
# n_nodes: the number of nodes
def create_grid_from_file_complex(nodeFilename, edgeFilename, n_nodes):
	nodeFile = open("../networks/" + nodeFilename + ".txt", "r")
	nodes = nodeFile.readlines()

	edgeFile = open("../networks/" + edgeFilename + ".txt", "r")
	edges = edgeFile.readlines()

	G1 = Graph()
	G2 = nx.Graph()
	G3 = nx.Graph()

	for node in nodes:
		nodeData = node.strip().split(' ')
		name = int(nodeData[0])
		prob_swap = float(nodeData[1])
		t_mem = float(nodeData[2])

		G1.add_node(Node(name=name,neighbours={},prob_swap=prob_swap,t_mem=t_mem))
		G2.add_node(name)
		G3.add_node(name)

	for edge in edges:
		edgeData = edge.strip().split(' ')
		s = int(edgeData[0])
		t = int(edgeData[1])
		prob_gen = float(edgeData[2])
		fid = float(edgeData[3])
		time = float(edgeData[4])

		linka = Link(source=G1.nodes[s], target=G1.nodes[t], fid=fid, prob_gen=prob_gen, time=time)
		G1.add_link(linka)
		G2.add_edge(s,t)
		G3.add_edge(s,t,weight=Capacity(cap=prob_gen*k(1-fid),length=1))

	return G1, G2, G3

## Sam Comment
# New function I created to create a uniform graph where all params are the same for each node/link
# filename: a file that contains the edges for the graph
# fid: the fidelity of each link
# prob: the probability of successful entanglement/(entang swapping)
# t_com: the communication time between nodes
# t_mem: the decoherence time of each nodes quantum memory
# n_nodes: number of nodes
def create_static_grid_from_file(filename, fid, prob, t_com, t_mem, n_nodes):
	file = open("../networks/"+filename+".txt","r")
	edges = file.readlines()

	## Sam Comment
	# G1 uses probability and fidelity for each link
	# G2 is just like from the graph paper (no prob, no fidelity, everything is perfect)
	# G3 just uses bounds-based: "where we perform the routing based on the maximum bound
	#                             on what a protocol could achieve on distributing
	#                             GHZ states. To each edge, one attributes a capacity that 
	#                             translates the loss in the channel, and the probabilistic 
	# 							  nature of entanglement generation. Then the capacity of 
	#                             distribution is bounded by the minimum capacity along the
	#                             links that make a given path."
	# nx.Graph is a default graph from a library
	# Graph is a custom implementation in Graphs.py that implements a graph with "quantum channels" as the links

	G1 = Graph()
	G2 = nx.Graph()
	G3 = nx.Graph()
    
	## Sam Comment
	# Create all of the nodes
	# Only G1 cares about the probability of successful entag swapping or the quantum memory
	for i in range(0,n_nodes):
		name = i
		G1.add_node(Node(name=name,neighbours={},prob_swap=prob,t_mem=t_mem))
		G2.add_node(name)
		G3.add_node(name)
		#print(G1.nodes[name])

	## Sam Comment
	# Add all of the edges (read from the passed in file)
	# Add (and generate) the details for the links for G1. (fidelity, probability of entang generation, comm time)
	# For G3 add the max capacity (loss of the channel) as the weight
	for edge in edges: 
		e = edge.strip().split(' ')
		s = int(e[0])
		t = int(e[1])

		linka = Link(source=G1.nodes[s], target=G1.nodes[t], fid=fid, prob_gen=prob, time=t_com)
		#print(linka)
		G1.add_link(linka)
		#G1.print()
		G2.add_edge(s,t)
		G3.add_edge(s,t,weight=Capacity(cap=prob*k(1-fid),length=1))

	return G1, G2, G3
		

## Sam Comment
# Same as create_grid_from_file, except all of the probabilities are 1, and decoherence in quantum memory doesn't happen
def create_grid_from_file_simple(filename,probs,n_nodes):
	t_mem = math.inf
	prob_swap = 1.
	time = 1
	fid = 1.

	file = open("../networks/"+filename+".txt","r")
	edges = file.readlines()

	G1 = Graph()
	G2 = nx.Graph()
	G3 = nx.Graph()
    
	for i in range(0,n_nodes):
		
		name = i
		G1.add_node(Node(name=name,neighbours={},prob_swap=prob_swap,t_mem=t_mem))
		G2.add_node(name)
		G3.add_node(name)
		#print(G1.nodes[name])

	for edge in edges: 
		e = edge.strip().split(' ')
		s = int(e[0])
		t = int(e[1])

		prob_gen = random.sample(probs,1)[0]
		linka = Link(source=G1.nodes[s], target=G1.nodes[t], fid=fid, prob_gen=prob_gen, time=time)
		G1.add_link(linka)
		G2.add_edge(s,t)
		G3.add_edge(s,t,weight=Capacity(cap=prob_gen*k(1-fid),length=1))

	return G1, G2, G3

def compare_stars(graph,star1,star2,star3):
	star_chosen_rate = None
	star_chosen_fid = None
	
	star_rate = 0
	star_fid = 0
	#print(star1[0])
	#print(star2)
	#print(star3)
	#Finding optimal rate and optimal fidelity
	for star in star1:
		if star.weight.rate > star_rate:
			star_rate = star.weight.rate
			star_chosen_rate = star

	for star in star1:
		if star.weight.fid > star_fid:
			star_fid = star.weight.fid
			star_chosen_fid = star

	point = {"n_nodes" : len(graph.nodes),
			 "rate_shortest" : 0,
			 "fid_shortest" : 0,
			 "rate_bound" : 0,
			 "fid_bound" : 0,
			 "n_links" : star2["dist"],
			 "n_links_bound" : sum([len(path) for path in star3["star"]])-3,
			 "valid_shortest" : True,
			 "valid_bound" : True}

	# NO PATH CONNECTING NODES
	for path in star2:
		if len(path) == 0:
			point["valid_shortest"] = False
			point["valid_bound"] = False

			return point

	star_rec_shortest = [reconstruct_path(graph,path) for path in star2["star"]]
	star_rec_bound = [reconstruct_path(graph,path) for path in star3["star"]]

	
	# FOR A FAIR COMPARISON - remember we are excluding paths with fid < (4*0.5**(1/3)-1)/3
	for path in star_rec_shortest:
		if path.weight.fid < (4*0.5**(1/3)-1)/3:
			point["valid_shortest"] = False

	for path in star_rec_bound:
		if path.weight.fid < (4*0.5**(1/3)-1)/3:
			point["valid_bound"] = False

	if star_rate == 0 or star_fid == 0:
		point["valid_shortest"] = False
		point["valid_bound"] = False

	if point["valid_shortest"]:
		star_compare_shortest = reconstruct_star(graph,list_of_paths=star_rec_shortest)
		point["rate_shortest"] = star_compare_shortest.weight.rate / star_rate
		point["fid_shortest"] = star_compare_shortest.weight.fid / star_fid

	if point["valid_bound"]:
		star_compare_bound = reconstruct_star(graph,list_of_paths=star_rec_bound)
		point["rate_bound"] = star_compare_bound.weight.rate / star_rate
		point["fid_bound"] = star_compare_bound.weight.fid / star_fid

	
	if point["rate_shortest"] > 1.0 or point["fid_shortest"]>1.0 or point["rate_bound"] > 1.0 or point["fid_bound"]>1.0:
		print("SOMETHING FUNNY IS GOING ON")
		print("BEST RATE STAR")
		print(star_chosen_rate)
		print("BEST FID STAR")
		print(star_chosen_fid)
		print("BEST DISTANCE STAR")
		print(star_compare)
		print(star_compare.get_center_node().print_all_paths())


	return point



