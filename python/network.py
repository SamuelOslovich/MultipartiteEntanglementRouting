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

def H2(p):
	if p == 0: 
		return 0
	else:
	#print(- p*math.log2(p) - (1-p)*math.log2(1-p))
		return - p*math.log2(p) - (1-p)*math.log2(1-p)

def k(p):
	return 1-H2(3*p/4)

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

	G1 = Graph()
	G2 = nx.Graph()
	G3 = nx.Graph()
    
	for i in range(0,n_nodes):
		t_mem = dist_type(t_mem_min,t_mem_max)
		prob_swap = dist_type(prob_min,prob_max)
		name = i
		G1.add_node(Node(name=name,neighbours={},prob_swap=prob_swap,t_mem=t_mem))
		G2.add_node(name)
		G3.add_node(name)
		#print(G1.nodes[name])

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



