import sys
import random as rd
import math
import numpy as np
import networkx as nx
from weights import *
from graphs import *
from heap import *
from queue import PriorityQueue

## --------------------------------------- ALGORITHMS --------------------------------------- ##


## ------------------------------------------ PATHS ----------------------------------------- ##


# SHORTEST_PATHS_M - shortest path algorithm with metric passed as an argument of algorithm

# NOTE: Using regular ordering >= where a maximisation is in place. Necessary to change for the case
# of a general ordering, by creating another class and defining the ordering relations in that class
# works with NETWORKX definition of graphs!
def shortest_paths_m(graph,source,metric=sum,special_signature=0):
    sol = {v:[special_signature,[]] for v in graph.nodes}
    
    graph_visited = {v:False for v in graph.nodes}
    sol[source][0] = math.inf
    sol[source][1] = [source]


    pq = PriorityQueue()
    pq.put((-sol[source][0], source))

    while not pq.empty():
        (dist, current_vertex) = pq.get()
        graph_visited[current_vertex] = True

        #print(current_vertex)
        #print(graph[current_vertex])
        for neighbor in graph[current_vertex]:
            distance = graph[current_vertex][neighbor]["weight"]
            if not graph_visited[neighbor]:
                old_cost = sol[neighbor][0]
                new_cost = metric(sol[current_vertex][0],distance)
                if new_cost > old_cost:
                	pq.put((-new_cost, neighbor))
                	sol[neighbor][0] = new_cost              
                	sol[neighbor][1] = sol[current_vertex][1].copy() + [neighbor]

                if new_cost == old_cost:
                	if len(sol[current_vertex][1]) < len(sol[neighbor][1]):
                		pq.put((-new_cost, neighbor))
                		sol[neighbor][0] = new_cost              
                		sol[neighbor][1] = sol[current_vertex][1].copy() + [neighbor]
                  
    return sol



# SHORTEST_PATHS - shortest path algorithm with user defined weights from graph construction

# Graph with weights defined by an Object - ALGEBRA
# Object must have defined at least a __lt__ method and a __add__ method (compare and concatenate paths respectively)
# Special_signature: impossible path signature (e.g. capacity = 0 , length = infty , ...)
# Neutral_signature: object + neutral = object (group neutral element) (e.g. capacity = infty since min(infty,a) = a, forall a finite)
def shortest_paths_general(graph,source,special_signature,neutral_signature):
    sol = {v:[special_signature,[]] for v in graph.nodes}
    #print(sol)
    
    graph_visited = {v:False for v in graph.nodes}
    sol[source][0] = neutral_signature
    sol[source][1] = [source]


    pq = PriorityQueue()
    pq.put((sol[source][0], source))

    while not pq.empty():
        (dist, current_vertex) = pq.get()
        graph_visited[current_vertex] = True

        for neighbor in graph[current_vertex]:
            distance = graph[current_vertex][neighbor]["weight"]
            if not graph_visited[neighbor]:
                old_cost = sol[neighbor][0]
                new_cost = sol[current_vertex][0] + distance
                if new_cost < old_cost:
                	pq.put((new_cost, neighbor))
                	sol[neighbor][0] = new_cost              
                	sol[neighbor][1] = sol[current_vertex][1].copy() + [neighbor]
                  
    return sol


# OPTIMAL_PATHS - multi-objective shortest-path algorithm for the setup presented over the paper
# each path has a fidelity (fid), a probability of success (prob), a time of comm. (time) and a rate of decoherence (sigma).
# the way weights are added can be seen from the weights.py
def optimal_paths(graph,source,weight_trunc=Weight(fid=0.3333, prob=0, time=math.inf, sigma=0)):
	heap = FibonacciHeap(source)
	heap.insert_node(source)
	for node in graph.nodes.values():
		node.paths[source] = []
	source.initialize()


	while heap.count != 0:
		#retrieve top of the heap
		node_iter = heap.extract_min()
		node_iter.visited[source] = True

		#go to each neighbour and update their list of paths and add them to the heap
		for n, link in node_iter.neighbours.items():
			paths_add = [path + link for path in node_iter.paths[source] if (path + link).check_possible(weight_trunc)]

			if len(paths_add) == 0:
				# print("No paths!")
				continue	

			if len(n.paths[source]) == 0:
				n.update_paths(paths_add,source=source)

				heap.insert_node(n)

			elif len(n.paths[source]) > 0:
				flag_change = n.update_paths(paths_add,source=source)
				if flag_change:
					if n.visited.get(source):
						heap.insert_node(n)
						n.visited[source] = False
					else:
						heap.consolidate()

	print("OPTIMAL PATHS ALGORITHM FINISHED RUNNING!####################################################")


#Finding all possible optimal paths from each node
def all_optimal_paths(graph):
	for node in graph.nodes:
		optimal_paths(graph,source=node)


## ------------------------------------------ TREES ----------------------------------------- ##

# SHORTEST_STAR - most simple shortest-star algorithm with each link having weight-1
def shortest_star(graph,terminal):
	shortest_paths = {}
	for node in terminal:
		shortest_paths[node] = nx.shortest_path(graph, source=node)
		#print(shortest_paths[node])
	
	solution = {"c":None, "max": math.inf, "dist":math.inf, "star":[]}
	for central_node in graph.nodes:
		distance = 0
		maximum = 0
		
		for node in terminal:
			if not shortest_paths[node].get(central_node):
				distance = math.inf
				maximum = math.inf
				break
			#print("testing " + str(central_node) + " " + str(node))
			maximum = max(maximum,len(shortest_paths[node][central_node])-1)
			distance += len(shortest_paths[node][central_node])-1
		
		if maximum < solution["max"]:
			solution = {"c":central_node,
						"max":maximum, 
						"dist":distance, 
						"star":[shortest_paths[node][central_node] for node in terminal]}
		elif maximum == solution["max"] and distance < solution["dist"]:
			solution = {"c":central_node,
						"max":maximum, 
						"dist":distance, 
						"star":[shortest_paths[node][central_node] for node in terminal]}
				
	#print(shortest_paths)
	return solution

# SHORTEST_STAR_M - shortest-star algorithm with input metric
def shortest_star_m(graph,terminal,metric=sum):
	shortest_paths = {}
	for node in terminal:
		shortest_paths[node] = shortest_paths_m(graph, source=node,metric=metric)
		print(shortest_paths[node])
	
	solution = {"c":None, "dist":0, "star":[]}
	for central_node in graph.nodes:
		distance = math.inf
		for node in terminal:
			if not shortest_paths[node].get(central_node):
				distance = math.inf
				break
			#print("testing " + str(central_node) + " " + str(node))
			distance = metric(distance,shortest_paths[node][central_node][0])
		
		if distance < solution["dist"]:
			solution = {"c":central_node, 
						"dist":distance, 
						"star":[shortest_paths[node][central_node][1] for node in terminal]}
	#print(shortest_paths)
	return solution


# SHORTEST_STAR_GENERAL - shortest-star algorithm with additional defined weights, in this case using capacities
def shortest_star_general(graph,terminal,special_signature,neutral_signature):
	shortest_paths = {}
	for node in terminal:
		shortest_paths[node] = shortest_paths_general(graph, source=node,special_signature=special_signature,neutral_signature=neutral_signature)
		#print(shortest_paths[node])
	
	solution = {"c":None, "dist":TreeCapacity(special_signature), "star":[]}
	for central_node in graph.nodes:
		distance = TreeCapacity(neutral_signature)
		for node in terminal:
			if not shortest_paths[node].get(central_node):
				distance = 0
				break
			#print("testing " + str(central_node) + " " + str(node))
			distance = distance + TreeCapacity(shortest_paths[node][central_node][0])
		
		if distance < solution["dist"]:
			solution = {"c":central_node, 
						"dist":distance, 
						"star":[shortest_paths[node][central_node][1] for node in terminal]}
	#print(shortest_paths)
	return solution

# OPTIMAL_STAR - multi-objective optimal-star algorithm 
def optimal_star(graph,terminal,fid_trunc=0.5):

	weight_trunc = Weight(fid=(4*fid_trunc**(1/len(terminal))-1)/3, prob=0, time=math.inf, sigma=0)
	tree_trunc = WeightTree([WeightPath(weight_trunc)]*len(terminal))
	print(weight_trunc)
	print(tree_trunc)

	terminal_nodes = [graph[t] for t in terminal]
	for node in terminal_nodes:
		if node.paths.get(node): #has found the optimal paths for such node
			continue
		optimal_paths(graph, source=node,weight_trunc=weight_trunc)

	for node in graph.nodes.values():
		for source in terminal_nodes:
			if node.reduced.get(source):
				continue
			node.reduction(source) # do a reduction in terms of the number of objectives - decrease search space for the stars

	# check if they are connected pairwise - necessary but not sufficient condition for existence


	solution = []
	for central_node in graph.nodes.values():
		for tree_paths in itertools.product(*[central_node.paths[s] for s in terminal_nodes]):
			tree = Star(terminal=terminal_nodes,paths=[p for p in tree_paths],weight=WeightTree([p.weight for p in tree_paths]))
			
			if len(solution) == 0:
				solution += [tree]
			else:
				tree_sub = []
				flag_add = False
				for t in solution:
					if t.D(tree):
						flag_add = False
						break

					elif tree.D(t):
						tree_sub += [t]
						flag_add = True
					
					else:
						flag_add = True

				solution = [t for t in solution if t not in tree_sub]
				if flag_add:
					solution += [tree]

	#print(shortest_paths)
	return solution
