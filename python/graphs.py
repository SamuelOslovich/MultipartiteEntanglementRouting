import sys
import random as rd
import math
import numpy as np
import itertools

from weights import *

## ------------------------------------- GRAPH METHODS ------------------------------------- ##


class Node():
	def __init__(self, name=None, neighbours={}, prob_swap=1, t_mem=math.inf):
		self.name = name
		self.neighbours = neighbours
		self.prob_swap = prob_swap
		self.t_mem = t_mem
		self.paths = {}
		self.visited = {}
		self.reduced = {}

	def __str__(self):
		output = " ---- NODE " + str(self.name) + " ---- \n"
		for n,l in self.neighbours.items():
			output += " - Neighbour: " + str(n.name) + " | " + str(l) + "\n"
		output+= " Parameters: PROB: " + str(self.prob_swap) + " | T_MEM: " + str(self.t_mem) + "\n"
		return output

	def clear_paths(self):
		self.paths.clear()
		self.visited.clear()
		self.reduced.clear()

	def print_paths(self,source):
		print(" ---- NODE " + str(self.name) + " ---- ")
		if self.paths.get(source):
			for p in self.paths[source]:
				print(p)
		else:
			print("Empty path list!")

	def print_all_paths(self):
		print(" ---- NODE " + str(self.name) + " ---- ALL PATHS -----")
		for key in self.paths:
			print(" ---- NODE " + str(key.name) + " ---- ")
			self.print_paths(key)


	def add_neighbour(self,node,link):
		#print("adding link " + str(link.source.name) + "-" + str(link.target.name) + " to node " + str(node.name))
		self.neighbours[node] = link

	def lessthan(self,other,source):
		# Priority queue ordering
		if len(self.paths[source]) != 0 and len(other.paths[source]) != 0:
			for path_i in self.paths[source]:
				for path_j in other.paths[source]:
					if path_i.D(path_j): 
						return True

		elif len(self.paths[source]) == 0:
			return True
		else:
			return False

	def greaterthan(self,other,source):
		# Priority queue ordering
		if len(self.paths[source]) != 0 and len(other.paths[source]) != 0:
			for path_i in self.paths[source]:
				for path_j in other.paths[source]:
					if path_j.D(path_i): 
						return True

		elif len(other.paths[source]) == 0:
			return True
		else:
			return False

	def initialize(self):
		self.paths[self] = [Path(source=self)]

	def update_paths(self,paths_add,source):
		if not self.paths.get(source): #empty list of paths
			self.paths[source] = paths_add

		else:
			return self.merge_paths(paths_add,source)

	def merge_paths(self,paths_add,source):
		flag_change = False

		paths_sub = []

		for path_i in paths_add:
			for path_j in self.paths[source]:

				if path_j in paths_sub:
					continue
				## some path dominates the new one, immediately disregard

				if path_j.D(path_i):
					paths_sub += [path_i]
					continue

				#new path dominates one or more of the paths, erase them

				elif path_i.D(path_j):
					self.paths[source].remove(path_j)
					flag_change = True


		# in the end, merge two lists, as self.paths[source] only has paths not dominated by any path in paths_add 
		# and paths_add only has paths non-dominated by any path in self.paths[source]
		paths_add = [p for p in paths_add if p not in paths_sub]
		if len(paths_add) > 0:
			flag_change = True
			self.paths[source] += paths_add

		return flag_change

	def reduction(self,source):
		if self.reduced.get(source):
			return True
		if not self.paths.get(source):
			#print("No paths to reduce")
			self.reduced[source] = True
			return True
		
		for p in self.paths[source]:
			p.weight = WeightPath(p.weight)

		paths_sub = []
		for p1, p2 in itertools.combinations(self.paths[source], 2):
		
			if p1 in paths_sub or p2 in paths_sub:
				continue

			if p1.D(p2):
				paths_sub += [p2]
				continue
			
			if p2.D(p1):
				paths_sub += [p1]
				continue

		self.paths[source] = [p for p in self.paths[source] if p not in paths_sub]
		self.reduced[source] = True



class Link():
	def __init__(self, source=None, target=None, fid=1, prob_gen=1, time=0):
		self.fid = fid
		self.prob_gen = prob_gen
		self.time = time
		self.source = source
		self.target = target

	def __str__(self):
		return "Link: " + str(self.source.name) + "<->" + str(self.target.name) +  " || FID: " + str(self.fid) + " || PROB: " + str(self.prob_gen) + " || TIME: " + str(self.time) 

	def switch_st(self):
		return Link(source=self.target, target=self.source, fid=self.fid, prob_gen=self.prob_gen, time=self.time)

	def has(self,node):
		if self.source == node or self.target == node:
			return True
		else:
			return False


class Graph():
	def __init__(self, nodes={}, links={}):
		self.nodes = nodes
		self.links = links

	def __getitem__(self,key):
		if type(key) == int:
			return self.nodes[key]
		elif type(key) == tuple:
			return self.links[key]
		else:
			print("Error, graph has no object type " + str(type(key)))
			return False


	def add_node(self, node=Node()):
		self.nodes[node.name] = node

	def add_link(self, link=Link()):
		s = link.source
		t = link.target
		self.links[(s.name,t.name)] = link
		self.links[(t.name,s.name)] = link.switch_st()
		self.nodes[s.name].add_neighbour(t,link)
		self.nodes[t.name].add_neighbour(s,link.switch_st())

	def print(self):
		for key, node in self.nodes.items():
			print(node)

	def print_paths(self,source):
		for key, node in self.nodes.items():
			node.print_paths(source)

	def average_degree(self):
		average = len(self.links)/len(self.nodes)
		sigma = math.sqrt(sum( [(len(node.neighbours)-average)**2 for node in self.nodes.values()] )/len(self.nodes))
		return average, sigma

	def clear(self):
		for node in self.nodes.values():
			node.clear_paths()



class Path():
	def __init__(self, source=None, target=None, links=[], weight=Weight()):
		self.source = source
		self.target = target
		self.links = links
		self.weight = weight

	def __add__(self,link):
		nsource = self.source
		ntarget = link.target
		nlinks = self.links.copy()+ [link]
		nweight = self.weight + link
		return Path(source=nsource,target=ntarget,links=nlinks,weight=nweight)

	def __str__(self):
		output = "Path: " + str(self.source.name) 
		for l in self.links:
			output += "->" + str(l.target.name)	
		output += "\n"
		output += str(self.weight)
		return output

	
	def D(self,other):
		if self.weight.D(other.weight):
			return True
		else:
			return False

	def check_possible(self,weight_trunc):
		if self.weight < weight_trunc and self.no_loops():
			return True
		else:
			return False

	def no_loops(self):
		for l in self.links[:-1]:
			if l.has(self.target):
				return False

		return True

def reconstruct_path(graph,list_of_nodes):
	source = graph[list_of_nodes[0]]
	target = graph[list_of_nodes[-1]]
	links = []
	if len(list_of_nodes)==1:
		return Path(source,target,links,WeightPath())
	else:
		for i in range(0,len(list_of_nodes)-1):
			links += [graph[(list_of_nodes[i],list_of_nodes[i+1])]]
		weight = Weight()
		for link in links:
			weight = weight + link
		nweight = WeightPath(weight)
		print(nweight)
		return Path(source,target,links,nweight)


class Star():
	def __init__(self, terminal=[], paths=[], weight=WeightTree()):
		self.terminal = terminal
		self.paths = paths
		self.weight = weight

	def __add__(self,path):
		nterminal = self.terminal + [path.source]
		npaths = self.paths +[path]
		nweight = self.weight + path.weight
		return Star(terminal=nterminal, paths=npaths, weight=nweight)

	def __str__(self):
		output = "###########################\n"
		output+= "TWINKLE TWINKLE LITTLE STAR\n"
		output+= "---------------------------\n"
		for p in self.paths:
			output += str(p) + "\n"	
		output+= "---------------------------\n"
		output+= "---------------------------\n"
		output+= "Star Weight: " + str(self.weight) + "\n"
		output+= "###########################\n"

		return output

	
	def D(self,other):
		if self.weight.D(other.weight):
			return True
		else:
			return False

	def check_possible(self,weight_trunc):
		if self.weight < weight_trunc and self.no_cycles():
			return True
		else:
			return False

	def no_cycles(self):
		## still need to code this, prolly grows factorially!

		return True


	def get_center_node(self):
		return self.paths[0].target


def reconstruct_star(graph,list_of_paths):
	terminal = [path.source for path in list_of_paths]

	if len(list_of_paths) == 0:
		return Star(terminal,list_of_paths,WeightTree())

	print([path.weight for path in list_of_paths])
	weight = WeightTree(weights=[path.weight for path in list_of_paths])
	return Star(terminal,list_of_paths,weight)



