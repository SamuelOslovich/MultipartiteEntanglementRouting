import sys
import random as rd
import math
import numpy as np

def prod(list_x):
	out = list_x[0]
	for x in list_x[1:]:
		out = out*x
	return out

## ------------------------------------- METRICS AND WEIGHTS ------------------------------------- ##


class Weight():   
	def __init__(self, fid=1, prob=1, time=0, sigma=math.inf):
		self.fid = fid
		self.prob = prob
		self.time = time
		self.sigma = sigma

	def __add__(self,link):
		probf = self.prob * link.prob_gen * link.source.prob_swap
		sigmaf = 0
		if (1/self.sigma+1/link.source.t_mem+1/link.target.t_mem) == 0:
			sigmaf = math.inf
		else:
			sigmaf = 1/(1/self.sigma+1/link.source.t_mem+1/link.target.t_mem)
		fidf = self.fid * link.fid
		timef = self.time + link.time
		return Weight(fidf,probf,timef,sigmaf)

	def __gt__(self,other):
		if self.fid < other.fid:
			return True
		elif self.fid == other.fid and self.prob < other.prob:
			return True
		elif self.fid == other.fid and self.prob < other.prob and time > other.time:
			return True
		elif self.fid == other.fid and self.prob < other.prob and time > other.time  and self.sigma < other.sigma:
			return True
		else:
			return False

	def __lt__(self,other):
		if self.fid >= other.fid and self.prob >= other.prob and self.time <= other.time and self.sigma >= other.sigma:
			return True
		else:
			return False

	def D(self,other):
		if self.fid == other.fid and self.prob == other.prob and self.time == other.time and self.sigma == other.sigma:
			return True
		elif  (( self.fid > other.fid and self.prob >= other.prob and self.time <= other.time and self.sigma >= other.sigma ) or 
			( self.fid >= other.fid and self.prob > other.prob and self.time <= other.time and self.sigma >= other.sigma ) or
			( self.fid >= other.fid and self.prob >= other.prob and self.time < other.time and self.sigma >= other.sigma ) or 
			( self.fid >= other.fid and self.prob >= other.prob and self.time <= other.time and self.sigma > other.sigma )) :
			return True
		else: 
			return False

	def __str__(self):
		return "Dist: " + str(self.fid) + " || PROB: " + str(self.prob) + " || TIME: " + str(self.time) + " || SIGMA: " + str(self.sigma) 


class WeightPath():
	def __init__(self, weight=Weight()):
		if weight.sigma == 0:
			self.fid = weight.fid
		else:
			self.fid = weight.fid*math.exp(-weight.time/weight.sigma)
		self.prob = weight.prob
		self.time = weight.time

	def __gt__(self,other):
		if self.fid < other.fid:
			return True
		elif self.fid == other.fid and self.prob < other.prob:
			return True
		elif self.fid == other.fid and self.prob < other.prob and time > other.time:
			return True
		else:
			return False

	def __lt__(self,other):
		if self.fid >= other.fid and self.prob >= other.prob and self.time <= other.time:
			return True
		else:
			return False

	def D(self,other):
		if self.fid == other.fid and self.prob == other.prob and self.time == other.time:
			return True
		elif  (( self.fid > other.fid and self.prob >= other.prob and self.time <= other.time ) or 
			( self.fid >= other.fid and self.prob > other.prob and self.time <= other.time  ) or
			( self.fid >= other.fid and self.prob >= other.prob and self.time < other.time  ) ): 
			return True
		else: 
			return False

	def __str__(self):
		return "FID: " + str(self.fid) + " || PROB: " + str(self.prob) + " || TIME: " + str(self.time) 



class WeightTree():
	def __init__(self, weights=[WeightPath()]):
		self.weights = weights
		self.fids = [prod([(1+w.fid)/2 for w in weights]),
					 prod([(1-w.fid)/2 for w in weights]),
					 prod([w.fid for w in weights])]
		self.time = max([w.time for w in weights])
		self.prob = prod([w.prob for w in weights])
		#metrics
		self.fid = 1/2 * (sum(self.fids))
		if self.time == 0:
			self.rate = math.inf
		else:
			self.rate = self.prob / self.time 

	def __add__(self,path):
		return WeightTree(weights=self.weights + [path.weight])

	def __gt__(self,other):
		if self.fid < other.fid:
			return True
		elif self.fid == other.fid and self.rate < other.rate:
			return True
		else:
			return False

	def __lt__(self,other):
		if self.fid >= other.fid and self.rate >= other.rate:
			return True
		else:
			return False

	def D(self,other):
		if self.fid == other.fid and self.rate == other.rate:
			return True
		elif  (( self.fid > other.fid and self.rate >= other.rate ) or 
			   ( self.fid >= other.fid and self.rate > other.rate ) ): 
			return True
		else: 
			return False

	def __str__(self):
		return "FID: " + str(self.fid) + " || RATE: " + str(self.rate) 



class Capacity():
	def __init__(self, cap=math.inf, length=0):
		self.capacity = cap
		self.length = length

	def __add__(self,other):
		return Capacity(cap=min(self.capacity,other.capacity),length=self.length+other.length)

	def __lt__(self,other):
		if self.capacity > other.capacity:
			return True
		elif self.capacity == other.capacity and self.length < other.length:
			return True
		else:
			return False

	def __str__(self):
		return f"<Capacity:{self.capacity} || Length:{self.length}>"


	def __repr__(self):
		return f"<Capacity:{self.capacity} | Length:{self.length}>"


class TreeCapacity():
	def __init__(self, cap=Capacity()):
		self.capacity = cap.capacity
		self.length = cap.length
		self.maximum = cap.length

	def __add__(self,other):
		out = TreeCapacity()
		out.capacity = min(self.capacity,other.capacity)
		out.length = self.length+other.length
		out.maximum = max(self.length,other.length)
		return out

	def __lt__(self,other):
		if self.capacity > other.capacity:
			return True
		elif self.capacity == other.capacity and self.maximum < other.maximum:
			return True
		elif self.capacity == other.capacity and self.maximum == other.maximum and self.length < other.length:
			return True
		else:
			return False

	def __str__(self):
		return f"<Capacity:{self.capacity} || Maximum:{self.maximum} || Length:{self.length}>"


	def __repr__(self):
		return f"<Capacity:{self.capacity} | Maximum:{self.maximum} | Length:{self.length}>"


