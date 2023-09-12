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
import matplotlib.pyplot as plt
import time
timestr = time.strftime("%Y%m%d")
print(timestr)

plot_type = sys.argv[1]
dist_type = sys.argv[2]

n_colors = 4
cmap = plt.get_cmap('twilight', n_colors)
colors = cmap(np.linspace(0,1,n_colors))
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 12

plt.figure(figsize=(8,6))

plt.xscale('linear')
plt.xlabel('Size of Star')
plt.xlim([2,15])

plt.hlines(y=1, xmin=2, xmax=15, linewidth=1, linestyle='-.',color="black",alpha=0.7)

#plt.ylim([0.08,1-10**(-6)])
#plt.subplots_adjust(left=0.15, bottom=None, right=0.98, top=0.98, wspace=None, hspace=None)
#plt.yscale('log')

name = ""

plt.yscale("linear")
if plot_type == "rate":
	plt.ylabel('$\text{Rate}_{bound}/\text{Rate}_{opt}$')
	name = "rate"
	if dist_type == "1":
		name += "1"
	elif dist_type == "2":
		name += "2"
#plt.ylim([0.9825, 1.001])


elif plot_type == "fid":
	plt.ylabel('\text{Fid}_{bound}/\text{Fid}_{opt}')
	name = "fid"
	if dist_type == "1":
		name += "1"
	elif dist_type == "2":
		name += "2"

plt.legend(loc=3)
# Show the major grid lines with dark grey lines
plt.grid(b=True, which='major', color='#666666', linestyle='-' , alpha=0.7)
# Show the minor grid lines with very faint and almost transparent grey lines
#plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
#plt.minorticks_on()


plt.savefig( sys.path[0] + "/Plots/" + name + "_" + timestr '.pdf',dpi=300,transparent=True)
plt.show()
