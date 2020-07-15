import os
import sys

from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt


near = ["MG", "cation", "N", "amino", "O", "water", "lysarg"]

COLORS = ['silver', 'green', 'red', 'darkcyan', 'orange', 'green', 'purple',
		  'maroon', 'navy', 'brown', 'olive', 'indigo', 'fuchsia',
		  'lime', 'aqua', 'teal', 'gray', 'black']*5

def read_file(filepath):
	lines = list()
	with open(filepath, 'r') as infile:
		for line in infile:
			if not line.startswith("\t"):
				line = [field.strip() for field in line.split("\t")]
				values = list(map(float, line[:3]))
				# values[1] = values[1] if values[1] > 0 else values[1]+360
				truths = list(map(bool, map(int, line[3:10])))
				site = line[9]
				site_class = "".join(line[4:10])
				lines.append((values, truths, site, site_class))
	return lines

def plot(filepath, case):
	lines = read_file(filepath)
	for index in case:
		print(near[index])
	for i, j in [(0, 1), (1, 2), (0, 2)]:
		classes = defaultdict(lambda: ([list(), list()]))
		labels = {0: "None"}
		for values, truths, site, site_class in lines:
			for index_index, truth_index in enumerate(case):
				if truths[truth_index]:
					site_class = index_index+1
					labels[site_class] = near[truth_index]
					break
			else:
				site_class = 0
			classes[site_class][0].append(values[i])
			classes[site_class][1].append(values[j])
		order = sorted(classes, key=lambda site_class: -len(classes[site_class][0]))
		for site_class in order:
			color = COLORS[site_class]
			if site_class:
				plt.scatter(classes[site_class][0], classes[site_class][1], color=color, marker='+', label=labels[site_class])
			else:
				plt.scatter(classes[site_class][0], classes[site_class][1], color=color, marker='+')
		plt.legend()
		plt.show()

plot(sys.argv[1], list(map(int, sys.argv[2:])))

