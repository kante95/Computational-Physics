from __future__ import division
import matplotlib.pyplot as plt
import matplotlib
# -*- coding: utf-8 -*-
import numpy as np
import math


#It sucks but i have to do it :(
def nthroot(x, n):
		r = 1
		for i in range(16):
		        r = (((n - 1) * r) + x / (r ** (n - 1))) / n
		return r

#This is a particle
class particle:
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
	def __str__(self):
		return "("+str(self.x)+","+str(self.y)+","+str(self.z)+")"

def initialize(N,L):
	i = np.floor(nthroot(N,3))
	#what?!! three for cycles in a single line O.O
	pr = [particle(j - L/2,k - L/2,m - L/2) for j in np.arange(0,i,L/i) for k in np.arange(0,i,L/i) for m in np.arange(0,i,L/i)]
	#RANDOM PARTICLES!
	for j in np.arange(N-i**3):
		pr.append(particle(np.random.uniform(-L/2, L/2),np.random.uniform(-L/2, L/2),np.random.uniform(-L/2, L/2)))
	return pr

def simulation(num_particle, density, time, step):
	particle = initialize(num_particle)

	#verlet algorithm

	#time avarage
	#P = 

	return P

a = initialize(125,5)
