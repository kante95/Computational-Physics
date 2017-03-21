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

def particle_rint(p):
	return particle(np.round(p.x),np.round(p.y),np.round(p.z))

#This is a particle
class particle:
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
	def __str__(self):
		return "("+str(self.x)+","+str(self.y)+","+str(self.z)+")"
	def __add__(self, other):
		return particle(self.x+other.x,self.y+other.y,self.z+other.z)
	def __sub__(self, other):
		return particle(self.x-other.x,self.y-other.y,self.z-other.z)
	def __mul__(self, other):
		return particle(self.x*other,self.y*other,self.z*other)
	def __truediv__(self, other):
		return particle(self.x/other,self.y/other,self.z/other)
	def distance(self,other,L):
		dist = self-other - L*particle_rint((self-other)/L)
		return (dist.x**2+dist.y**2+dist.z**2)**0.5

def initialize(N,L):
	i = np.floor(nthroot(N,3))
	#what?!! three for cycles in a single line O.O
	pr = [particle(j - L/2,k - L/2,m - L/2) for j in np.arange(0,i+1,L/(i-1)) for k in np.arange(0,i+1,L/(i-1)) for m in np.arange(0,i+1,L/(i-1))]
	#RANDOM PARTICLES!
	for j in np.arange(N-i**3):
		pr.append(particle(np.random.uniform(-L/2, L/2),np.random.uniform(-L/2, L/2),np.random.uniform(-L/2, L/2)))
	return pr

# shouldn't I take the derivate :P?
def lennard_jones(r):
	return 4*(1/r**12-1/r**6)

def potential(particle,particle_istance,L):
	#shit, this will be challenging
	for i in particle_istance:
		dist = particle.distance(i,L)
		if dist<L/2:
			return lennard_jones(dist)


def simulation(num_particle,L, density, temperature time, step):
	particle = initialize(num_particle, L)
	
	#verlet algorithm
	mesh = np.arange(0,time,step)
	supermegaipermatrix = [initialize(num_particle,L) for i in mesh] #oh yeah mothafucka!
	np.insert(supermegaipermatrix,[0],particle)

	for t in range(len(mesh)-1):
		for j in range(len(particle)):
			new_coordinate = 2*supermegaipermatrix[j][t]-supermegaipermatrix[j][t-1]-(1/m)*potential(j,supermegaipermatrix[:][t],L) #will this create a new supermegaipermatrix in my memory? I hope it doesn't
			supermegaipermatrix[j][t+1] = new_coordinata - L*particle_rint(new_coordinata/L) #boundary condition


	#time avarage
	#P = 

	#return P

simulation(125,5,1,273,10,1)
