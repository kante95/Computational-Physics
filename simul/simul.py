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
	def __init__(self, x, y, z, vx, vy, vz):
		self.position = np.array([x,y,z])
		self.speed = np.array([vx,vy,vz])
	def __str__(self):
		return self.position

def initialize(N,L):
	i = np.floor(nthroot(N,3))
	#what?!! three for cycles in a single line O.O
	pr = [particle(j - L/2,k - L/2,m - L/2,0,0,0) for j in np.arange(0,i+1,L/(i-1)) for k in np.arange(0,i+1,L/(i-1)) for m in np.arange(0,i+1,L/(i-1))]
	#RANDOM PARTICLES!
	for j in np.arange(N-i**3):
		pr.append(particle(np.random.uniform(-L/2, L/2),np.random.uniform(-L/2, L/2),np.random.uniform(-L/2, L/2),0,0,1))
	return pr

def lennard_jones_derivative(r):
	den = r[0]**2+r[1]**2 +r[2]**2
	return np.array([4*(-14*r[0]/den**7 + 8*r[0]/den**4),4*(-14*r[1]/den**7 + 8*r[1]/den**4),4*(-14*r[2]/den**7 + 8*r[2]/den**4)])

def forza(particle,state,L):
	total = 0
	for i in range(len(state)):
		if particle == i:
			continue
		r = (state[particle].position-state[i].position) - np.round((state[particle].position-state[i].position)/L)
		dist = (r[0]**2+r[1]**2 +r[2]**2)**0.5
		if dist<L/2:
			total += lennard_jones_derivative(r)
	return total 

def initialize_empty(N):
	return [particle(0,0,0,0,0,0) for i in range(N)]
	
def simulation(num_particles,L, density, temperature, time, step, integration="verlet"):
	print "Inizializzo lo stato iniziale...."
	state = initialize(num_particles, L)

	average = 0
	#velocity verlet algorithm
	mesh = np.arange(0,time,step)
	if integration=="velocity_verlet":
		for t in mesh:
			print "Simulando secondo "+str(t)+" s"
			newstate = initialize_empty(num_particles)
			accelleration = []
			for j in range(len(state)):
				accelleration.append(forza(j,state,L))
				newstate[j].position = state[j].position + state[j].speed*step + 0.5*accelleration[-1]*step**2
				newstate[j].position = newstate[j].position - L*np.round(newstate[j].position/L)
			for j in range(len(state)):
				newstate[j].speed = state[j].speed + 0.5*(accelleration[j] + forza(j,newstate,L))*step
			state = newstate
	elif integration=="verlet":
		
		newstate = initialize_empty(num_particles)
		current_state = state
		for i in range(len(state)):
			newstate[i].position = state[i].position + state[i].speed*step + 0.5*forza(i,state,L) 
		
		for t in mesh[1:]:
			print "Simulando secondo "+str(t)+" s"
			previous_state = current_state
			current_state = newstate
			newstate = initialize_empty(num_particles)
			for i in range(len(current_state)):
				newstate[i].position = 2*current_state[i].position - previous_state[i].position + forza(i,current_state,L)
				newstate[i].position = newstate[i].position - L*np.round(newstate[i].position/L)
				#print newstate[i].position



		#average =  what da fuck have i to put here??

	#time average
	#P = 

	#return P

simulation(125,5,1,273,3,0.1)
