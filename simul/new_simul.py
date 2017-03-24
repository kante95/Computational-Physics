# -*- coding: utf-8 -*-
import numpy as np

def nthroot(x, n):
		r = 1
		for i in range(16):
		        r = (((n - 1) * r) + x / (r ** (n - 1))) / n
		return r

def initialize(N,L):
	i = np.floor(nthroot(N,3))
	positions = np.array([[j - L/2,k - L/2,m - L/2] for j in np.arange(0,i+1,L/(i-1)) for k in np.arange(0,i+1,L/(i-1)) for m in np.arange(0,i+1,L/(i-1))])
	#RANDOM PARTICLES!
	for j in np.arange(N-i**3):
		positions.append([np.random.uniform(-L/2, L/2),np.random.uniform(-L/2, L/2),np.random.uniform(-L/2, L/2)])
	velocities = np.array([[0,0,0] for i in range(N)])
	return positions,velocities

def lennard_jones_derivative(r):
	den = r[0]**2+r[1]**2 +r[2]**2
	return -np.array([4*(-14*r[0]*den**-7 + 8*r[0]*den**-4),4*(-14*r[1]*den**-7 + 8*r[1]*den**-4),4*(-14*r[2]*den**-7 + 8*r[2]*den**-4)])

def calculate_accellerations(positions,N,L):
	accelerations = np.zeros((125,3))
	for i in range(N):
		for j in range(i+1,N):
			r = positions[i]-positions[j]-np.round((positions[i]-positions[j])/L)
			dist = (r[0]**2+r[1]**2 +r[2]**2)**0.5
			if dist<L/2:
				accelerations[i] += lennard_jones_derivative(r)
				accelerations[j] -= lennard_jones_derivative(r)
	return accelerations

def simulation(num_particles,L, density, temperature, time, step):
	print("Inizializzo lo stato iniziale....")
	positions,velocities = initialize(num_particles, L)

	#velocity verlet algorithm
	current_accelerations = calculate_accellerations(positions,num_particles,L)
	for t in np.arange(step,time,step):
		#print(positions)
		print("Simulando secondo "+str(t)+" s")
		#calcolo le nuove posizioni
		positions = positions + velocities*step + 0.5*current_accelerations*step**2
		positions = positions - L*np.round(positions/L)
		#accelerazioni dovute alle nuove posizioni
		new_accelerations = calculate_accellerations(positions,num_particles,L)
		#nuove velocità che si calcolano con le vecchie e le nuove accelerazioni
		velocities = velocities + 0.5*(current_accelerations + new_accelerations)*step
		#aggiorno le accelerazioni per il prossimo ciclo
		current_accelerations = new_accelerations

		#average =  what da fuck have i to put here??

	#time average
	#P = 
	#return P

simulation(125,5,1,273,3,0.01)
