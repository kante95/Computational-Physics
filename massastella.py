from __future__ import division
import matplotlib.pyplot as plt
import matplotlib
# -*- coding: utf-8 -*-
import numpy as np
import math

matplotlib.rc('text', usetex=True)

n_range = np.arange(0,4,0.5) 

h = 0.001

mesh = np.arange(0,10,h)

xi1 = np.zeros(len(n_range))

temp =0

for n in n_range:
	theta = np.zeros(len(mesh))
	phi = np.zeros(len(mesh))

	#primo step
	phi[0] = 0
	theta[0] = 1
	#secondo step per evitare divisioni per zero
	k1 = 0
	l1 = h*(-(1**n))
	k2 = h*(l1/2)
	l2 = h*(-2*(l1/2)/(h*0.5)-(1+k1/2)**n)
	k3 = h*(l2/2)
	l3 = h*(-2*(l2/2)/(h*0.5)-(1+k2/2)**n)
	k4 = h*(l3)
	l4 = h*(-2*(l3)/(h) - (1+k3)**n)
	theta[1] = 1 + (k1+2*k2+2*k3+k4)/6
	phi[1] = (l1+2*l2+2*l3+l4)/6

	for i in range(1,len(mesh)-1):
		k1 = h*phi[i]
		l1 = h*(-2*phi[i]/(mesh[i])-theta[i]**n)
		k2 = h*(phi[i]+l1/2)
		l2 = h*(-2*(phi[i]+l1/2)/(mesh[i]+h*0.5)-(theta[i]+k1/2)**n)
		k3 = h*(phi[i]+l2/2)
		l3 = h*(-2*(phi[i]+l2/2)/(mesh[i]+h*0.5)-(theta[i]+k2/2)**n)
		k4 = h*(phi[i]+l3)
		l4 = h*(-2*(phi[i]+l3)/(mesh[i]+h) - (theta[i]+k3)**n)
		theta[i+1] = theta[i] + (k1+2*k2+2*k3+k4)/6
		phi[i+1] = phi[i] + (l1+2*l2+2*l3+l4)/6
		if theta[i+1]<0 or np.isnan(theta[i+1]) or np.isnan(phi[i+1]):
			break
		
	theta = theta[theta!=0]
	phi = phi[phi!=0]
	xi1[temp] = mesh[i+1]
	temp +=1
	#print phi,theta

	plt.figure(1)
	plt.plot(mesh[:len(theta)],theta)

	plt.figure(2)
	plt.plot(mesh[:len(phi)],phi)


plt.figure(1)
leg = ["n = "+ str(i) for i in n_range]
plt.legend(leg)
plt.xlabel(r"$\xi$")
plt.ylabel(r"$\theta(\xi)$")
plt.savefig("theta.png")

plt.figure(2)
leg = ["n = "+ str(i) for i in n_range]
plt.legend(leg)
plt.xlabel(r"$\xi$")
plt.ylabel(r"$\varphi(\xi)$")
plt.savefig("phi.png")


plt.figure(3)
plt.plot(n_range,xi1)

plt.show()