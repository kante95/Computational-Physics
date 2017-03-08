from __future__ import division
import matplotlib.pyplot as plt
import matplotlib
# -*- coding: utf-8 -*-
import numpy as np
import math

def f1(x,en):
	return 4*(1/(x**12) - 1/(x**6))-en

def f2(x,en):
	r = (4*(-1/(x**12) + 1/(x**6))+en)**(0.5)
	if(np.iscomplex(r) or np.isnan(r)):
		return 0
	else:
		return r

def Newton(start,energy,err):
	x = start
	while abs(0-f1(x,energy))>err:
		x = x - f1(x,energy)/(24*(1/(x**7) - 2/(x**13))) 
		
	return x

def F(E,n):
	gamma = 21.7
	xstart = 1.1224620483
	err_newton = 0.0001
	int_step = 0.0001

	xmin = Newton(xstart-0.1,E,err_newton) 
	xmax = Newton(xstart+0.1,E,err_newton)

	mesh = np.arange(xmin,xmax,int_step)

	tempsum = 0
	i = 1
	while i <len(mesh)-1:
		tempsum += 4*f2(mesh[i],E) + 2*f2(mesh[i+1],E)
		i+=2

	I = (int_step/3)*(tempsum+f2(mesh[0],E)+f2(mesh[-1],E))
	return gamma*I-np.pi*(n+0.5)

def secante(x1start,x2start,n,err):
	x1 = x1start
	x2 = x2start
	tempf1 = F(x1,n)
	tempf2 = F(x2,n)
	while abs(0-tempf2)>err:
		x = x2 - (x2-x1)*tempf2/(tempf2-tempf1)
		tempf1 = tempf2
		tempf2 = F(x,n)
		x1=x2
		x2=x
	return x2

E = -0.9
dE = 0.01
n = 5
err_secante = 0.0001

for i in range(n):
	print "Cerchiamo l'energia numero " + str(i)
	Etemp = E+dE
	xsec = F(E,i)
	xtry = F(Etemp,i)
	while xtry*xsec>0:
		Etemp = Etemp+dE
		xtry = F(Etemp,i)
	result = secante(E,Etemp,i,err_secante)
	print result
	E = result


