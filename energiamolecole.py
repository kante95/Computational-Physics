from __future__ import division
import matplotlib.pyplot as plt
import matplotlib
# -*- coding: utf-8 -*-
import numpy as np
import math


#funzione f per il metodo di newton
def f1(x,en):
	return 4*(1/(x**12) - 1/(x**6))-en

#funzione all'interno dell'integrale
def f2(x,en):
	return (4*(-1/(x**12) + 1/(x**6))+en)**(0.5)
	

#implementazione metodo di Newton
def Newton(start,energy,err):
	x = start
	while abs(0-f1(x,energy))>err:
		x = x - f1(x,energy)/(24*(1/(x**7) - 2/(x**13))) 		
	return x

#calcolo dell F
def F(E,n):
	gamma = 21.7
	xstart = 1.1224620483
	err_newton = 0.0001
	int_step = 0.0001

	# cerchiamo gli estremi di integrazione, uno a sinistra del minimo del potenziale, uno a destra
	xmin = Newton(xstart-0.1,E,err_newton) 
	xmax = Newton(xstart+0.1,E,err_newton)

	#calcolo dell'integrale
	mesh = np.arange(xmin,xmax,int_step)
	tempsum = 0
	i = 1
	while i <len(mesh)-1:
		tempsum += 4*f2(mesh[i],E) + 2*f2(mesh[i+1],E)
		i+=2

	I = (int_step/3)*(tempsum) #formula di Simpson, dove abbiamo tolto il valore iniziale e finale perche dovrebbero essere nulli
	return gamma*I-np.pi*(n+0.5)

#implementazione metodo della secante
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
err_secante = 0.00001
result = np.zeros(n)

#ciclo sui valori valori energetici
for i in range(n):
	Etemp = E+dE #tentativo di energia
	xsec = F(E,i) #energia di riferimento
	xtry = F(Etemp,i)
	while xtry*xsec>0:
		Etemp = Etemp+dE
		xtry = F(Etemp,i)
	result[i] = secante(E,Etemp,i,err_secante)
	print result[i]
	E = result[i]


#plot dello spettro energetico
plt.figure()

x = np.arange(0.95,2,0.01)
plt.plot(x,f1(x,0))

for i in result:
	plt.plot(x,[i]*len(x))

legend = ["Potenziale","n = 0","n = 1","n = 2","n = 3","n = 4"]

plt.ylabel("Potenziale V(x)")
plt.xlabel("x")
plt.legend(legend)
plt.savefig("spettro.png")
#plt.show()