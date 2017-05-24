import numpy as np
import matplotlib.pyplot as plt

def read_data(file):
	#Open the file and read from line 2
	data = np.genfromtxt(file, delimiter=",",usecols = (0,1))
	#eliminate rows with nan
	data = data[~np.isnan(data).any(axis=1)]
	#extract columns
	pression = data[:,0]
	volume = data[:,1]
	# return data between 350 and 750 nm
	return pression,volume


for i in [1,2,3,4,5]:
	pression,volume = read_data("./ultimi/simulation"+str(i)+".txt")
	plt.plot(volume,pression)
	plt.xscale('log')
	plt.yscale('log')
plt.show()

