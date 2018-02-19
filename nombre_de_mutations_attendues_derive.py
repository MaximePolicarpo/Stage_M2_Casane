#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas
import seaborn
import numpy
from matplotlib.backends.backend_pdf import PdfPages
import random


Sequence=("X", "X")*2



Nb_simulation=10
Ne=500.0
Nb_generation=200000


TOT=list()
for nb_simulation in range(Nb_simulation):
	Seq=list()
	for i in Sequence:
		Seq.append({"Base":i,"freq":0.0})
	for j in range(Nb_generation):
		for i in range(len(Seq)):
			Prob_muta=1e-08*48000.0*(1.0-Seq[i]["freq"])
			if random.random()<Prob_muta:
				Seq[i]["freq"]+=(1.0/(2*float(Ne)))
			Seq[i]["freq"]=numpy.random.binomial(2*float(Ne), Seq[i]["freq"])/float(2*Ne)
	b = [el["freq"] for el in Seq]
	TOT.append((b.count(1)))		


print(Seq)		
print(sum(TOT)/len(TOT))


#b = [el["freq"] for el in Seq]
#print(b.count(1))
