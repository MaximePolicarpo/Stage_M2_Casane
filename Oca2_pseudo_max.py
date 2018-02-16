#!/usr/bin/env python3


#On voit grace a ce script que c'est carrement stochastique le fait qu'un codon STOP apparaisse et atteigne la fixation

#Import of necessary modules

import sys
import random
import numpy as np
import numpy.random as npr
import time
import matplotlib as mpl
#mpl.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import matplotlib.pyplot as plt


seqLength= 2163 #Oca2_CDS=2163 bases

parser = argparse.ArgumentParser(description="Simulate gene loss by stop mutations")
parser.add_argument('-g', '--generations', type=int, required=True, help="Number of generations")        #Number of generations that will be simulated
parser.add_argument('-r', '--repeats', type=int, required=True, help="Number of repeats")		 #Number of simulations that will be runned
parser.add_argument('-m', '--mu', type=float, required=True, help="Mutation rate")			 #Mutation rate : estimated to be 1e-08 per base per generation for humans
parser.add_argument('-N', '--Ne', type=float, required=True, help="Effective population size")		#Effective population size (estimated to be between 100 and 500 (Fumey and Casane 2017, BioxRiv))
parser.add_argument('-s', '--stop', type=float, required=True, help="Frequency of stop mutations")       #Frequency of stop mutation (3,6% taking into account the codon composition of zebrafish opsins)
parser.add_argument('-l', '--tauxmigra', type=float, required=True,help="Combination of frequency of migration event and of migrant fish", default=0)    #migration frequency
parser.add_argument('-S', '--fitness', nargs='+', type=float, required=True, default=0, help="Fitness of mutated alleles")            #Define the selective value of the STOP allele
parser.add_argument('-H', '--dominance', type=float, required=True, default=0, help="Dominance of mutated alleles")	    #Define the dominance/recessive model
parser.add_argument('-T', '--frequenceP', type=float, required=True, default=0, help="Starting frequency of the pseudogene")
args = parser.parse_args()




X=args.fitness
fitness_list=list()
for i in X:
	fitness_list.append(float(i))




param = {"nbGeneration":args.generations,
"nbRepeats":args.repeats,
"txMut":args.mu,
"popSize":args.Ne,
"freqMutaStop":args.stop,
"txmigra":args.tauxmigra,
"h":args.dominance,
"freq":args.frequenceP
}

Matrice=list()


for fitness in fitness_list:
	generation_pseudogenisation=[0]*param["nbGeneration"]
	for j in range(param["nbRepeats"]):
		sys.stderr.write("\rDoing repeat {} on {}".format(j,param["nbRepeats"]))
	        sys.stderr.flush()
		freq=param["freq"]
		suivi_freq=list()
		suivi_nbGen=list()
		Liste_Evolution=list()
		Probabilite=0
		for n in range(param["nbGeneration"]):

				probMuta = param["txMut"] * seqLength * param["freqMutaStop"] * 2*param["popSize"] * (1-freq)
				if(random.random() < probMuta):
					freq += 1/(2*param["popSize"])


				averageSelection = ((1-freq)**2)*(1-param["txmigra"]) + param["txmigra"] + 2*freq*(1-freq)*(1-param["h"]*fitness)*(1-param["txmigra"]) + (freq**2)*(1-fitness)*(1-param["txmigra"])
				newFreq = ((freq**2)*(1-fitness)*(1-param["txmigra"]) + freq*(1-freq)*(1-param["h"]*fitness)*(1-param["txmigra"]))/averageSelection
				freq = npr.binomial(2*param["popSize"],newFreq)/(2*param["popSize"])
		
				if freq==1:
					for k in range((n-1), (param["nbGeneration"]+1)):
						generation_pseudogenisation[k-1]+=1.0
					break
	
		Probabilite=[x / float(param["nbRepeats"]) for x in generation_pseudogenisation]
	Matrice.append(Probabilite)



#Probabilite=[x / float(param["nbRepeats"]) for x in generation_pseudogenisation]
#plt.plot(range(param["nbGeneration"]), Probabilite)
#plt.show()

colors_used=("black", "red", "yellowgreen", "royalblue", "violet", "tomato", "chocolate", "lawngreen", "gray", "saddlebrown", "turquoise")
labels_list=list()
for fitness in fitness_list:
	labels_list.append(str(fitness))

labels_list=["S=" + x for x in labels_list]


with PdfPages("Probabilite_Oca2_Fixation.pdf") as pdf:
		for N in range((len(Matrice))):
			plt.plot(range(param["nbGeneration"]), Matrice[N], color=colors_used[N], label=labels_list[N])
		plt.legend()
		plt.ylim([-0.05, 1.05])
		plt.title("Prability for Oca2 to get fixed as pseudogene")
		pdf.savefig()
		plt.close()







#Si on veut tracer l'evolution de la frequence en fonction du nb de generation , ajouter cela au script:
			#Pour suivre l'evolution de la frequence dans le temps
#			suivi_freq.append(freq)
#			suivi_nbGen.append(n)
#			Liste_Evolution.append(suivi_freq)

			#Pour suivre la proba d'etre pseudogeniser en fonction de la generation
#plt.plot(suivi_nbGen, suivi_freq)











