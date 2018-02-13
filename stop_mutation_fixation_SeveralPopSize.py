#!/usr/bin/env python3

#Import of necessary modules

from cycler import cycler
import sys
import random
import numpy as np
import numpy.random as npr
import time
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse


#Define the argument that the user have to write down on the terminal.

parser = argparse.ArgumentParser(description="Simulate gene loss by stop mutations")
parser.add_argument('-g', '--generations', type=int, required=True, help="Number of generations")        #Number of generations that will be simulated
parser.add_argument('-r', '--repeats', type=int, required=True, help="Number of repeats")		 #Number of simulations that will be runned
parser.add_argument('-m', '--mu', type=float, required=True, help="Mutation rate")			 #Mutation rate : estimated to be 2,5e-08 per base per generation for humans
#parser.add_argument('-N', '--Ne', type=float, required=True, help="Effective population size")		#Effective population size (estimated to be between 100 and 500 (Fumey and Casane 2017, BioxRiv))
parser.add_argument('-s', '--stop', type=float, required=True, help="Frequency of stop mutations")       #Frequency of stop mutation (3,6% taking into account the codon composition of zebrafish opsins)
parser.add_argument('-f', '--freqmigra', type=float, required=False, help="Frequency of migration event", default=0)   #See below
parser.add_argument('-F', '--migrant', type=float, required=False, help="Frequency of migrant fish", default=0)        #See below
parser.add_argument('-l', '--tauxmigra', type=float, required=False,help="Combination of frequency of migration event and of migrant fish", default=0)    #See below
parser.add_argument('-S', '--fitness', type=float, required=False, default=0, help="Fitness of mutated alleles")            #Define if the broken gene (stop codon) brings an advantage or not
parser.add_argument('-H', '--dominance', type=float, required=False, default=0, help="Dominance of mutated alleles")	    #Define if a broken gene is dominant to the functionnal one

args = parser.parse_args()



s = raw_input("Liste des tailles de populations voulues, separees par des espaces: ")
Nes=map(int, s.split())
Ne_list=list()
for i in Nes:
	Ne_list.append(float(i))



#The user have the choice to either define a migration rate (=tauxmigra) or to set a frequency of events (for example there is 5% chance that a migration event occurs at each generation) coupled
#with a frequency of migrant fish (

if args.tauxmigra == 0 and ((args.freqmigra == 0 and args.migrant != 0) or (args.freqmigra != 0 and args.migrant == 0)):
    args.freqmigra = 0
    args.migrant = 0
    print("Assuming no migration")

if args.tauxmigra != 0 and (args.freqmigra != 0 or args.migrant != 0):
    print("ERROR : option -l can't be used together with -F and -f. Quitting")
    sys.exit()



#If everything is alright, we create a dictionnary called param that regroup all the arguments defined by the user

param = {"nbGeneration":args.generations,
"nbRepeats":args.repeats,
"txMut":args.mu,
"freqMutaStop":args.stop,
"probMigra":args.freqmigra,
"indMigra":args.migrant,
"txmigra":args.tauxmigra,
"fitness":args.fitness,
"h":args.dominance
}



#Write here the length of your genes. (CDS length). More genes mean more chance to get a stop mutation reaching fixation.

seqLength = (1050,1011,1056,1158,1059,1218,1056,966,861,945,1101,1221,1071,2802,1653,1227,2520,1062,1011,1146,1164,1227,2958,1197,1047,1881,891,882,1743)


#Create a dictionnary called listPos that will, for each gene, report its length, the mutation frequency and the number of mutation on it during the simulation. Youwill just have to write down
#initData(seqLength) to put the above list of genes in this new function creating a dictionnary.

def initData(seqLength):
    listePos = list()
    for i in seqLength:
        listePos.append({"seqSize":i,"freq":0.0,"nb":0}) #size,freq of the mutated alllele,nb muta
    return(listePos)


#The genetic drift is define as a binomial law. 2 possibles results : Picking the mutated allele (1) or picking the wt allele (0). 1 = Success 0= fail. 
#So the new frequency of the mutated allele will be the binomial result divided by 2*Ne (2 alleles * Population size)

def drift(freq,Ne):
    return npr.binomial(2*Ne,freq)/(2*Ne)



#Now we add the effect of natural selection and migration on our allele frquencies. freq=frequence of the STOP mutation allele.

def selection(freq,param,freqMigra):
    averageSelection = ((1-freq)**2)*(1-freqMigra) + freqMigra + 2*freq*(1-freq)*(1-param["h"]*param["fitness"])*(1-freqMigra) + (freq**2)*(1-param["fitness"])*(1-freqMigra)
    newFreq = ((freq**2)*(1-param["fitness"])*(1-freqMigra) + freq*(1-freq)*(1-param["h"]*param["fitness"])*(1-freqMigra))/averageSelection
    return(newFreq)


#Now we define the evolution of our allele -> First, the selection and mirgation act on our allele frequencies. Then we put a probability of STOP apparition which is equal to the probability of
#mutation * probability of themutation iducing a stop codon * 2N (2gene per individuals) * (1-freq of the non mutated allele). Finally those frequencies evolve thanks to genetic drift. 

def evolve(pos, param, freqMigra, Ne):
    pos["freq"] = selection(pos["freq"], param, freqMigra)
	
    probMuta = param["txMut"] * pos["seqSize"] * param["freqMutaStop"] * 2*Ne * (1-pos["freq"])
    if(random.random() < probMuta):
        pos["freq"] += 1/(2*Ne)
        pos["nb"] += 1

    pos["freq"] = drift(pos["freq"],Ne)
    return(pos)


#We put a counter that will count the number of alleles with a STOP mutation that got fixed. 

def mutationCounter(listeFreq):
    return listeFreq.count(1)



#We create a matrix which will, every 10 generations, store the mutated allele frequency (for each of genes)
countMutaMatrixF=list()

for Ne in Ne_list:
	dimension = (((param["nbGeneration"]//10)),(len(seqLength)+1))
	countMutaMatrix = np.zeros(dimension)
	tps1 = time.clock()
	for j in range(param["nbRepeats"]):
 	    sys.stderr.write("\rDoing repeat {} on {}".format(j,param["nbRepeats"]))
	    sys.stderr.flush()
	    listePos = initData(seqLength)
	    for i in range(param["nbGeneration"]):
	        freqMigra = 0.0
	        if param["txmigra"] != 0:
	            freqMigra = npr.binomial(2*Ne,param["txmigra"])/(2*Ne)
	        elif random.random() < param["probMigra"]:
	            freqMigra = param["indMigra"]

	        for pos in listePos:
	            pos = evolve(pos, param, freqMigra, Ne)

	        if(i % 10 == 0):
	            #print(i)
	            b = [el["freq"] for el in listePos]
	            countMutaMatrix[(i//10),mutationCounter(b)] += 1
	countMutaMatrixF.append(countMutaMatrix)


print(countMutaMatrixF)


#Make a graph with the matrix and save it. 

#with PdfPages("graphe_simu_stop_Several_Ne_{}_generations_{}_repeats_{}_mu_{}_migra_{}_migrants_{}_mig_{}_h_{}_s.pdf".format(param["nbGeneration"], param["nbRepeats"], param["txMut"],param["probMigra"], param["indMigra"],param["txmigra"], param["h"], param["fitness"])) as pdf:
#    for n in range((len(countMutaMatrixF))):
#	for i in range((len(seqLength)+1)):
#	        plt.plot(countMutaMatrixF[n][:,i]/param["nbRepeats"], 'k-')
#		plt.plot()
#		plt.ylim([-0.05,1.05])
#		plt.title("Probability of {} pseudogenes".format(i))
#		pdf.savefig()
#		plt.close()
	



#METTRE AUTANT DE COULEUR ET LABELS QUIL Y A DE TAILLES DE POPULATION DIFFERENTES !!!!!!



colors_used=("black", "red", "yellowgreen", "royalblue", "violet", "tomato", "chocolate", "lawngreen", "gray", "saddlebrown", "turquoise")
labels_list=list()
for ne in Ne_list:
	labels_list.append(str(ne))

labels_list=["Ne=" + x for x in labels_list]


with PdfPages("graphe_simu_stop_Several_Ne_{}_generations_{}_repeats_{}_mu_{}_migra_{}_migrants_{}_mig_{}_h_{}_s.pdf".format(param["nbGeneration"], param["nbRepeats"], param["txMut"],param["probMigra"], param["indMigra"],param["txmigra"], param["h"], param["fitness"])) as pdf:
	for i in range((len(seqLength)+1)):
		for N in range((len(countMutaMatrixF))):
			plt.plot(countMutaMatrixF[N][:,i]/param["nbRepeats"], 'k-', color=colors_used[N], label=labels_list[N])
#		plt.plot(countMutaMatrixF[1][:,i]/param["nbRepeats"], 'k-')
#		plt.plot(countMutaMatrixF[1][:,i]/param["nbRepeats"], 'k-')
#		plt.plot(countMutaMatrixF[1][:,i]/param["nbRepeats"], 'k-')
		plt.legend()
		plt.ylim([-0.05, 1.05])
		plt.title("Probability of {} pseudogenes".format(i))
		pdf.savefig()
		plt.close()





np.savetxt("out_simu_stop_{}_generations_{}_repeats_{}_mu_{}_migra_{}_migrants_{}_mig_{}_h_{}_s.csv".format(param["nbGeneration"], param["nbRepeats"], param["txMut"], param["probMigra"], param["indMigra"],param["txmigra"], param["h"], param["fitness"]), countMutaMatrixF, delimiter=",")

