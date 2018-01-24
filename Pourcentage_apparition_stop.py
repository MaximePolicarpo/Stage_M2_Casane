#!/usr/bin/env python3

#Commandes bash:
#se mettre dans le dossier ou sont toutes les opsines du zebrafish :
# for n in * ; do  grep -v "^>" $n | tr -d '\n' | sed '/^$/d' | sed 's/ //g' > ../linearisation/$n.lin ; done
# cat *.lin > total
#Permet de lineariser toutes les séquences d'opsines en une sequence "geante"

with open("total" , "r") as mesopsines:
	seq=mesopsines.read()

print(len(seq))

n=0
x=0

codon=seq[n:n+3]

moyenne_seq=list()

while n!=len(seq):
	if codon == "CAA" or codon == "GAA" or codon == "AAA" or codon == "AGA" or codon == "CGA" or codon == "GGA" or codon == "TGG" or codon == "TGT" or codon == "TGC" or codon == "AAG" or codon == "GAG" or codon == "CAG" or codon == "TTG" or codon == "TCG":
		x+=1.0/3.0
		n+=3
		codon=seq[n:n+3]
	elif codon == "TTA" or codon == "TCA" or codon == "TAC" or codon == "TAT" or codon == "TCA" or codon == "TTA" or codon == "TGG" or codon == "TAC" or codon == "TAT":
		x+=2.0/3.0
		n+=3
		codon=seq[n:n+3]
	else :
		x+=0
		n+=3
		codon=seq[n:n+3]

z=x/n
print(n)
print(z)  #z= pourcentage des mutations qui peuvent faire apparaître un stop dans mes opsines.
