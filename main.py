#!/usr/bin/env  python
# -*- coding:utf8 -*-

"""
Auteur : Arthur ROBIEUX
E-mail : arthur.robieux@gmail.com
Date : 12/04/2017

Objectif : Utilisation hiérarchique de nos fonctions afin de répondre au problème initial.
"""

from parsePDB   import *
from RMSD       import *
from interface  import *

import matplotlib.pyplot as plt

#########################################
# Parsage de la conformation de référence
#########################################

parse_ref = parsePDBfile("pab21_ref.pdb")

###############################################
# Parsage des 500 conformations de la structure
###############################################

parse_frames = parsePDBmulti("pab21_500frames.pdb")



############################
# Vérifications des parsages
############################

#len(parse_ref[" "].keys()) = 10 000 résidus
#len(parse_frames["5000"][" "].keys()) = 321 résidus

# commun = 0
#
# for a in parse_frames["5000"][" "].keys():
#     if a in parse_ref[" "].keys():
#         commun += 1
#
# print commun





################################################
# Calcul des RMSD pour chacune des conformations
################################################


RMSD = {}           #Clés : Model utilisé contre la référence, valeurs : RMSD

for model in parse_frames.keys():

    RMSD[model] = {}

    RMSD[model]["global"] = RMSDglobal(parse_ref, parse_frames[model])          # Calcul RMSD Global pour chaque conformation
    RMSD[model]["domainA1"] = RMSDlocal(parse_ref, parse_frames[model], "A1")   # Calcul RMSD local à A1
    RMSD[model]["domainA2"] = RMSDlocal(parse_ref, parse_frames[model], "A2")   # Calcul RMSD local à A2
    RMSD[model]["domainA3"] = RMSDlocal(parse_ref, parse_frames[model], "A3")   # Calcul RMSD local à A3
    RMSD[model]["domainA4"] = RMSDlocal(parse_ref, parse_frames[model], "A4")   # Calcul RMSD local à A4
    RMSD[model]["domainB"] = RMSDlocal(parse_ref, parse_frames[model], "B")     # Calcul RMSD local à B




##########################################
# Affichage RMSD dans un fichier de sortie
##########################################


# En tête

f = open("RMSD.txt", "w")
f.write("Model    Global    Domaine1    Domaine2    Domaine3    Domaine4    Domaine5\n\n")
f.close()

#Contenu

f = open("RMSD.txt", "a")

x = []
yG = []
yA1 = []
yA2 = []
yA3 = []
yA4 = []
yB = []

for model in sorted(map(int, RMSD.keys())):      #Pour chaque modèle (ordre croissant)

    #Affichage RMSD Globaux & locaux (pour chaque domaine)

    f.write(str(model))
    f.write("   ")
    f.write(str(RMSD[str(model)]["global"]))
    f.write("   ")
    f.write(str(RMSD[str(model)]["domainA1"]))
    f.write("   ")
    f.write(str(RMSD[str(model)]["domainA2"]))
    f.write("   ")
    f.write(str(RMSD[str(model)]["domainA3"]))
    f.write("   ")
    f.write(str(RMSD[str(model)]["domainA4"]))
    f.write("   ")
    f.write(str(RMSD[str(model)]["domainB"]))
    f.write("\n")

    #Stockage des valeurs pour faire des graphiques

    x.append(model)

    yG.append(RMSD[str(model)]["global"])
    yA1.append(RMSD[str(model)]["domainA1"])
    yA2.append(RMSD[str(model)]["domainA2"])
    yA3.append(RMSD[str(model)]["domainA3"])
    yA4.append(RMSD[str(model)]["domainA4"])
    yB.append(RMSD[str(model)]["domainB"])


f.close()


###Affichage des graphiques pour chaque RMSD (global & local)

plt.title("RMSD Global")
plt.plot(x, yG)
plt.xlabel('Model (Temps)')
plt.ylabel('RMSD Global')
plt.show()

plt.title("RMSD Domaine A1")
plt.plot(x, yA1)
plt.xlabel('Model (Temps)')
plt.ylabel('RMSD Domaine A1')
plt.show()

plt.title("RMSD Domaine A2")
plt.plot(x, yA2)
plt.xlabel('Model (Temps)')
plt.ylabel('RMSD Domaine A2')
plt.show()

plt.title("RMSD Domaine A3")
plt.plot(x, yA3)
plt.xlabel('Model (Temps)')
plt.ylabel('RMSD Domaine A3')
plt.show()

plt.title("RMSD Domaine A4")
plt.plot(x, yA4)
plt.xlabel('Model (Temps)')
plt.ylabel('RMSD Domaine A4')
plt.show()

plt.title("RMSD Domaine B")
plt.plot(x, yB)
plt.xlabel('Model (Temps)')
plt.ylabel('RMSD Domaine B')
plt.show()




#####################################
# Changements Conformationnels locaux
#####################################

#Test sur un fichier ayant que 2 conformations

parse_2frames = parsePDBmultiInt("2frames.pdb")

seuil = 9           #Seuil d'appartenance à l'interface
interface = {}

for model in parse_2frames.keys():

    print model

    for domain1 in parse_2frames[model][" "].keys():            #On compare 2 domaines

            domainARN = "B"

            if domain1 != domainARN:        #Pour éviter ARN vs ARN

                for residu1 in parse_2frames[model][" "][domain1]["reslist"]:   # On compare 2 résidus

                    inter = False

                    for residu2 in parse_2frames[model][" "][domainARN]["reslist"]:

                        if inter == False:

                            #Calcul de la distance entre 2 résidus
                            dist = distanceResidus(parse_2frames[model][" "][domain1][residu1], parse_2frames[model][" "][domainARN][residu2])

                            if dist < seuil:
                                if residu1 not in interface.keys():
                                    interface[residu1] = 1
                                    inter = True
                                else:
                                    interface[residu1] += 1
                                    inter = True



for residu in interface.keys():
    interface[residu] = float(interface[residu]) / 2



print interface







