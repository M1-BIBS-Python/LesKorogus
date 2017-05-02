#!/usr/bin/env  python
# -*- coding:utf8 -*-

"""
Auteur : Arthur ROBIEUX - Julien ROZIERE
E-mail : arthur.robieux@gmail.com, julien.roziere@u-psud.fr
Date : 12/04/2017

Objectif : Utilisation hiérarchique de nos fonctions afin de répondre au problème initial.
"""

from parsePDB   import *
from RMSD       import *
from interface  import *


#########################################
# Parsage de la conformation de référence
#########################################

parse_ref = parsePDBfile("pab21_ref.pdb")

###############################################
# Parsage des 500 conformations de la structure
###############################################

parse_frames = parsePDBmulti("pab21_500frames.pdb")


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
# Obtention de graphique
##########################################

sortieRMSD(RMSD)

"""

#####################################
# Changements Conformationnels locaux
#####################################

#Parsage spécifique pour cette partie (changement de hiérarchie)

parse_framesInt = parsePDBmultiInt("pab21_500frames.pdb")

seuil = 9           #Seuil d'appartenance à l'interface
interface = {}      #Clé = Résidu appartenant à l'interface, Valeur : Fréquence
domainARN = "B"

for model in parse_framesInt.keys():
    for domain1 in parse_framesInt[model][" "].keys():            #On compare 2 domaines

            if domain1 != domainARN:        #Pour éviter ARN vs ARN

                for residu1 in parse_framesInt[model][" "][domain1]["reslist"]:   # On compare 2 résidus

                    inter = False           #Booléen permettant de savoir si le résidu est deja dans l'interface

                    for residu2 in parse_framesInt[model][" "][domainARN]["reslist"]:

                        if inter == False:          #Si le résidu n'est pas deja dans l'interface

                            #Calcul de la distance entre 2 résidus
                            dist = distanceResidus(parse_framesInt[model][" "][domain1][residu1], parse_framesInt[model][" "][domainARN][residu2])

                            if dist < seuil:                            #Si distance inférieur au seuil, résidu dans interface
                                if residu1 not in interface.keys():     #puis on change de résidu
                                    interface[domain1] = {}
                                    interface[domain1][residu1] = 1
                                    inter = True
                                else:
                                    interface[domain1][residu1] += 1
                                    inter = True

for domain in interface.keys():
    for residu in interface[domain].keys():                                   #On calcul la fréquence de la présence dans l'interface
        interface[domain][residu] = float(interface[domain][residu]) / 500    #pour chaque résidu présent 1x au moins


###Affichage des fréquences des résidus dans l'interface

f = open("interface.txt", "w")
for domain in interface.keys():
    f.write(domain)
    f.write("\n")
    for residu in sorted(map(int, interface[domain].keys())):
        f.write(str(residu))
        f.write("   ")
        f.write(str(interface[domain][str(residu)]))
        f.write("\n")
f.close()



####################################################################
# Calcul des temps de contacts pour les paires de résidus d'intérêts
####################################################################

tpsContact = {}             #Clé = paire de résidus, Valeurs = Temps de Contact
tpsContact["26-38"] = 0     #Initialisation du dico
tpsContact["25-46"] = 0
tpsContact["33-34"] = 0
tpsContact["31-100"] = 0

for model in parse_framesInt.keys():      #Calcul des Temps de Contact pour chaque paire pour chaque conformation
    if(tempsContact(parse_framesInt[model][" "]["B"]["26"], parse_framesInt[model][" "]["A1"]["38"], seuil)):
        tpsContact["26-38"] += 1

    if(tempsContact(parse_framesInt[model][" "]["B"]["25"], parse_framesInt[model][" "]["A4"]["46"], seuil)):
        tpsContact["25-46"] += 1

    if(tempsContact(parse_framesInt[model][" "]["B"]["33"], parse_framesInt[model][" "]["A3"]["34"], seuil)):
        tpsContact["33-34"] += 1

    if(tempsContact(parse_framesInt[model][" "]["B"]["31"], parse_framesInt[model][" "]["A4"]["100"], seuil)):
        tpsContact["31-100"] += 1



###Affichage des temps de Contact

f = open("contact.txt", "w")
for paire in tpsContact.keys():
    f.write(paire)
    f.write("   ")
    f.write(str(tpsContact[paire]))
    f.write("\n")
f.close()
"""