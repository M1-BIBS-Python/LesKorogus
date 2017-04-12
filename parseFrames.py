#!/usr/bin/env  python
# -*- coding:utf8 -*-

"""
Auteur : Arthur ROBIEUX
E-mail : arthur.robieux@gmail.com
Date : 12/04/2017

Objectif : Permet d'obtenir toutes les conformations sous un format lisible par Python grâce à parsePDB.
"""

from parsePDB import parsePDBchain
from parsePDB import parsePDBfile

#######################################################
#Récupération des données de la protéine de référence
#######################################################

parse_ref = parsePDBfile("pab21_ref.pdb")


##############################################
#Parsage des 500 conformations de la structure
##############################################

#Récupération des données du fichier contenant les 500 conformations

frames = open("pab21_500frames.pdb", "r")
lines = frames.readlines()
frames.close()

#Parsage

parse_frames = {}       #Dico contenant les sorties du parseur pour chaque conformations (clé : n°MODEL, valeur : sortie parsePDB())
nb_frames = 0           #Contient le nombre de conformations

conformations = []      #Contient le contenu du .pdb pour un seul model
model = ""              #Numéro du modèle

for line in lines:

    if "MODEL" in line:                                         #Si on commence une nouvelle conformations
        if conformations != "":                                 #On parse la conformation actuelle
            parse_frames[model] = parsePDBchain(conformations)

        conformations = []                                      #On supprime la conformation après l'avoir parsé
        model = line[10:14]                                     #Determination du nouveau "MODEL"
        nb_frames += 1

    else:
        conformations.append(line)                                 #Sinon on ajoute les lignes correspondant à la conformation

parse_frames[model] = parsePDBchain(conformations)                   #Dernière conformation


print parse_ref
print parse_frames
print nb_frames


