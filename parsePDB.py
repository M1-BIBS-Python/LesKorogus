#!/usr/bin/env  python
# -*- coding:utf8 -*-

"""
Auteur : Arthur ROBIEUX
E-mail : arthur.robieux@gmail.com
Date : 11/04/2017

Objectif : Fonction permettant de parser une protéine au format .pdb afin
de pouvoir l'utiliser dans des sripts Python.
"""

import math, string


def parsePDBfile(infile) :
    """
    Permet de parser un fichier PDB.

    Input : Fichier pdb représentant une protéine.
    Output : Protéine parsée, lisible par Python.

    """

    # lecture du fichier PDB 
    f = open(infile, "r")
    lines = f.readlines()
    f.close()


    # var init
    chaine = True
    firstline = True
    prevres = None
    dPDB = {}
    dPDB["reslist"] = []
    dPDB["chains"] = []
    
    # parcoure le PDB   
    for line in lines :
        if line[0:4] == "ATOM" :
            chain = line[21]
            if not chain in dPDB["chains"] :
                dPDB["chains"].append(chain)
                dPDB[chain] = {}
                dPDB[chain]["reslist"] = []
            curres = "%s"%(line[22:27]).strip()
            if not curres in dPDB[chain]["reslist"] :
                dPDB[chain]["reslist"].append(curres)
                dPDB[chain][curres] = {}
                dPDB[chain][curres]["resname"] = string.strip(line[17:20])
                dPDB[chain][curres]["atomlist"] = []
            atomtype = string.strip(line[12:16])
            dPDB[chain][curres]["atomlist"].append(atomtype)
            dPDB[chain][curres][atomtype] = {}
            #print "cures ", curres
            #print dPDB[chain][curres]
 
            dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
            dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
            dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
            dPDB[chain][curres][atomtype]["domain"] = line[72:74].strip()

    return dPDB



def parsePDBstring(lines):
    """
    Permet de parser une chaine de caractère contenant le contenu d'un fichier pdb.

    Input : Chaine de caractère contenant une protéine au format pdb.
    Output : Protéine parsée, lisible par Python.

    """

    # var init
    chaine = True
    firstline = True
    prevres = None
    dPDB = {}
    dPDB["reslist"] = []
    dPDB["chains"] = []

    # parcoure le PDB
    for line in lines:
        if line[0:4] == "ATOM":
            chain = line[21]
            if not chain in dPDB["chains"]:
                dPDB["chains"].append(chain)
                dPDB[chain] = {}
                dPDB[chain]["reslist"] = []
            curres = "%s" % (line[22:27]).strip()
            if not curres in dPDB[chain]["reslist"]:
                dPDB[chain]["reslist"].append(curres)
                dPDB[chain][curres] = {}
                dPDB[chain][curres]["resname"] = string.strip(line[17:20])
                dPDB[chain][curres]["atomlist"] = []
            atomtype = string.strip(line[12:16])
            dPDB[chain][curres]["atomlist"].append(atomtype)
            dPDB[chain][curres][atomtype] = {}
            # print "cures ", curres
            # print dPDB[chain][curres]

            dPDB[chain][curres][atomtype]["x"] = float(line[30:38])
            dPDB[chain][curres][atomtype]["y"] = float(line[38:46])
            dPDB[chain][curres][atomtype]["z"] = float(line[46:54])
            dPDB[chain][curres][atomtype]["domain"] = line[72:74].strip()

    return dPDB



def parsePDBmulti(infile):
    """
    Permet à partir d'un fichier contenant plusieurs proteines au format pdb, de parser chacune d'entre elles.

    Input : Fichier contenant plusieurs protéines au format pdb.
    Output : Dictionnaire contenant chaque protéine (clé) et son format parsé (lisible par Python).

    """
    frames = open(infile, "r")
    lines = frames.readlines()
    frames.close()

    parse_frames = {}       #Dico contenant les sorties du parseur pour chaque conformations (clé : n°MODEL, valeur : sortie parsePDB())
    nb_frames = 0           #Contient le nombre de conformations

    conformations = []      #Contient le contenu du .pdb pour un seul model
    model = ""              #Numéro du modèle (temps T de l'observation)

    for line in lines:

        if "MODEL" in line:                                         #Si on commence une nouvelle conformations
            if model != "":                                         #On parse la conformation actuelle
                parse_frames[model] = parsePDBstring(conformations)

            conformations = []                                      #On supprime la conformation après l'avoir parsé
            model = line[10:14].strip()                             #Determination du nouveau "MODEL"
            nb_frames += 1

        else:
            conformations.append(line)                              #Sinon on ajoute les lignes correspondant à la conformation

    parse_frames[model] = parsePDBstring(conformations)             #Dernière conformation

    return parse_frames

