#!/usr/bin/env  python
# -*- coding:utf8 -*-

"""
Auteur : Arthur ROBIEUX
E-mail : arthur.robieux@gmail.com
Date : 11/04/2017

Objectif : Fonction permettant de parser une protéine au format .pdb et de
la tester avec un fichier test

"""

import math, string


def parsePDB(infile) :
    """
    Permet de parser un fichier PDB.
    Prend en entrée un fichier pdb et renvoi en sortie le fichier parsé
    (organisé de façon à être utilisable dans d'autre programme).

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
            curres = "%s"%(line[22:26]).strip()
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
            dPDB[chain][curres][atomtype]["id"] = line[6:11].strip()

    return dPDB


def usage():
    """
    Permet de tester la fonction parsePDB avec un fichier test
    à entrer dans la fonction.
    """
    print parsePDB("arginine.pdb")


usage()