#!/usr/bin/env  python
# -*- coding:utf8 -*-

"""
Auteur : Arthur ROBIEUX
E-mail : arthur.robieux@gmail.com
Date : 12/04/2017

Objectif : Permet de calculer la RMSD entre 2 protéines.
"""

from math import sqrt


def RMSDglobal(parse1, parse2):
    """
    Fonction permettant de calculer la RMSD entre 2 protéines.

    Input : 2 protéines parsées
    Output : RMSD correspondant à ces 2 protéines
    """

    N = float(0)
    somme = float(0)

    for residu1 in parse1[" "].keys():              #Pour chaque résidus dans les 2 protéines
        for residu2 in parse2[" "].keys():

            if (residu1 == residu2 and residu2 != 'reslist'):                #Si on est sur le même résidus

                for atom in parse2[" "][residu2]["atomlist"]:

                    somme += (sqrt((parse2[" "][residu2][atom]['x'] - parse1[" "][residu1][atom]['x']) ** 2 +
                             (parse2[" "][residu2][atom]['y'] - parse1[" "][residu1][atom]['y']) ** 2 +
                             (parse2[" "][residu2][atom]['z'] - parse1[" "][residu1][atom]['z']) ** 2)) ** 2
                    N += 1

    RMSD = sqrt(1 / N * (somme))

    return RMSD





def RMSDlocal(parse1, parse2, domain):
    """
    Fonction permettant de calculer la RMSD entre 2 protéines en un domaine.

    Input : 2 protéines parsées et le domaine souhaité
    Output : RMSD correspondant à ces 2 protéines dans ce domaine
    """

    N = float(0)
    somme = float(0)

    for residu1 in parse1[" "].keys():              #Pour chaque résidus dans les 2 protéines
        for residu2 in parse2[" "].keys():

            if (residu1 == residu2 and residu2 != 'reslist'):               #Si on est sur le même résidus

                for atom in parse2[" "][residu2]["atomlist"]:

                    if parse2[" "][residu2][atom]['domain'] == domain:      #Si l'atome est dans le domaine voulu

                        somme += sqrt((parse2[" "][residu2][atom]['x'] - parse1[" "][residu1][atom]['x']) ** 2 +
                                 (parse2[" "][residu2][atom]['y'] - parse1[" "][residu1][atom]['y']) ** 2 +
                                 (parse2[" "][residu2][atom]['z'] - parse1[" "][residu1][atom]['z']) ** 2)
                        N += 1



    RMSD = sqrt(1 / N * (somme))

    return RMSD




