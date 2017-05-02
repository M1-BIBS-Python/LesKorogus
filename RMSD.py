#!/usr/bin/env  python
# -*- coding:utf8 -*-

"""
Auteur : Arthur ROBIEUX
E-mail : arthur.robieux@gmail.com
Date : 12/04/2017
Objectif : Permet de calculer la RMSD entre 2 protéines.
"""

from math import sqrt
import matplotlib.pyplot as plt



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

                    somme += ((parse2[" "][residu2][atom]['x'] - parse1[" "][residu1][atom]['x']) ** 2 +
                             (parse2[" "][residu2][atom]['y'] - parse1[" "][residu1][atom]['y']) ** 2 +
                             (parse2[" "][residu2][atom]['z'] - parse1[" "][residu1][atom]['z']) ** 2)
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

                        somme += ((parse2[" "][residu2][atom]['x'] - parse1[" "][residu1][atom]['x']) ** 2 +
                                 (parse2[" "][residu2][atom]['y'] - parse1[" "][residu1][atom]['y']) ** 2 +
                                 (parse2[" "][residu2][atom]['z'] - parse1[" "][residu1][atom]['z']) ** 2)
                        N += 1



    RMSD = sqrt(1 / N * (somme))

    return RMSD




def sortieRMSD(RMSD):
    """

    :param RMSD: Dictionnaire contenant les RMSD
    :return: 0
    """
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

    graphRMSD(x, yG, yA1, yA2, yA3, yA4, yB)



def graphRMSD(x, yG, yA1, yA2, yA3, yA4, yB):
    """
    Affiche les graphiques des RMSD.

    :param x: abscisse
    :param yG: ordonnée global
    :param yA1: ordonnée domaine A1
    :param yA2: ordonnée domaine A2
    :param yA3: ordonnée domaine A3
    :param yA4: ordonnée domaine A4
    :param yB: ordonnée domaine B
    :return: 0
    """

    ###Affichage des graphiques pour RMSD global

    plt.plot(x, yG, 'r')
    plt.xlabel('Conformations')
    plt.ylabel('RMSD')
    plt.title('RMSD Global de la proteine')
    plt.show()

    ###Affichage des graphiques pour RMSD locaux

    plt.title("RMSD Local de chaque domaine")

    plt.plot(x, yA1, 'g', label = "A1")
    plt.plot(x, yA2, 'r', label = "A2")
    plt.plot(x, yA3, 'y', label = "A3")
    plt.plot(x, yA4, 'b', label = "A4")

    plt.show()