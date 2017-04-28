#!/usr/bin/env  python
# -*- coding:utf8 -*-

"""
Auteur : Arthur ROBIEUX - Julien ROZIERE
E-mail : arthur.robieux@gmail.com, julien.roziere@u-psud.fr
Date : 24/04/2017

Objectif : Utilisation hiérarchique de nos fonctions afin de répondre au problème initial.
"""

from math import sqrt


def distanceResidus(res1, res2):
    """
    Fonction qui calcul la distance entre 2 résidus. Celle ci correspond à la distance
    minimale entre chacun des atomes pris 2 à 2 de chacun des ces résidus.

    :param res1: 1er résidu
    :param res2: 2eme résidu
    :return: Distance minimale entre les 2 atomes les plus proches de ces 2 résidus
    """

    distMin = 10000000

    for atom1 in res1["atomlist"]:
        coord1 = [res1[atom1]["x"], res1[atom1]["y"], res1[atom1]["z"]]
        for atom2 in res2["atomlist"]:
            coord2 = [res2[atom2]["x"], res2[atom2]["y"], res2[atom2]["z"]]

            dist = sqrt((coord2[0] - coord1[0]) ** 2 +
                        (coord2[1] - coord1[1]) ** 2 +
                        (coord2[2] - coord1[2]) ** 2)

            if dist < distMin:
                distMin = dist


    return distMin


