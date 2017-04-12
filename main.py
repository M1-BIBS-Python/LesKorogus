#!/usr/bin/env  python
# -*- coding:utf8 -*-

"""
Auteur : Arthur ROBIEUX
E-mail : arthur.robieux@gmail.com
Date : 12/04/2017

Objectif : Utilisation hiérarchique de nos fonctions afin de répondre au problème initial.
"""

from parsePDB import *


########################################
#Parsage de la conformation de référence
########################################

parse_ref = parsePDBfile("pab21_ref.pdb")

f = open("parse_ref.txt", "w")
f.write(str(parse_ref))
f.close()

##############################################
#Parsage des 500 conformations de la structure
##############################################

parse_frames = parsePDBmulti("pab21_500frames.pdb")

f = open("parse_frames.txt", "w")
f.write(str(parse_frames))
f.close()


###########################
#Vérifications des parsages
###########################

#len(parse_ref[" "].keys()) = 10 000 résidus
#len(parse_frames["5000"][" "].keys()) = 321 résidus

#Pourquoi la protéine de ref a 10 000 clés, et à un temps T que 321 clés ? (clés = residue seq number)
#Vérification : Les 321 clés aux temps T sont dans la protéine de référence (elle est plus complète).

commun = 0

for a in parse_frames["5000"][" "].keys():
    if a in parse_ref[" "].keys():
        commun += 1

print commun



###############################################
#Calcul des RMSD pour chacune des conformations
###############################################


