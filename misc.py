"""
Author : José CUNHA TEIXEIRA
License : SNCF Réseau, UMR 7619 METIS
Date : November 30, 2023
"""

from math import sqrt



def verify_expected(kwargs, list):
    for key in kwargs:
        if key not in list:
            diag_print("ERROR", "", "Argument {} not expected".format(key))
            raise SystemExit
        


def distance(pos1, pos2):
    if len(pos1) != len(pos2) or len(pos1) != 2:
        diag_print("ERROR", "distance", "Size of both elements is not (2,)")
    x1 = pos1[0]
    y1 = pos1[1]
    x2 = pos2[0]
    y2 = pos2[1]
    return sqrt( (x2-x1)**2 + (y2-y1)**2 )
