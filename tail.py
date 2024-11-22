"""
Group member:
- Leo Di Puma           - Louis Belboom
- Emile de Lamalle      - Diego Schyns
- Amos David            - Antoine Van Hoye
- Guillaume Vaillot
"""

import numpy as np
import matplotlib as plt
###Variables a changer car pas sur des données#############

V_h = 0.74 #m^3 horizontal tail volume
V_v=0.071 #m^3 vertical tail volume
S = 164.9 #m² area of the wing 
MAC = 4.622 #m mean chord of the line 
b = 29 #m span

#NACA 0012 for vertical parts of the tail
def tail(x_h,x_v):

    #calcul des surfaces d une tail conventionnel 
    S_h=V_h*S*MAC/(x_h)
    S_v=V_v*S*b/x_v
    S_tot_tail=np.sqrt(S_v**2+S_h**2)
    #calcul du butterfly angle
    gamma_h=np.arctan(S_v/S_h)
    #Calcul de l AR de la tail

    b_tail=5 #m longueur de la queue
    AR_tail=b_tail**2/S_tot_tail
    
    taper_ratio_tail=0.5
    c_tail_root=2.5 # m chord a la racine de la tail horizontal 
    c_tail_tip=1.5 #m chord a la fin de la tail horizontal 
    sweep_angle_tail=0.726632927 #radians 
    c_v_tail=2*S/(b(1+taper_ratio_tail)) #tail volume coefficient

    return gamma_h
