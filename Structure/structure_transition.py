# -*- coding: utf-8 -*-
"""
Created on Wed May 14 16:26:24 2025

@author: Emile
"""

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import math


# Geometry wings
beta = np.radians(0) # Setting angle [rad] 
(dx_ac,dy_ac,dz_ac) = ((7.788-4.5), (8.977-5.363), 0) # Coord AC of the regular tapered + transition wing 
(dx_w_cg,dy_w_cg,dz_w_cg) = (4.5969564, 4.7624667, -0.0069718) # Placement of the center of gravity of the wing to the root of the regular wing


# ---- Structural loads ---- (will be imported from another fct, here just expl values)

# Here, import all the n and all the alpha from the envelope : 
n =  [2.5, 2.5, 0, -1, -1] # Load factor
alpha =  [0.06151096, -0.04670625, -0.10033783, -0.13220202, -0.24113608]  # [rad]
L_overall = [1057331.80668917, 1095922.26154146, 79107.31234472, -336009.15840703, -347684.38494123]  # [N]
D_wing =   [5313.09423855, 15816.89997096, 15816.89997096, 10122.81598142, 2123.01255461]   # [N]
M_wing =  [-303561.349929, -903691.7644038, -903691.7644038, -578362.72921843, -121297.40751009]  # [N.m]

Surf_aile_equiv = 119.85 #[m²]
Surf_demi_wing_studied = 21.86460271713637 # [m²]
Surf_demi_trans_wing = 13.78039728282926 # [m²]
Coeff = (2*(Surf_demi_trans_wing + Surf_demi_trans_wing))/Surf_aile_equiv

Lwt = np.array(L_overall)*Coeff # Extract the correct lift
Dwt = np.array(D_wing)          # Will be divided by 2 in the structural loads
Wwt = (4735.9488633)*9.81       # [kg] mass -> to force [N] 
Mw = np.array(M_wing)

# ---- Material ----

# CFRP
# sigma_x_0 = 1500*10**6
# tau_max = 94*10**6
# safety_factor = 1.5

# Alu moyen Granta
sigma_x_0 = (797+683)/2*10**6  
tau_max = (402+345)/2*10**6
safety_factor = 1.5

# print('')
# print('Lift transition wing = ', Lwt/2)
# print('')
# print('Drag transition wing =', Dwt/2)
# print('')
# print('Moment transition wing =', Mw/2)
# print('')
# print('Weight transition + regular tapered wing (N) =', Wwt/2)
# print('')


# ---- Structural loads (function) ----
def structural_loads_regular_wing (n, alpha, Lwt, Dwt, Wwt, Mw):
    T_x = 0
    T_y = ((n*Wwt/2-Lwt/2)*np.sin(alpha+beta) + Dwt/2*np.cos(alpha+beta))
    T_z = ((-n*Wwt/2+Lwt/2)*np.cos(alpha+beta) + Dwt/2*np.sin(alpha+beta))
    
    M_x = np.cos(alpha+beta)/2 * (-n*Wwt* dy_w_cg + Lwt * dy_ac - Dwt * dz_ac) + np.sin(alpha+beta)/2 * (-n*Wwt * dz_w_cg + Lwt * dz_ac + Dwt * dy_ac) - Mw/2
    M_y = np.cos(alpha+beta)/2 * (n*Wwt * dx_w_cg - Lwt * dx_ac) + np.sin(alpha+beta)/2 * (Dwt * dx_ac) 
    M_z = np.sin(alpha+beta)/2 * (n*Wwt * dx_w_cg - Lwt * dx_ac) + np.cos(alpha+beta)/2 * (Dwt * dx_ac)  
    return (T_x,T_y,T_z,M_x,M_y,M_z)

# -------------------------------------------------------------------

def dist_2_pts(x1, y1, x2, y2):
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

# -------------------------------------------------------------------
def dist_2_booms(x_coord, y_coord): 
    dist = [] 
    max_dist = 0
    max_pair = (None, None)  # Will hold the points of the max distance

    n = len(x_coord)
    for i in range(1, n):
        dx = x_coord[i] - x_coord[i-1]
        dy = y_coord[i] - y_coord[i-1]
        d = np.sqrt(dx**2 + dy**2)
        dist.append(d)

        if d > max_dist:
            max_dist = d
            max_pair = ((x_coord[i-1], y_coord[i-1]), (x_coord[i], y_coord[i]))

    # Add distance between last and first point to close the shape
    dx = x_coord[0] - x_coord[-1]
    dy = y_coord[0] - y_coord[-1]
    d = np.sqrt(dx**2 + dy**2)
    dist.append(d)

    if d > max_dist:
        max_dist = d
        max_pair = ((x_coord[-1], y_coord[-1]), (x_coord[0], y_coord[0]))

    return dist, max_dist, max_pair

# -------------------------------------------------------------------
def boom_area(z_booms_ordered_centroid, y_booms_ordered_centroid, M_y, M_z, sigma_x_0, safety_factor) : 
    
    # Inertia per unit area (the boom area 'B' is considered to be the same everywhere)
    I_yy_over_B = np.sum(np.array(z_booms_ordered_centroid)**2) 
    I_zz_over_B = np.sum(np.array(y_booms_ordered_centroid)**2) 
    I_yz_over_B = np.sum(np.array(y_booms_ordered_centroid) * np.array(z_booms_ordered_centroid)) 
    
    # Direct stress in each boom (check formula)
    denom = (I_yy_over_B * I_zz_over_B - I_yz_over_B**2)  # Denominator of the stress equation
    # Compute stress for each boom
    B_sigma_xx = []
    for y, z in zip(y_booms_ordered_centroid, z_booms_ordered_centroid):
        stress = ((I_zz_over_B * M_y + I_yz_over_B * M_z) * z - (I_yz_over_B * M_y + I_yy_over_B * M_z) * y) / denom
        B_sigma_xx.append(stress)
    
    # Find maximum stress
    max_B_stress = max(B_sigma_xx)  # Find maximum stress value
    max_B_stress_index = B_sigma_xx.index(max_B_stress)  # Find the index of the maximum stress
    
    # ---- Minimum Area -----
    sigma_xx_max_material = sigma_x_0 / safety_factor
    B_min = max_B_stress/sigma_xx_max_material
    
    # Real stress in the boom with the new value of B 
    sigma_xx = np.array(B_sigma_xx) / B_min
    
    return B_min, sigma_xx


# -------------------------------------------------------------------

def polygon_area(x_coords, y_coords): # Fct to compute the area of the cells
    n = len(x_coords)
    area = 0.0

    for i in range(n):
        j = (i + 1) % n  # Wraps around to connect last to first
        area += x_coords[i] * y_coords[j]
        area -= x_coords[j] * y_coords[i]

    return abs(area) / 2.0

# -------------------------------------------------------------------

def swept_area_from_center(y_points, z_points, y_centroid, z_centroid):
    n = len(y_points)
    total_area = 0.0
    individual_areas = []

    for i in range(n):
        y1, z1 = y_points[i], z_points[i]
        y2, z2 = y_points[(i + 1) % n], z_points[(i + 1) % n]  # wrap around

        # Area computation
        area = 0.5 * abs(y1 * z2 + y2 * z_centroid + y_centroid * z1 - z1 * y2 - z2 * y_centroid - z_centroid * y1)
        individual_areas.append(area)
        total_area += area

    return total_area, individual_areas



# -------------------------------------------------------------------
def skin_thickness(B, sigma_xx, delta_x, delta_y, delta_z, T_y, T_z, M_x, y_booms_ordered_centroid, z_booms_ordered_centroid, y_booms_cell_1, z_booms_cell_1, tau_max, y_booms_cell_2, z_booms_cell_2, y_centroid, z_centroid, y_booms_cell_3, z_booms_cell_3, y_booms_cell_4, z_booms_cell_4) :
    
    
    # Taper effect (suite) :
    T_y_web = T_y - B* np.sum(sigma_xx * (delta_y/delta_x))
    T_z_web = T_z - B* np.sum(sigma_xx * (delta_z/delta_x))
    print('T_y', T_y)
    print('T_y_web', T_y_web)    
    print('T_z', T_z)
    print('T_z_web', T_z_web)  
    
    # Inertia per unit area (the boom area 'B' is considered to be the same everywhere)
    I_yy_over_B = np.sum(np.array(z_booms_ordered_centroid)**2) 
    I_zz_over_B = np.sum(np.array(y_booms_ordered_centroid)**2) 
    I_yz_over_B = np.sum(np.array(y_booms_ordered_centroid) * np.array(z_booms_ordered_centroid)) 
    I_yy = I_yy_over_B * B
    I_zz = I_zz_over_B * B
    I_yz = I_yz_over_B * B
    
    # ---- Open shear flow ----
    denom = (I_yy*I_zz - I_yz**2)
    factor_1 = (I_zz*T_z_web - I_yz*T_y_web)/denom
    factor_2 = (I_yy*T_y_web - I_yz*T_z_web)/denom
    # print('factor_1 :', factor_1)
    # print('factor_2 :', factor_2)
    
    
    # -- Cell 1 -- (anticlockwise : q > 0)
    q_0_cell_1 = [0] # Cut in cell 1, the open shear flow between booms 7 and 6 is zero                                              (q_0_7,6) (Cut)
    q_0_cell_1.append(q_0_cell_1[0] - factor_1 * z_booms_ordered_centroid[6-1]*B - factor_2 * y_booms_ordered_centroid[6-1]*B)     # (q_0_6,5)
    q_0_cell_1.append(q_0_cell_1[1] - factor_1 * z_booms_ordered_centroid[5-1]*B - factor_2 * y_booms_ordered_centroid[5-1]*B)     # (q_0_5,4)
    q_0_cell_1.append(q_0_cell_1[2] - factor_1 * z_booms_ordered_centroid[4-1]*B - factor_2 * y_booms_ordered_centroid[4-1]*B)     # (q_0_4,3)
    q_0_cell_1.append(q_0_cell_1[3] - factor_1 * z_booms_ordered_centroid[3-1]*B - factor_2 * y_booms_ordered_centroid[3-1]*B)     # (q_0_3,2)
    q_0_cell_1.append(q_0_cell_1[4] - factor_1 * z_booms_ordered_centroid[2-1]*B - factor_2 * y_booms_ordered_centroid[2-1]*B)     # (q_0_2,1)
    q_0_cell_1.append(q_0_cell_1[5] - factor_1 * z_booms_ordered_centroid[1-1]*B - factor_2 * y_booms_ordered_centroid[1-1]*B)     # (q_0_1,72)
    q_0_cell_1.append(q_0_cell_1[6] - factor_1 * z_booms_ordered_centroid[72-1]*B - factor_2 * y_booms_ordered_centroid[72-1]*B)   # (q_0_72,71)
    q_0_cell_1.append(q_0_cell_1[7] - factor_1 * z_booms_ordered_centroid[71-1]*B - factor_2 * y_booms_ordered_centroid[71-1]*B)   # (q_0_71,70)
    q_0_cell_1.append(q_0_cell_1[8] - factor_1 * z_booms_ordered_centroid[70-1]*B - factor_2 * y_booms_ordered_centroid[70-1]*B)   # (q_0_70,69)
    q_0_cell_1.append(q_0_cell_1[9] - factor_1 * z_booms_ordered_centroid[69-1]*B - factor_2 * y_booms_ordered_centroid[69-1]*B)   # (q_0_69,68)
    q_0_cell_1.append(q_0_cell_1[10] - factor_1 * z_booms_ordered_centroid[68-1]*B - factor_2 * y_booms_ordered_centroid[68-1]*B)  # (q_0_68,67) 
    q_0_cell_1.append(q_0_cell_1[11] - factor_1 * z_booms_ordered_centroid[67-1]*B - factor_2 * y_booms_ordered_centroid[67-1]*B)  # (q_0_67,66) 
    q_spar_1_2 = q_0_cell_1[12] - factor_1 * z_booms_ordered_centroid[66-1]*B - factor_2 * y_booms_ordered_centroid[66-1]*B        # (q_0_66,7)     => q_spar_1_2 !!
    q_0_cell_1.append(q_spar_1_2)                                                                                                  # (q_0_66,7)     => q_spar_1_2 !!
    
    # Remettre en ordre 
    # Étape 1 : trouver l'indice du point de départ : max(y), puis max(z) parmi les égaux
    idx_start_candidates = np.where(y_booms_cell_1 == np.max(y_booms_cell_1))[0]
    if len(idx_start_candidates) > 1:
        idx_start = idx_start_candidates[np.argmax(z_booms_cell_1[idx_start_candidates])]
    else:
        idx_start = idx_start_candidates[0]
    
    # Étape 2 : trier tous les points selon z décroissant
    sorted_indices = np.argsort(-z_booms_cell_1)  # z décroissant
    
    # Étape 3 : faire tourner le tableau pour que le point de départ soit en premier
    start_in_sorted = np.where(sorted_indices == idx_start)[0][0]
    sorted_indices = np.roll(sorted_indices, -start_in_sorted)
    
    # Points réordonnés
    y_booms_cell_1_sorted = y_booms_cell_1[sorted_indices]
    z_booms_cell_1_sorted = z_booms_cell_1[sorted_indices]
    # print("y sorted :", y_booms_cell_1_sorted)
    # print("z sorted :", z_booms_cell_1_sorted)
    
    
    # -- Cell 2 -- 
    
    # - Upper - (clockwise : q < 0)
    q_0_cell_2_upper = [q_spar_1_2]                                                                                                           # (q_0_66,7)     => q_spar_1_2 !!
    q_0_cell_2_upper.append(q_0_cell_2_upper[0] - factor_1 * z_booms_ordered_centroid[7-1]*B - factor_2 * y_booms_ordered_centroid[7-1]*B)    # (q_0_7,8)
    q_0_cell_2_upper.append(q_0_cell_2_upper[1] - factor_1 * z_booms_ordered_centroid[8-1]*B - factor_2 * y_booms_ordered_centroid[8-1]*B)    # (q_0_8,9)
    q_0_cell_2_upper.append(q_0_cell_2_upper[2] - factor_1 * z_booms_ordered_centroid[9-1]*B - factor_2 * y_booms_ordered_centroid[9-1]*B)    # (q_0_9,10)
    q_0_cell_2_upper.append(q_0_cell_2_upper[3] - factor_1 * z_booms_ordered_centroid[10-1]*B - factor_2 * y_booms_ordered_centroid[10-1]*B)  # (q_0_10,11)
    q_0_cell_2_upper.append(q_0_cell_2_upper[4] - factor_1 * z_booms_ordered_centroid[11-1]*B - factor_2 * y_booms_ordered_centroid[11-1]*B)  # (q_0_11,12)
    q_0_cell_2_upper.append(q_0_cell_2_upper[5] - factor_1 * z_booms_ordered_centroid[12-1]*B - factor_2 * y_booms_ordered_centroid[12-1]*B)  # (q_0_12,13)
    q_0_cell_2_upper.append(q_0_cell_2_upper[6] - factor_1 * z_booms_ordered_centroid[13-1]*B - factor_2 * y_booms_ordered_centroid[13-1]*B)  # (q_0_13,14)
    q_0_cell_2_upper.append(q_0_cell_2_upper[7] - factor_1 * z_booms_ordered_centroid[14-1]*B - factor_2 * y_booms_ordered_centroid[14-1]*B)  # (q_0_14,15)
    q_0_cell_2_upper.append(q_0_cell_2_upper[8] - factor_1 * z_booms_ordered_centroid[15-1]*B - factor_2 * y_booms_ordered_centroid[15-1]*B)  # (q_0_15,16)
    q_0_cell_2_upper.append(q_0_cell_2_upper[9] - factor_1 * z_booms_ordered_centroid[16-1]*B - factor_2 * y_booms_ordered_centroid[16-1]*B)  # (q_0_16,17)
    q_0_cell_2_upper.append(q_0_cell_2_upper[10] - factor_1 * z_booms_ordered_centroid[17-1]*B - factor_2 * y_booms_ordered_centroid[17-1]*B) # (q_0_17,18)
    q_spar_2_3 = q_0_cell_2_upper[11] - factor_1 * z_booms_ordered_centroid[18-1]*B - factor_2 * y_booms_ordered_centroid[18-1]*B             # (q_0_18,55)  => q_spar_2_3 !!
    q_0_cell_2_upper.append(q_spar_2_3)                                                                                                       # (q_0_18,55)  => q_spar_2_3 !!
    q_0_18_55 = q_spar_2_3                                                                                                                    # (q_0_18,55)  (need to know the seperate value)

    q_0_cell_2_upper_flipped = -np.array(q_0_cell_2_upper)         # to put anticlockwise   (invert all the values)                                                                                                                # (q_0_18,55)  (need to know the seperate value)
    q_0_cell_2_upper = q_0_cell_2_upper_flipped[::-1]    # invert the order
    
    # - Lower - (anticlockwise : q > 0)
    q_0_cell_2_lower = [0]                                                                                                                    # (q_0_66,65) (Cut)
    q_0_cell_2_lower.append(q_0_cell_2_lower[0] - factor_1 * z_booms_ordered_centroid[65-1]*B - factor_2 * y_booms_ordered_centroid[65-1]*B)  # (q_0_65,64)
    q_0_cell_2_lower.append(q_0_cell_2_lower[1] - factor_1 * z_booms_ordered_centroid[64-1]*B - factor_2 * y_booms_ordered_centroid[64-1]*B)  # (q_0_64,63)
    q_0_cell_2_lower.append(q_0_cell_2_lower[2] - factor_1 * z_booms_ordered_centroid[63-1]*B - factor_2 * y_booms_ordered_centroid[63-1]*B)  # (q_0_63,62)
    q_0_cell_2_lower.append(q_0_cell_2_lower[3] - factor_1 * z_booms_ordered_centroid[62-1]*B - factor_2 * y_booms_ordered_centroid[62-1]*B)  # (q_0_62,61)
    q_0_cell_2_lower.append(q_0_cell_2_lower[4] - factor_1 * z_booms_ordered_centroid[61-1]*B - factor_2 * y_booms_ordered_centroid[61-1]*B)  # (q_0_61,60)
    q_0_cell_2_lower.append(q_0_cell_2_lower[5] - factor_1 * z_booms_ordered_centroid[60-1]*B - factor_2 * y_booms_ordered_centroid[60-1]*B)  # (q_0_60,59)
    q_0_cell_2_lower.append(q_0_cell_2_lower[6] - factor_1 * z_booms_ordered_centroid[59-1]*B - factor_2 * y_booms_ordered_centroid[59-1]*B)  # (q_0_59,58)
    q_0_cell_2_lower.append(q_0_cell_2_lower[7] - factor_1 * z_booms_ordered_centroid[58-1]*B - factor_2 * y_booms_ordered_centroid[58-1]*B)  # (q_0_58,57)
    q_0_cell_2_lower.append(q_0_cell_2_lower[8] - factor_1 * z_booms_ordered_centroid[57-1]*B - factor_2 * y_booms_ordered_centroid[57-1]*B)  # (q_0_57,56)
    q_0_56_55 = q_0_cell_2_lower[9] - factor_1 * z_booms_ordered_centroid[56-1]*B - factor_2 * y_booms_ordered_centroid[56-1]*B               # (q_0_56,55) (need to know the seperate value)
    q_0_cell_2_lower.append(q_0_56_55)                                                                                                        # (q_0_56,55)
      
    q_0_cell_2 = np.concatenate((q_0_cell_2_lower, q_0_cell_2_upper)) # assemble the openshear flow for cell 2
    y_booms_cell_2_sorted = y_booms_cell_2[::-1]  # reorder the coordinates of the booms to match the order of the openshear flow
    z_booms_cell_2_sorted = z_booms_cell_2[::-1] 
  
    
    # -- Cell 3 --
    
    # - Upper - (clockwise : q < 0)
    q_0_cell_3_upper = [0]                                                                                                                    # (q_0_18,19) (Cut)
    q_0_cell_3_upper.append(q_0_cell_3_upper[0] - factor_1 * z_booms_ordered_centroid[19-1]*B - factor_2 * y_booms_ordered_centroid[19-1]*B)  # (q_0_19,20)
    q_0_cell_3_upper.append(q_0_cell_3_upper[1] - factor_1 * z_booms_ordered_centroid[20-1]*B - factor_2 * y_booms_ordered_centroid[20-1]*B)  # (q_0_20,21)
    q_0_cell_3_upper.append(q_0_cell_3_upper[2] - factor_1 * z_booms_ordered_centroid[21-1]*B - factor_2 * y_booms_ordered_centroid[21-1]*B)  # (q_0_21,22)
    q_0_cell_3_upper.append(q_0_cell_3_upper[3] - factor_1 * z_booms_ordered_centroid[22-1]*B - factor_2 * y_booms_ordered_centroid[22-1]*B)  # (q_0_22,23)
    q_0_cell_3_upper.append(q_0_cell_3_upper[4] - factor_1 * z_booms_ordered_centroid[23-1]*B - factor_2 * y_booms_ordered_centroid[23-1]*B)  # (q_0_23,24)
    
    q_0_cell_3_upper = (-np.array(q_0_cell_3_upper))[::-1] # invert the signs (to be anticlokwise), then invert all the positions of all the elements (last becomes first, etc)
    
    # - Lower - (anticlockwise : q > 0)
    q_0_cell_3_lower = [q_0_18_55]                                                                                                                # (q_0_18,55) => q_spar_2_3 !!
    q_0_cell_3_lower.append(q_0_56_55 +  q_0_18_55   - factor_1 * z_booms_ordered_centroid[55-1]*B - factor_2 * y_booms_ordered_centroid[55-1]*B) # (q_0_55,54) = q_0_56,55 + q_0_18,55 + contribution of the boom 55
    q_0_cell_3_lower.append(q_0_cell_3_lower[1] - factor_1 * z_booms_ordered_centroid[54-1]*B - factor_2 * y_booms_ordered_centroid[54-1]*B)      # (q_0_54,53)
    q_0_cell_3_lower.append(q_0_cell_3_lower[2] - factor_1 * z_booms_ordered_centroid[53-1]*B - factor_2 * y_booms_ordered_centroid[53-1]*B)      # (q_0_53,52)
    q_0_cell_3_lower.append(q_0_cell_3_lower[3] - factor_1 * z_booms_ordered_centroid[52-1]*B - factor_2 * y_booms_ordered_centroid[52-1]*B)      # (q_0_52,51)
    q_0_cell_3_lower.append(q_0_cell_3_lower[4] - factor_1 * z_booms_ordered_centroid[51-1]*B - factor_2 * y_booms_ordered_centroid[51-1]*B)      # (q_0_51,50)
    q_0_cell_3_lower.append(q_0_cell_3_lower[5] - factor_1 * z_booms_ordered_centroid[50-1]*B - factor_2 * y_booms_ordered_centroid[50-1]*B)      # (q_0_50,49)
    q_spar_3_4 = q_0_cell_3_lower[6] - factor_1 * z_booms_ordered_centroid[49-1]*B - factor_2 * y_booms_ordered_centroid[49-1]*B                  # (q_0_49,24) => q_spar_3_4 !!
    q_0_cell_3_lower.append(q_spar_3_4)                                                                                                           # (q_0_49,24) => q_spar_3_4 !!
    q_0_49_24 = q_spar_3_4                                                                                                                        # (q_0_49,24) (need to know the seperate value)

    q_0_cell_3 = np.concatenate((q_0_cell_3_lower, q_0_cell_3_upper)) # assemble the openshear flow for cell 3
    q_0_cell_3 = np.concatenate((q_0_cell_3[1:], q_0_cell_3[:1])) # Décaler le premier élément (q_spar_2,3) à la fin du tableau
    y_booms_cell_3_sorted = y_booms_cell_3[::-1]  # reorder the coordinates of the booms to match the order of the openshear flow
    z_booms_cell_3_sorted = z_booms_cell_3[::-1] 
    
    
    # -- Cell 4 -- (anticlockwise : q > 0)
    
    q_0_cell_4 = [0]                                                                                                               # (q_0_49,48) (Cut)
    q_0_cell_4.append(q_0_cell_4[0] - factor_1 * z_booms_ordered_centroid[48-1]*B - factor_2 * y_booms_ordered_centroid[48-1]*B)   # (q_0_48,47)
    q_0_cell_4.append(q_0_cell_4[1] - factor_1 * z_booms_ordered_centroid[47-1]*B - factor_2 * y_booms_ordered_centroid[47-1]*B)   # (q_0_47,46)
    q_0_cell_4.append(q_0_cell_4[2] - factor_1 * z_booms_ordered_centroid[46-1]*B - factor_2 * y_booms_ordered_centroid[46-1]*B)   # (q_0_46,45)
    q_0_cell_4.append(q_0_cell_4[3] - factor_1 * z_booms_ordered_centroid[45-1]*B - factor_2 * y_booms_ordered_centroid[45-1]*B)   # (q_0_45,44)
    q_0_cell_4.append(q_0_cell_4[4] - factor_1 * z_booms_ordered_centroid[44-1]*B - factor_2 * y_booms_ordered_centroid[44-1]*B)   # (q_0_44,43)
    q_0_cell_4.append(q_0_cell_4[5] - factor_1 * z_booms_ordered_centroid[43-1]*B - factor_2 * y_booms_ordered_centroid[43-1]*B)   # (q_0_43,42)
    q_0_cell_4.append(q_0_cell_4[6] - factor_1 * z_booms_ordered_centroid[42-1]*B - factor_2 * y_booms_ordered_centroid[42-1]*B)   # (q_0_42,41)
    q_0_cell_4.append(q_0_cell_4[7] - factor_1 * z_booms_ordered_centroid[41-1]*B - factor_2 * y_booms_ordered_centroid[41-1]*B)   # (q_0_41,40)
    q_0_cell_4.append(q_0_cell_4[8] - factor_1 * z_booms_ordered_centroid[40-1]*B - factor_2 * y_booms_ordered_centroid[40-1]*B)   # (q_0_40,39)
    q_0_cell_4.append(q_0_cell_4[9] - factor_1 * z_booms_ordered_centroid[39-1]*B - factor_2 * y_booms_ordered_centroid[39-1]*B)   # (q_0_39,38)
    q_0_cell_4.append(q_0_cell_4[10] - factor_1 * z_booms_ordered_centroid[38-1]*B - factor_2 * y_booms_ordered_centroid[38-1]*B)  # (q_0_38,37) 
    q_0_cell_4.append(q_0_cell_4[11] - factor_1 * z_booms_ordered_centroid[37-1]*B - factor_2 * y_booms_ordered_centroid[37-1]*B)  # (q_0_37,36) 
    q_0_cell_4.append(q_0_cell_4[12] - factor_1 * z_booms_ordered_centroid[36-1]*B - factor_2 * y_booms_ordered_centroid[36-1]*B)  # (q_0_36,35)
    q_0_cell_4.append(q_0_cell_4[13] - factor_1 * z_booms_ordered_centroid[35-1]*B - factor_2 * y_booms_ordered_centroid[35-1]*B)  # (q_0_35,34)
    q_0_cell_4.append(q_0_cell_4[14] - factor_1 * z_booms_ordered_centroid[34-1]*B - factor_2 * y_booms_ordered_centroid[34-1]*B)  # (q_0_34,33)
    q_0_cell_4.append(q_0_cell_4[15] - factor_1 * z_booms_ordered_centroid[33-1]*B - factor_2 * y_booms_ordered_centroid[33-1]*B)  # (q_0_33,32)
    q_0_cell_4.append(q_0_cell_4[16] - factor_1 * z_booms_ordered_centroid[32-1]*B - factor_2 * y_booms_ordered_centroid[32-1]*B)  # (q_0_32,31)
    q_0_cell_4.append(q_0_cell_4[17] - factor_1 * z_booms_ordered_centroid[31-1]*B - factor_2 * y_booms_ordered_centroid[31-1]*B)  # (q_0_31,30)
    q_0_cell_4.append(q_0_cell_4[18] - factor_1 * z_booms_ordered_centroid[30-1]*B - factor_2 * y_booms_ordered_centroid[30-1]*B)  # (q_0_30,29)
    q_0_cell_4.append(q_0_cell_4[19] - factor_1 * z_booms_ordered_centroid[29-1]*B - factor_2 * y_booms_ordered_centroid[29-1]*B)  # (q_0_29,28)
    q_0_cell_4.append(q_0_cell_4[20] - factor_1 * z_booms_ordered_centroid[28-1]*B - factor_2 * y_booms_ordered_centroid[28-1]*B)  # (q_0_28,27)
    q_0_cell_4.append(q_0_cell_4[21] - factor_1 * z_booms_ordered_centroid[27-1]*B - factor_2 * y_booms_ordered_centroid[27-1]*B)  # (q_0_27,26)
    q_0_cell_4.append(q_0_cell_4[22] - factor_1 * z_booms_ordered_centroid[26-1]*B - factor_2 * y_booms_ordered_centroid[26-1]*B)  # (q_0_26,25) 
    q_0_cell_4.append(q_0_cell_4[23] - factor_1 * z_booms_ordered_centroid[25-1]*B - factor_2 * y_booms_ordered_centroid[25-1]*B)  # (q_0_25,24) 
    q_0_cell_4.append(-q_0_49_24)                                                                                                  # (q_0_24,49) => -q_spar_3_4 !!
    
    y_booms_cell_4_sorted = y_booms_cell_4[::-1] # reorder the coordinates of the booms to match the order of the openshear flow
    z_booms_cell_4_sorted = z_booms_cell_4[::-1]
    
    
    
    # ---- Lengths of the spars and the cells ----
    # - Spars -
    l_spar_1_2 = dist_2_pts(y_booms_ordered_centroid[7-1], z_booms_ordered_centroid[7-1], y_booms_ordered_centroid[66-1], z_booms_ordered_centroid[66-1])
    l_spar_2_3 = dist_2_pts(y_booms_ordered_centroid[18-1], z_booms_ordered_centroid[18-1], y_booms_ordered_centroid[55-1], z_booms_ordered_centroid[55-1])
    l_spar_3_4 =  dist_2_pts(y_booms_ordered_centroid[24-1], z_booms_ordered_centroid[24-1], y_booms_ordered_centroid[49-1], z_booms_ordered_centroid[49-1])
    # print('l_spar_1_2 =', l_spar_1_2)
    # print('l_spar_2_3 =', l_spar_2_3)
    # print('l_spar_3_4 =', l_spar_3_4)

    # - Cell 1 -
    lengths_c_1 = dist_2_booms(y_booms_cell_1_sorted, z_booms_cell_1_sorted)[0]
    l_cell_1 = np.sum(lengths_c_1)
    #print('lengths_c_1 = ', lengths_c_1)
    #print(l_cell_1)
    
    # - Cell 2 -
    lengths_c_2 = dist_2_booms(y_booms_cell_2_sorted, z_booms_cell_2_sorted)[0]
    l_cell_2 = np.sum(lengths_c_2)
    #print('lengths_c_2 = ', lengths_c_2)
    #print(l_cell_2)
    
    # - Cell 3 -
    lengths_c_3 = dist_2_booms(y_booms_cell_3_sorted, z_booms_cell_3_sorted)[0]
    l_cell_3 = np.sum(lengths_c_3)
    # print('lengths_c_3 = ', lengths_c_3)
    # print(l_cell_3)
    
    # - Cell 4 -
    lengths_c_4 = dist_2_booms(y_booms_cell_4_sorted, z_booms_cell_4_sorted)[0]
    l_cell_4 = np.sum(lengths_c_4)
    # print('lengths_c_4 = ', lengths_c_4)
    # print(l_cell_4)
    
    
    # ---- Cell Area ----
    A_h_c_1 = polygon_area(y_booms_cell_1_sorted, z_booms_cell_1_sorted)
    A_h_c_2 = polygon_area(y_booms_cell_2_sorted, z_booms_cell_2_sorted)
    A_h_c_3 = polygon_area(y_booms_cell_3_sorted, z_booms_cell_3_sorted)
    A_h_c_4 = polygon_area(y_booms_cell_4_sorted, z_booms_cell_4_sorted)
    #print(f"Cell Areas: A_h_c_1 = {A_h_c_1:.4f} m², A_h_c_2 = {A_h_c_2:.4f} m², A_h_c_3 = {A_h_c_3:.4f} m², A_h_c_4 = {A_h_c_4:.4f} m².")
    

    # ---- Open shear flow integration ----
    # - Cell 1 -
    q_0_cell_1 = np.array(q_0_cell_1) # Convert it to NumPy arrays
    lengths_c_1 = np.array(lengths_c_1)
    # print('q_0_cell_1 =', q_0_cell_1)
    # print('lengths_c_1 =', lengths_c_1)
    open_shear_flux_1 = np.sum(q_0_cell_1 * lengths_c_1)
    # print('open_shear_flux_1 =', open_shear_flux_1)
    
    # - Cell 2 -
    q_0_cell_2 = np.array(q_0_cell_2) # Convert the lists to NumPy arrays
    lengths_c_2 = np.array(lengths_c_2)
    # print('q_0_cell_2 =', q_0_cell_2)
    # print('lengths_c_2 =', lengths_c_2)
    open_shear_flux_2 = np.sum(q_0_cell_2 * lengths_c_2)
    # print('open_shear_flux_2 =', open_shear_flux_2)
    
    # - Cell 3 -
    q_0_cell_3 = np.array(q_0_cell_3) # Convert it to NumPy arrays
    lengths_c_3 = np.array(lengths_c_3)
    # print('q_0_cell_3 =', q_0_cell_3)
    # print('lengths_c_3 =', lengths_c_3)
    open_shear_flux_3 = np.sum(q_0_cell_3 * lengths_c_3)
    # print('open_shear_flux_3 =', open_shear_flux_3)
    
    # - Cell 4 -
    q_0_cell_4 = np.array(q_0_cell_4) # Convert the lists to NumPy arrays
    lengths_c_4 = np.array(lengths_c_4)
    # print('q_0_cell_4 =', q_0_cell_4)
    # print('lengths_c_4 =', lengths_c_4)
    open_shear_flux_4 = np.sum(q_0_cell_4 * lengths_c_4)
    # print('open_shear_flux_4 =', open_shear_flux_4)
    
    
    # ---- Swept Areas ----
    # Cell 1
    swept_area_cell_1 = swept_area_from_center(y_booms_cell_1_sorted, z_booms_cell_1_sorted, y_centroid = 0, z_centroid = 0)[1] # it is an array
    # Cell 2
    swept_area_cell_2 = swept_area_from_center(y_booms_cell_2_sorted, z_booms_cell_2_sorted, y_centroid = 0, z_centroid = 0)[1]
    # Cell 3
    swept_area_cell_3 = swept_area_from_center(y_booms_cell_3_sorted, z_booms_cell_3_sorted, y_centroid = 0, z_centroid = 0)[1]
    # Cell 4
    swept_area_cell_4 = swept_area_from_center(y_booms_cell_4_sorted, z_booms_cell_4_sorted, y_centroid = 0, z_centroid = 0)[1]
    # Term in the momentum equilibrium (term swpet area est le terme ): 
    term_swept_area = 2 * (np.sum(q_0_cell_1 * swept_area_cell_1) + np.sum(q_0_cell_2 * swept_area_cell_2) + np.sum(q_0_cell_3 * swept_area_cell_3) + np.sum(q_0_cell_4 * swept_area_cell_4))
    # print('term_swept_area =', term_swept_area)
    
    
    # ---- Correction term computation : solve a 5 eqns syst (because 4 cells and 5 unknowns with the twist rate) ----
    # X^T = [Inc, q_1_corr, q_2_corr, q_3_corr, q_4_corr]
    
    # A = matrix of the coefficents
    
    A = np.array([[A_h_c_1, -l_cell_1, l_spar_1_2, 0, 0],
                  [A_h_c_2, l_spar_1_2, -l_cell_2, l_spar_2_3, 0],
                  [A_h_c_3, 0, l_spar_2_3, -l_cell_3, l_spar_3_4],
                  [A_h_c_4, 0, 0, l_spar_3_4, -l_cell_4],
                  [0, 2*A_h_c_1, 2*A_h_c_2, 2*A_h_c_3, 2*A_h_c_4]])
    
    # B = matrix of the independent terms
    term_P_z = B * np.sum(y_booms_ordered_centroid * sigma_xx * delta_z/delta_x)
    term_P_y = B * np.sum(y_booms_ordered_centroid * sigma_xx * delta_y/delta_x)
    ter_ind_mom = - term_swept_area - term_P_z + term_P_y - (y_centroid * T_z - z_centroid * T_y)
    B = np.array([open_shear_flux_1, open_shear_flux_2, open_shear_flux_3, open_shear_flux_4, ter_ind_mom])
    
    # Solve syst : 
    sol = np.linalg.solve(A, B)
    Inc, q_1_corr, q_2_corr, q_3_corr, q_4_corr = sol
    print(f"q_1_corr = {q_1_corr:.4f}, q_2_corr = {q_2_corr:.4f}, q_3_corr = {q_3_corr:.4f}, q_4_corr = {q_4_corr:.4f}")
    
    # ---- Shear flow due to torsion ----
    
    # Solve syst :    
    A_torsion = np.array([[A_h_c_1, -l_cell_1, l_spar_1_2, 0, 0],
                  [A_h_c_2, l_spar_1_2, -l_cell_2, l_spar_2_3, 0],
                  [A_h_c_3, 0, l_spar_2_3, -l_cell_3, l_spar_3_4],
                  [A_h_c_4, 0, 0, l_spar_3_4, -l_cell_4],
                  [0, 2*A_h_c_1, 2*A_h_c_2, 2*A_h_c_3, 2*A_h_c_4]])
    
    B_torsion = np.array([0, 0, 0, 0, M_x])
    
    sol_torsion = np.linalg.solve(A_torsion, B_torsion)
    Unkn, q_1_torsion, q_2_torsion, q_3_torsion, q_4_torsion = sol_torsion
    print(f"q_1_torsion = {q_1_torsion:.4f}, q_2_torsion = {q_2_torsion:.4f}, q_3_torsion = {q_3_torsion:.4f}, q_4_torsion = {q_4_torsion:.4f}")
    
    # ---- Shear flow (closed) ----
    q_closed_cell_1 = np.array(q_0_cell_1) + q_1_corr + q_1_torsion # Apply the correction and the shear flow due to torsion
    q_closed_cell_2 = np.array(q_0_cell_2) + q_2_corr + q_2_torsion
    q_closed_cell_3 = np.array(q_0_cell_3) + q_3_corr + q_3_torsion 
    q_closed_cell_4 = np.array(q_0_cell_4) + q_4_corr + q_4_torsion 
    
    
    # Take into account the spars which has the correction added in both sense (localise their position in the different arrays)
    # -- Spar 1_2 -- (I decide to count "> 0" if go from 66 to 7)
    q_closed_spar_1_2 = q_closed_cell_1[14-1] - q_2_corr - q_2_torsion
    q_closed_cell_1[14-1] = q_closed_spar_1_2
    q_closed_cell_2[24-1] = -q_closed_spar_1_2
    
    # -- Spar 2_3 --
    q_closed_spar_2_3 = q_closed_cell_2[12-1] - q_3_corr - q_3_torsion
    q_closed_cell_2[12-1] = q_closed_spar_2_3
    q_closed_cell_3[14-1] = -q_closed_spar_2_3
    
    # -- Spar 3_4 --
    q_closed_spar_3_4 = q_closed_cell_3[7-1] - q_4_corr - q_4_torsion
    q_closed_cell_3[7-1] = q_closed_spar_3_4
    q_closed_cell_4[26-1] = -q_closed_spar_3_4
    
    q_closed = np.concatenate([q_closed_cell_1, q_closed_cell_2, q_closed_cell_3, q_closed_cell_4]) # Note that q_closed contains all the spar shear flows' twice but it is not a problem since their absolute value is the same 
    #print('q_closed = ',q_closed)
    # ---- Thickness computation ----   
    print("q_max =", np.max(np.abs(q_closed)) )
    print("")
    #print('tau_max', tau_max)
    thickness = np.max(np.abs(q_closed)) / (tau_max/safety_factor) 
    
    return thickness


def plotAirfoil(plot_airfoil, n_booms):
    
    # ---- Geometry ----
    transition_semi_span = 3.10921985816 #[m]
    x_tip = transition_semi_span
    x_root = 0
    x_arrondi = 3.10921985816 #
    chord_root = 11.437 #[m]
    chord_tip = 4.018 #[m]
    sweep_angle = np.radians(31.599) #[°]
    delta_x = (x_arrondi - x_root)
    y_dist_root_tip = delta_x *np.tan(sweep_angle)
    #print('delta_x : ', delta_x)
    #print('y_dist_root_tip : ', y_dist_root_tip)
    #print('')
    
    # Load airfoil data
    data_wing = pd.read_csv("NACA45118_XYZ.csv")
    y_c_all = data_wing.iloc[:, 0].values
    z_c_all = data_wing.iloc[:, 1].values  

    # Separate and sort upper and lower surfaces
    upper_surface = data_wing[data_wing.iloc[:, 1] >= 0].copy().sort_values(by=data_wing.columns[0])  # z >= 0
    lower_surface = data_wing[data_wing.iloc[:, 1] < 0].copy().sort_values(by=data_wing.columns[0])  # z < 0

    
    # Add two symmetric booms at the leading edge (x = 0)
    #z_le_offset = 0.012  # vertical distance from centerline
    y_c_booms = [0.003, 0.0025]
    z_c_booms = [0.018, -0.01]
    
    # Ensure even number of remaining booms for symmetric pairs
    if (n_booms - 2) % 2 != 0:
        raise ValueError("Number of booms must be even to allow symmetric pairs around leading edge.")
    
    # Generate y positions from 0 to y_c_max, skipping the leading edge (already added)
    n_pairs = (n_booms - 2) // 2
    y_c_positions = np.linspace(0, 0.6, n_pairs + 1)[1:]

    # Interpolate z from upper and lower surfaces (find the z corresponding to the y selected)
    z_c_upper = np.interp(y_c_positions, upper_surface.iloc[:, 0], upper_surface.iloc[:, 1])
    z_c_lower = np.interp(y_c_positions, lower_surface.iloc[:, 0], lower_surface.iloc[:, 1])

    
    for y, z_u, z_l in zip(y_c_positions, z_c_upper, z_c_lower):
        z_center = 0.5 * (z_u + z_l)  # local "midline"
        z_offset = abs(z_u - z_center)
        y_c_booms.extend([y, y])
        z_c_booms.extend([z_center + z_offset, z_center - z_offset])
    

    # ---- Centroid location ----
    y_c_centroid = np.sum(y_c_booms) / len(y_c_booms)
    z_c_centroid = np.sum(z_c_booms) / len(z_c_booms)
    # Multiply by the chord
    y_booms = np.array(y_c_booms) * chord_root
    z_booms = np.array(z_c_booms) * chord_root
    y_centroid = y_c_centroid * chord_root
    z_centroid = z_c_centroid * chord_root
    print('')
    #print('Centroid location:')
    #print(f'y = {y_centroid} m')
    #print(f'z = {z_centroid} m')
    # Conversion to feet
    print(f'y_centroid = {y_centroid*3.28084} ft')
    print(f'z_centroid = {z_centroid*3.28084} ft')
    
    
    # ---- Verify that the distance between 2 succesive booms is between 0.1 and 0.2 m (Limits suggested by M.Noels) ----
    # Separate booms on upper and lower surfaces (skipping the first LE point)
    upper_booms = [(y, z) for i, (y, z) in enumerate(zip(y_booms[1:], z_booms[1:])) if z > 0]
    lower_booms = [(y, z) for i, (y, z) in enumerate(zip(y_booms[1:], z_booms[1:])) if z < 0]

    # Sort upper by increasing y (already is), lower by decreasing y to follow airfoil profile
    upper_booms.sort(key=lambda pt: pt[0])
    lower_booms.sort(key=lambda pt: -pt[0])

    # Full reordered list, starting at LE, around upper, then lower
    reordered_booms = [(y_booms[0], z_booms[0])] + upper_booms + lower_booms

    # Separate reordered y and z for distance calculation
    y_booms_ordered = [pt[0] for pt in reordered_booms]
    z_booms_ordered = [pt[1] for pt in reordered_booms]
    # print("y_booms_ordered :", y_booms_ordered)
    # print("z_booms_ordered :", z_booms_ordered)
    
    # Compute distances in correct order
    dist, max_dist, max_pair = dist_2_booms(y_booms_ordered, z_booms_ordered)
    
    # Reorder also the y_c_booms and z_c_booms in the same way
    # (normalised coordinates corresponding to reordered_booms)
    reordered_booms_c = [(y / chord_root, z / chord_root) for y, z in reordered_booms]

    y_c_booms_ordered = [pt[0] for pt in reordered_booms_c]
    z_c_booms_ordered = [pt[1] for pt in reordered_booms_c]
    #print('y_c_booms_ordered  =', y_c_booms_ordered)
    #print('z_c_booms_ordered  =', z_c_booms_ordered)
    
    
    # print('dist =', dist)
    # print("Max distance between successive booms: {:.3f} m".format(max_dist))

    # if not all(0.1 <= d <= 0.2 for d in dist):
    #     print("Warning: Some boom distances are outside the allowed range [0.1 m, 0.20 m], except if that is the vertical distance at the TE between cell 2 and 3.")
    #     print("Booms with max distance are at:")
    #     print(f"  Point a: y = {max_pair[0][0]:.3f}, z = {max_pair[0][1]:.3f}")
    #     print(f"  Point b: y = {max_pair[1][0]:.3f}, z = {max_pair[1][1]:.3f}")
    # else:
    #         print("All boom distances are within the required range.")
    
    
    # ---- Coordinate transformation to centroid-based system (for plotting) ----
    y_c_booms_centroid = [y - y_c_centroid for y in y_c_booms]
    z_c_booms_centroid = [z - z_c_centroid for z in z_c_booms]
    # Multiply by the chord to remove the normalization (centroid based frame of reference)
    y_booms_centroid =  np.array(y_c_booms_centroid) * chord_root
    z_booms_centroid =  np.array(z_c_booms_centroid) * chord_root
    # Multiply by the chord to remove the normalization (LE root based frame of reference)
    y_all = y_c_all*chord_root
    z_all = z_c_all*chord_root
    
    # ---- Coordinate transformation to centroid-based system (for the calculations) ----
    y_booms_ordered_centroid = [y - y_centroid for y in y_booms_ordered]
    z_booms_ordered_centroid = [z - z_centroid for z in z_booms_ordered]
    # print('')
    # print('y_booms_ordered_centroid =', y_booms_ordered_centroid)
    # print('')
    # print('z_booms_ordered_centroid =', z_booms_ordered_centroid)
    # print('')
    
    
    # ---- Delimitation of the cells ----
        
    # Define the y/c targets to delimit the different cells
    y_c_targets = [0.1, 0.3, 0.4, 0.6]
    
    # Convert booms to numpy arrays for easier processing
    y_c_booms_np = np.array(y_c_booms)
    z_c_booms_np = np.array(z_c_booms)
    
    # Store the real coordinates of the closest boom pair for each y/c target
    division_lines_yc = []
    division_lines_zc_upper = []
    division_lines_zc_lower = []
    division_lines_y = []
    division_lines_z_upper = []
    division_lines_z_lower = []
    
    # Find vertical lines based on closest booms
    for yc_target in y_c_targets:
        # Find all booms at that y/c (should be two: upper and lower)
        diffs = np.abs(y_c_booms_np - yc_target)
        idxs = np.where(diffs == np.min(diffs))[0]  # indices of closest boom(s)
    
        # Among the closest, separate upper and lower
        z_vals = z_c_booms_np[idxs]
        if len(z_vals) >= 2:
            z_upper = max(z_vals)
            z_lower = min(z_vals)
        else:
            # In case only one boom is found (rare), assume symmetry
            z_upper = z_vals[0]
            z_lower = -z_vals[0]
    
        y_closest = y_c_booms_np[idxs[0]]  # same for both booms
    
        # Append normalized values
        division_lines_yc.append(y_closest)
        division_lines_zc_upper.append(z_upper)
        division_lines_zc_lower.append(z_lower)
    
        # Convert to real units
        division_lines_y.append(y_closest * chord_root)
        division_lines_z_upper.append(z_upper * chord_root)
        division_lines_z_lower.append(z_lower * chord_root)

    # ---- Group booms by cell in centroid-based coordinate system ----
    
    y_booms_ordered = np.array(y_booms_ordered)
    z_booms_ordered = np.array(z_booms_ordered)
    
    yc_limits = [0.1, 0.3, 0.4, 0.6]
    # y_cutoffs = [0] + [chord_root * yc for yc in yc_limits] + [np.max(y_booms_ordered)]
    
    # Étape 1 : générer les coupures cibles
    raw_cutoffs = [0] + [chord_root * yc for yc in yc_limits] + [np.max(y_booms_ordered)]
    
    # Étape 2 : associer chaque cutoff à la valeur de boom la plus proche
    def get_closest_y(y_targets, y_available):
        return [y_available[np.argmin(np.abs(y_available - y))] for y in y_targets]
    
    y_cutoffs = get_closest_y(raw_cutoffs, y_booms_ordered)

    # Initialisation des sorties
    y_booms_cell_1 = []
    z_booms_cell_1 = []
    y_booms_cell_2 = []
    z_booms_cell_2 = []
    y_booms_cell_3 = []
    z_booms_cell_3 = []
    y_booms_cell_4 = []
    z_booms_cell_4 = []
    
    # Liste de référence pour traitement automatique
    y_outputs = [y_booms_cell_1, y_booms_cell_2, y_booms_cell_3, y_booms_cell_4]
    z_outputs = [z_booms_cell_1, z_booms_cell_2, z_booms_cell_3, z_booms_cell_4] 
    
    # Traitement par cellule
    for i in range(4):
        y_min, y_max = y_cutoffs[i], y_cutoffs[i + 1]
    
        indices = np.where((y_booms_ordered >= y_min) & (y_booms_ordered <= y_max))[0]
        
        y_c = y_booms_ordered[indices] - y_centroid
        z_c = z_booms_ordered[indices] - z_centroid
    
        y_outputs[i][:] = y_c
        z_outputs[i][:] = z_c
    
    # Conversion en np.array
    y_booms_cell_1 = np.array(y_booms_cell_1)
    z_booms_cell_1 = np.array(z_booms_cell_1)
    y_booms_cell_2 = np.array(y_booms_cell_2)
    z_booms_cell_2 = np.array(z_booms_cell_2)
    y_booms_cell_3 = np.array(y_booms_cell_3)
    z_booms_cell_3 = np.array(z_booms_cell_3)
    y_booms_cell_4 = np.array(y_booms_cell_4)
    z_booms_cell_4 = np.array(z_booms_cell_4)
    
    # print('y_booms_cell_1', y_booms_cell_1)
    # print('z_booms_cell_1', z_booms_cell_1)
    # print('y_booms_cell_2', y_booms_cell_2)
    # print('z_booms_cell_2', z_booms_cell_2)
    # print('y_booms_cell_3', y_booms_cell_3)
    # print('z_booms_cell_3', z_booms_cell_3)
    # print('y_booms_cell_4', y_booms_cell_4)
    # print('z_booms_cell_4', z_booms_cell_4)
 

    # ---- Boom area ---- (different pts of the envelope) 
    
    # Point A
    T_x_A = structural_loads_regular_wing(n[0], alpha[0], Lwt[0], Dwt[0], Wwt, Mw[0])[0]
    T_y_A = structural_loads_regular_wing(n[0], alpha[0], Lwt[0], Dwt[0], Wwt, Mw[0])[1]
    T_z_A = structural_loads_regular_wing(n[0], alpha[0], Lwt[0], Dwt[0], Wwt, Mw[0])[2]
    M_x_A = structural_loads_regular_wing(n[0], alpha[0], Lwt[0], Dwt[0], Wwt, Mw[0])[3]
    M_y_A = structural_loads_regular_wing(n[0], alpha[0], Lwt[0], Dwt[0], Wwt, Mw[0])[4]
    M_z_A = structural_loads_regular_wing(n[0], alpha[0], Lwt[0], Dwt[0], Wwt, Mw[0])[5]
    B_pt_A = boom_area(z_booms_ordered_centroid, y_booms_ordered_centroid, M_y_A, M_z_A,sigma_x_0, safety_factor)[0]
    sigma_xx_A = boom_area(z_booms_ordered_centroid, y_booms_ordered_centroid, M_y_A, M_z_A,sigma_x_0, safety_factor)[1]
    
    # Point B
    T_x_B = structural_loads_regular_wing(n[1], alpha[1], Lwt[1], Dwt[1], Wwt, Mw[1])[0]
    T_y_B = structural_loads_regular_wing(n[1], alpha[1], Lwt[1], Dwt[1], Wwt, Mw[1])[1]
    T_z_B = structural_loads_regular_wing(n[1], alpha[1], Lwt[1], Dwt[1], Wwt, Mw[1])[2]
    M_x_B = structural_loads_regular_wing(n[1], alpha[1], Lwt[1], Dwt[1], Wwt, Mw[1])[3]
    M_y_B = structural_loads_regular_wing(n[1], alpha[1], Lwt[1], Dwt[1], Wwt, Mw[1])[4]
    M_z_B = structural_loads_regular_wing(n[1], alpha[1], Lwt[1], Dwt[1], Wwt, Mw[1])[5]
    B_pt_B = boom_area(z_booms_ordered_centroid, y_booms_ordered_centroid, M_y_B, M_z_B,sigma_x_0, safety_factor)[0]
    sigma_xx_B = boom_area(z_booms_ordered_centroid, y_booms_ordered_centroid, M_y_B, M_z_B,sigma_x_0, safety_factor)[1]
    
    # Point C
    T_x_C = structural_loads_regular_wing(n[2], alpha[2], Lwt[2], Dwt[2], Wwt, Mw[2])[0]
    T_y_C = structural_loads_regular_wing(n[2], alpha[2], Lwt[2], Dwt[2], Wwt, Mw[2])[1]
    T_z_C = structural_loads_regular_wing(n[2], alpha[2], Lwt[2], Dwt[2], Wwt, Mw[2])[2]
    M_x_C = structural_loads_regular_wing(n[2], alpha[2], Lwt[2], Dwt[2], Wwt, Mw[2])[3]
    M_y_C = structural_loads_regular_wing(n[2], alpha[2], Lwt[2], Dwt[2], Wwt, Mw[2])[4]
    M_z_C = structural_loads_regular_wing(n[2], alpha[2], Lwt[2], Dwt[2], Wwt, Mw[2])[5]
    B_pt_C = boom_area(z_booms_ordered_centroid, y_booms_ordered_centroid, M_y_C, M_z_C,sigma_x_0, safety_factor)[0]
    sigma_xx_C = boom_area(z_booms_ordered_centroid, y_booms_ordered_centroid, M_y_C, M_z_C,sigma_x_0, safety_factor)[1]
    
    # Point D
    T_x_D = structural_loads_regular_wing(n[3], alpha[3], Lwt[3], Dwt[3], Wwt, Mw[3])[0]
    T_y_D = structural_loads_regular_wing(n[3], alpha[3], Lwt[3], Dwt[3], Wwt, Mw[3])[1]
    T_z_D = structural_loads_regular_wing(n[3], alpha[3], Lwt[3], Dwt[3], Wwt, Mw[3])[2]
    M_x_D = structural_loads_regular_wing(n[3], alpha[3], Lwt[3], Dwt[3], Wwt, Mw[3])[3]
    M_y_D = structural_loads_regular_wing(n[3], alpha[3], Lwt[3], Dwt[3], Wwt, Mw[3])[4]
    M_z_D = structural_loads_regular_wing(n[3], alpha[3], Lwt[3], Dwt[3], Wwt, Mw[3])[5]
    B_pt_D = boom_area(z_booms_ordered_centroid, y_booms_ordered_centroid, M_y_D, M_z_D,sigma_x_0, safety_factor)[0]
    sigma_xx_D = boom_area(z_booms_ordered_centroid, y_booms_ordered_centroid, M_y_D, M_z_D,sigma_x_0, safety_factor)[1]
    
    # Point E
    T_x_E = structural_loads_regular_wing(n[4], alpha[4], Lwt[4], Dwt[4], Wwt, Mw[4])[0]
    T_y_E = structural_loads_regular_wing(n[4], alpha[4], Lwt[4], Dwt[4], Wwt, Mw[4])[1]
    T_z_E = structural_loads_regular_wing(n[4], alpha[4], Lwt[4], Dwt[4], Wwt, Mw[4])[2]
    M_x_E = structural_loads_regular_wing(n[4], alpha[4], Lwt[4], Dwt[4], Wwt, Mw[4])[3]
    M_y_E = structural_loads_regular_wing(n[4], alpha[4], Lwt[4], Dwt[4], Wwt, Mw[4])[4]
    M_z_E = structural_loads_regular_wing(n[4], alpha[4], Lwt[4], Dwt[4], Wwt, Mw[4])[5]
    B_pt_E = boom_area(z_booms_ordered_centroid, y_booms_ordered_centroid, M_y_E, M_z_E,sigma_x_0, safety_factor)[0]
    sigma_xx_E = boom_area(z_booms_ordered_centroid, y_booms_ordered_centroid, M_y_E, M_z_E,sigma_x_0, safety_factor)[1]

    # Print the values : 
    N_to_lbf = 0.224809
    Nm_to_lbft = 0.737562

    points = ['A', 'B', 'C', 'D', 'E']

    for point in points:
        print(f"\n--- Point {point} ---")
        print(f"T_x_{point} = {eval(f'T_x_{point}'):.2f} N  = {eval(f'T_x_{point}') * N_to_lbf:,.2f} lbf")
        print(f"T_y_{point} = {eval(f'T_y_{point}'):.2f} N  = {eval(f'T_y_{point}') * N_to_lbf:,.2f} lbf")
        print(f"T_z_{point} = {eval(f'T_z_{point}'):.2f} N  = {eval(f'T_z_{point}') * N_to_lbf:,.2f} lbf")
        print(f"M_x_{point} = {eval(f'M_x_{point}'):.2f} N·m = {eval(f'M_x_{point}') * Nm_to_lbft:,.2f} lbf·ft")
        print(f"M_y_{point} = {eval(f'M_y_{point}'):.2f} N·m = {eval(f'M_y_{point}') * Nm_to_lbft:,.2f} lbf·ft")
        print(f"M_z_{point} = {eval(f'M_z_{point}'):.2f} N·m = {eval(f'M_z_{point}') * Nm_to_lbft:,.2f} lbf·ft")
        #print(f"B_pt_{point} = {eval(f'B_pt_{point}')}")


    # Max value (most critical)
    B_values = [B_pt_A, B_pt_B, B_pt_C, B_pt_D, B_pt_E]
    #print('')
    print('B_values =', B_values)
    sigma_xx_values = [sigma_xx_A, sigma_xx_B, sigma_xx_C, sigma_xx_D, sigma_xx_E]

    B = max(B_values)
    index_of_max = B_values.index(B)
    sigma_xx = sigma_xx_values[index_of_max]
    sigma_xx = np.array(sigma_xx)
      
    print(f"Index de la contrainte max : {index_of_max}")
    #print(f"sigma_xx associée (à B max) : {sigma_xx}")
    
    # Constants for conversion
    m2_to_mm2 = 1e6           # 1 m² = 1,000,000 mm²
    m2_to_in2 = 1550.0031     # 1 m² ≈ 1550.0031 in²
    print("")
    print(f"Boom area : {B:.9f} m² | {B * m2_to_mm2:.3f} mm² | {B * m2_to_in2:.3f} in²")
    print("")
    # print('The associated sigma_xx are :', sigma_xx) # should be an array for each boom
    
    
    # ---- Skin thickness ---- 
    
    # Taper effect along y (local taper of the transition wing, not the overall taper of the outerwing)
    y_pos_booms_local_tip = np.array(y_c_booms_ordered) * chord_tip + y_dist_root_tip 
    y_pos_booms_root = np.array(y_c_booms_ordered) * chord_root
    delta_y = y_pos_booms_local_tip - y_pos_booms_root
    delta_y = np.array(delta_y)
    # print('y_dist_root_tip  :', y_dist_root_tip )
    # print('')
    # print('delta_y :', delta_y) 
    # print('')
    
    # Taper effect along z
    z_pos_booms_local_tip = np.array(z_c_booms_ordered) * chord_tip 
    z_pos_booms_arrondi = np.array(z_c_booms_ordered) * chord_root # where chord_length_arrondi is the chord from the tip, just before the 'arrondi'
    delta_z = z_pos_booms_local_tip - z_pos_booms_arrondi 
    delta_z = np.array(delta_z)
    # print('delta_z :', delta_z)
    # print('')
    
    # # Point A
    thickness_pt_A = skin_thickness(B, sigma_xx, delta_x, delta_y, delta_z, T_y_A, T_z_A, M_x_A, y_booms_ordered_centroid, z_booms_ordered_centroid, y_booms_cell_1, z_booms_cell_1, tau_max, y_booms_cell_2, z_booms_cell_2, y_centroid, z_centroid, y_booms_cell_3, z_booms_cell_3, y_booms_cell_4, z_booms_cell_4)
    
    # # Point B
    thickness_pt_B = skin_thickness(B, sigma_xx, delta_x, delta_y, delta_z, T_y_B, T_z_B, M_x_B, y_booms_ordered_centroid, z_booms_ordered_centroid, y_booms_cell_1, z_booms_cell_1, tau_max, y_booms_cell_2, z_booms_cell_2, y_centroid, z_centroid, y_booms_cell_3, z_booms_cell_3, y_booms_cell_4, z_booms_cell_4)
    
    # # Point C
    thickness_pt_C = skin_thickness(B, sigma_xx, delta_x, delta_y, delta_z, T_y_C, T_z_C, M_x_C, y_booms_ordered_centroid, z_booms_ordered_centroid, y_booms_cell_1, z_booms_cell_1, tau_max, y_booms_cell_2, z_booms_cell_2, y_centroid, z_centroid, y_booms_cell_3, z_booms_cell_3, y_booms_cell_4, z_booms_cell_4)
    
    # # Point D
    thickness_pt_D = skin_thickness(B, sigma_xx, delta_x, delta_y, delta_z, T_y_D, T_z_D, M_x_D, y_booms_ordered_centroid, z_booms_ordered_centroid, y_booms_cell_1, z_booms_cell_1, tau_max, y_booms_cell_2, z_booms_cell_2, y_centroid, z_centroid, y_booms_cell_3, z_booms_cell_3, y_booms_cell_4, z_booms_cell_4)
    
    # # Point E
    thickness_pt_E = skin_thickness(B, sigma_xx, delta_x, delta_y, delta_z, T_y_E, T_z_E, M_x_E, y_booms_ordered_centroid, z_booms_ordered_centroid, y_booms_cell_1, z_booms_cell_1, tau_max, y_booms_cell_2, z_booms_cell_2, y_centroid, z_centroid, y_booms_cell_3, z_booms_cell_3, y_booms_cell_4, z_booms_cell_4)
    
    # # Take the maximum value
    thickness_values = [thickness_pt_A, thickness_pt_B, thickness_pt_C, thickness_pt_D, thickness_pt_E]
    print('Thickness_values =', thickness_values)
    thickness = max(thickness_values)

    # # Conversion constants
    m_to_mm = 1000           # 1 m = 1000 mm
    m_to_in = 39.3701        # 1 m ≈ 39.3701 in
    print('')
    print(f"Skin thickness : {thickness:.9f} m | {thickness * m_to_mm:.6f} mm | {thickness * m_to_in:.6f} in")
    print('')
    
    
    # ---- Plotting ----
    if not plot_airfoil:
        return
    
    fig, ax = plt.subplots(figsize=(11, 5),dpi=300)
    airfoil_line, = ax.plot(y_c_all, z_c_all)  # Airfoil plot
    ax.scatter(y_c_booms, z_c_booms, color='red', zorder=3)

    # Get the color of the airfoil plot to use for the dotted line
    airfoil_color = airfoil_line.get_color()    
    
    # Plot all vertical cell delimiters
    for yc, z_top, z_bot in zip(division_lines_yc, division_lines_zc_upper, division_lines_zc_lower):
        ax.plot([yc, yc], [z_top, z_bot], '--', color=airfoil_color)

    
    # Plot the centroid location as a blue point
    ax.scatter(y_c_centroid, z_c_centroid, color='blue', zorder=4)
    
    # Add "C" label next to the centroid
    ax.text(y_c_centroid - 0.05, z_c_centroid - 0.025, 'C', color='blue', fontsize=12)
    
    # Plot coordinate axes at the centroid with arrows
    ax.arrow(y_c_centroid, z_c_centroid, 0.05, 0, head_width=0.02, head_length=0.02, fc='blue', ec='blue')  # X-axis arrow
    ax.text(y_c_centroid + 0.08, z_c_centroid - 0.01, 'y', color='blue', fontsize=12) # Y-axis name
    ax.arrow(y_c_centroid, z_c_centroid, 0, 0.05, head_width=0.02, head_length=0.02, fc='blue', ec='blue')  # Z-axis arrow
    ax.text(y_c_centroid - 0.01, z_c_centroid + 0.08, 'z', color='blue', fontsize=12) # Z-axis name
    
    # Plot max distance segment in orange 
    ax.plot([max_pair[0][0]/chord_root, max_pair[1][0]/chord_root],
            [max_pair[0][1]/chord_root, max_pair[1][1]/chord_root],
            linestyle='--', color='orange', linewidth=2)
        
    # Set equal scaling for both axes (same units for y and z axes)
    ax.set_aspect('equal')
    
    # Labels and title
    ax.set_xlabel('$y/c$')
    ax.set_ylabel('$z/c$')
    ax.set_title(f'SC(2)-0710 Airfoil with {n_booms} booms')
    ax.grid()
    #fig.savefig(fname='airfoil_transition_y_c_centroid_location.pdf')
    plt.show()

    
    # ---- Plotting (Real Size) ----
    fig2, ax2 = plt.subplots(figsize=(11, 5), dpi=300)
    ax2.plot(y_all, z_all, label='Airfoil')  # Airfoil in meters
    ax2.scatter(y_booms, z_booms, color='red', zorder=3)

    # Plot all vertical delimiters in real coordinates
    for y, z_top, z_bot in zip(division_lines_y, division_lines_z_upper, division_lines_z_lower):
        ax2.plot([y, y], [z_top, z_bot], '--', color=airfoil_color)
    

    # Plot centroid in real coordinates
    ax2.scatter(y_centroid, z_centroid, color='blue', zorder=4)
    ax2.text(y_centroid - 0.25, z_centroid - 0.125, 'C', color='blue', fontsize=12)

    # Plot coordinate axes at the centroid
    ax2.arrow(y_centroid, z_centroid, 0.05 * chord_root, 0, head_width=0.02 *chord_root, head_length=0.02 * chord_root, fc='blue', ec='blue')
    ax2.text(y_centroid + 0.4, z_centroid - 0.05, 'y', color='blue', fontsize=12)
    ax2.arrow(y_centroid, z_centroid, 0, 0.05 * chord_root, head_width=0.02 * chord_root, head_length=0.02 * chord_root, fc='blue', ec='blue')
    ax2.text(y_centroid - 0.05 , z_centroid + 0.4, 'z', color='blue', fontsize=12)

    ax2.set_aspect('equal')
    ax2.set_xlabel('y [m]')
    ax2.set_ylabel('z [m]')
    ax2.set_title(f'SC(2)-0710 Airfoil in Real Scale ({n_booms} booms)')
    ax2.grid()
    #ax2.legend()
    fig2.savefig(fname='airfoil_transition_y_centroid_location.pdf')
    plt.show()


    # ---- Third Plot: Minimal Clean Version ----
    fig3, ax3 = plt.subplots(figsize=(11, 5), dpi=300)
    ax3.plot(y_all, z_all, color='#ff914d')  # Airfoil outline in black
    ax3.scatter(y_booms, z_booms, color='black', s=11, zorder=3)  # Booms as small black dots

    # Plot all division lines in minimal plot
    for y, z_top, z_bot in zip(division_lines_y, division_lines_z_upper, division_lines_z_lower):
        ax3.plot([y, y], [z_top, z_bot], '--', color='orange', linewidth=1)

    # Plot the centroid as a small point only
    ax3.scatter(y_centroid, z_centroid, color='black', s=8, zorder=2)

    # Remove extras
    ax3.set_aspect('equal')
    ax3.axis('off')  # No axis, ticks, grid, labels, or title

    # Save the clean plot
    fig3.savefig('airfoil_transition_minimal_clean.pdf', bbox_inches='tight', pad_inches=0)
    plt.show()

# -------------------------------------------------------------------

plotAirfoil(True, n_booms=72) 


