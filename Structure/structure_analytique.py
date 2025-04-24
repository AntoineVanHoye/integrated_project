

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import math


#geometry wings
beta = np.radians(4.3) # Setting angle [rad] 
(dx_ac,dy_ac,dz_ac) = (13.193-11.27477,9.871-(4.5+3.109),0) # Coord AC of the wing ()
(dx_w_cg,dy_w_cg,dz_w_cg) = (3.2418235, 2.6913599, 0.1539403) #y placement of the center of gravity of the wing to the root of the regular wing
#(dxa_emp,dya_emp,dza_emp) = (0,0,0) #distance of the aerodynamic center to the root wing of the empenage


# ---- Structural loads ---- (will be imported from another fct, here just expl values)

# Here, import all the n and all the alpha from the envelope : 
n =  [2.5, 2.5, 0, -1, -1] # Load factor
alpha =  [0.15358704, -0.01836998, -0.09427013, -0.13817235, -0.29559977]   # [rad]
L_overall = [1080014.79567327, 1159506.87070659, 137839.72949363, -289992.70212547, -319343.31241878] # [N]
D_wing =  [4827.34564335, 15528.48125026, 15528.48125026, 9938.22800017, 1927.89480254] # [N]

# Trouver un moyen de calculer L_wing via la formule analytique avec C_L = 0.433 via pondération des surfaces 
# Drag : pondérer par les surfaces aussi 
Lwt = np.array(L_overall) # Extract the correct lift
Dwt = np.array(D_wing)
Wwt = (2858.58)*9.81 # [kg] mass -> to force [N] 
Mw = [-248806.10076367, -800353.06276004, -800353.06276004, -512225.96016643, -99365.577678] # [N.m]


# ---- Material ----

sigma_y_0 = 1500*10**6 # CHOOSE THE MATERIAL
tau_max = 94*10**6 # maximum shear stress
safety_factor = 1.5


# -------------------------------------------------------------------

def structural_loads_regular_wing (n, alpha, Lwt, Dwt, Wwt, Mw): # There will be more parameters as the lift also varies, etc
    
    T_x = (n*Wwt/2-Lwt/2)*np.sin(alpha+beta) + Dwt/2*np.cos(alpha+beta) #alpha is the angle of attack, beta is the setting angle
    T_y = 0
    T_z = (-n*Wwt/2+Lwt/2)*np.cos(alpha+beta) + Dwt/2*np.sin(alpha+beta)
    M_x = 1/2*(-n*(Wwt*dy_w_cg)+Lwt*dy_ac)*np.cos(alpha+beta) - Dwt/2*dy_ac*np.sin(alpha+beta)#dyw is the y placement of the center of gravity to the root of the regular wing , dya is the distance of the aerodynamic center to the root wing
    M_y = 1/2*(-n*Wwt*dx_w_cg + Lwt*dx_ac -Dwt*dz_ac)*np.cos(alpha+beta) + 1/2*(-n*Wwt*dz_w_cg +Lwt*dz_ac + Dwt*dx_ac)*np.sin(alpha+beta) + Mw/2 
    M_z = 1/2*(+n*(Wwt*dy_w_cg)-Lwt*dy_ac)*np.sin(alpha+beta) + Dwt/2*dy_ac*np.cos(alpha+beta)
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

def polygon_area(x_coords, y_coords): # Fct to compute the area of the cells
    n = len(x_coords)
    area = 0.0

    for i in range(n):
        j = (i + 1) % n  # Wraps around to connect last to first
        area += x_coords[i] * y_coords[j]
        area -= x_coords[j] * y_coords[i]

    return abs(area) / 2.0

# -------------------------------------------------------------------

def swept_area_from_center(x_points, z_points, x_centroid, z_centroid):
    n = len(x_points)
    total_area = 0.0
    individual_areas = []

    for i in range(n):
        x1, z1 = x_points[i], z_points[i]
        x2, z2 = x_points[(i + 1) % n], z_points[(i + 1) % n]  # wrap around

        # Area computation
        area = 0.5 * abs(x1 * z2 + x2 * z_centroid + x_centroid * z1 - z1 * x2 - z2 * x_centroid - z_centroid * x1)
        individual_areas.append(area)
        total_area += area

    return total_area, individual_areas

# -------------------------------------------------------------------

def boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x, M_z, sigma_y_0, safety_factor) : 
    
    # Inertia per unit area (the boom area 'B' is considered to be the same everywhere)
    I_xx_over_B = np.sum(np.array(z_booms_ordered_centroid)**2) 
    I_zz_over_B = np.sum(np.array(x_booms_ordered_centroid)**2) 
    I_xz_over_B = np.sum(np.array(x_booms_ordered_centroid) * np.array(z_booms_ordered_centroid)) 
    
    # Direct stress in each boom (check formula)
    denom = (I_xx_over_B * I_zz_over_B - I_xz_over_B**2)  # Denominator of the stress equation
    # Compute stress for each boom
    B_sigma_yy = []
    for x, z in zip(x_booms_ordered_centroid, z_booms_ordered_centroid):
        stress = ((I_zz_over_B * M_x + I_xz_over_B * M_z) * z - (I_xz_over_B * M_x + I_xx_over_B * M_z) * x) / denom
        B_sigma_yy.append(stress)
    
    # Find maximum stress
    max_stress = max(B_sigma_yy)  # Find maximum stress value
    max_stress_index = B_sigma_yy.index(max_stress)  # Find the index of the maximum stress
    
    
    # ---- Minimum Area -----
    sigma_yy_max_material = sigma_y_0 / safety_factor
    B_min = max_stress/sigma_yy_max_material
    
    # Real stress in the boom with the new value of B 
    sigma_yy = np.array(B_sigma_yy) / B_min
    
    return B_min, sigma_yy
    
# -------------------------------------------------------------------

def skin_thickness(B, sigma_yy, delta_y, delta_x, delta_z, T_x, T_z, M_y, x_booms_ordered_centroid, z_booms_ordered_centroid, x_booms_cell_1, z_booms_cell_1, tau_max, x_booms_cell_2, z_booms_cell_2, x_centroid, z_centroid) :
    
    # Taper effect (suite) :
    T_x_web = T_x - B* np.sum(sigma_yy * (delta_x/delta_y))
    T_z_web = T_z - B* np.sum(sigma_yy * (delta_z/delta_y))
    
    # Inertia per unit area (the boom area 'B' is considered to be the same everywhere)
    I_xx_over_B = np.sum(np.array(z_booms_ordered_centroid)**2) 
    I_zz_over_B = np.sum(np.array(x_booms_ordered_centroid)**2) 
    I_xz_over_B = np.sum(np.array(x_booms_ordered_centroid) * np.array(z_booms_ordered_centroid)) 
    
    # ---- Open shear flow ----
    denom = (I_xx_over_B*I_zz_over_B - I_xz_over_B**2)
    factor_1 = (I_zz_over_B*T_z_web - I_xz_over_B*T_x_web)/denom
    factor_2 = (I_xx_over_B*T_x_web - I_xz_over_B*T_z_web)/denom
    
    
    # -- Cell 1 --
    q_0_cell_1 = [0] # Cut in cell 1, the open shear flow between booms 7 and 6 is zero                                          (q_0_7,6)
    q_0_cell_1.append(q_0_cell_1[0] - factor_1 * z_booms_ordered_centroid[6-1] - factor_2 * x_booms_ordered_centroid[6-1])     # (q_0_6,5)
    q_0_cell_1.append(q_0_cell_1[1] - factor_1 * z_booms_ordered_centroid[5-1] - factor_2 * x_booms_ordered_centroid[5-1])     # (q_0_5,4)
    q_0_cell_1.append(q_0_cell_1[2] - factor_1 * z_booms_ordered_centroid[4-1] - factor_2 * x_booms_ordered_centroid[4-1])     # (q_0_4,3)
    q_0_cell_1.append(q_0_cell_1[3] - factor_1 * z_booms_ordered_centroid[3-1] - factor_2 * x_booms_ordered_centroid[3-1])     # (q_0_3,2)
    q_0_cell_1.append(q_0_cell_1[4] - factor_1 * z_booms_ordered_centroid[2-1] - factor_2 * x_booms_ordered_centroid[2-1])     # (q_0_2,1)
    q_0_cell_1.append(q_0_cell_1[5] - factor_1 * z_booms_ordered_centroid[1-1] - factor_2 * x_booms_ordered_centroid[1-1])     # (q_0_1,37)
    q_0_cell_1.append(q_0_cell_1[6] - factor_1 * z_booms_ordered_centroid[37-1] - factor_2 * x_booms_ordered_centroid[37-1])   # (q_0_37,36)
    q_0_cell_1.append(q_0_cell_1[7] - factor_1 * z_booms_ordered_centroid[36-1] - factor_2 * x_booms_ordered_centroid[36-1])   # (q_0_36,35)
    q_0_cell_1.append(q_0_cell_1[8] - factor_1 * z_booms_ordered_centroid[35-1] - factor_2 * x_booms_ordered_centroid[35-1])   # (q_0_35,34)
    q_0_cell_1.append(q_0_cell_1[9] - factor_1 * z_booms_ordered_centroid[34-1] - factor_2 * x_booms_ordered_centroid[34-1])   # (q_0_34,33)
    q_0_cell_1.append(q_0_cell_1[10] - factor_1 * z_booms_ordered_centroid[33-1] - factor_2 * x_booms_ordered_centroid[33-1])  # (q_0_33,32) 
    q_spar_1_2 = - factor_1 * z_booms_ordered_centroid[7-1] - factor_2 * x_booms_ordered_centroid[7-1]  # (q_0_7,32) et pas (q_0_32,7) !
    q_0_cell_1.append(-q_spar_1_2)  # (q_0_32,7) mettre dans le meme sens que les autres
    
    # Mettre le tableau dans le meme ordre que x_booms_cell_1 et z_booms_cell_1
    n = 6  # nombre d'éléments à déplacer à la fin 
    q_0_cell_1 = q_0_cell_1[n:] + q_0_cell_1[:n] # Cut the array between case n=6 and 7 and invert the 2 part
    # Invert the order of all elements
    q_0_cell_1 = q_0_cell_1[::-1] 
    # The array "q_0_cell_1_t" is now in the same order than x_booms_cell_1 and z_booms_cell_1
    
    
    # -- Cell 2 (t) -- 
    q_0_cell_2 = [0]                                                                                                            # (q_0_7,8)
    q_0_cell_2.append(q_0_cell_2[0] - factor_1 * z_booms_ordered_centroid[8-1] - factor_2 * x_booms_ordered_centroid[8-1])    # (q_0_8,9)
    q_0_cell_2.append(q_0_cell_2[1] - factor_1 * z_booms_ordered_centroid[9-1] - factor_2 * x_booms_ordered_centroid[9-1])    # (q_0_9,10)
    q_0_cell_2.append(q_0_cell_2[2] - factor_1 * z_booms_ordered_centroid[10-1] - factor_2 * x_booms_ordered_centroid[10-1])  # (q_0_10,11)
    q_0_cell_2.append(q_0_cell_2[3] - factor_1 * z_booms_ordered_centroid[11-1] - factor_2 * x_booms_ordered_centroid[11-1])  # (q_0_11,12)
    q_0_cell_2.append(q_0_cell_2[4] - factor_1 * z_booms_ordered_centroid[12-1] - factor_2 * x_booms_ordered_centroid[12-1])  # (q_0_12,13)
    q_0_cell_2.append(q_0_cell_2[5] - factor_1 * z_booms_ordered_centroid[13-1] - factor_2 * x_booms_ordered_centroid[13-1])  # (q_0_13,14)
    q_0_cell_2.append(q_0_cell_2[6] - factor_1 * z_booms_ordered_centroid[14-1] - factor_2 * x_booms_ordered_centroid[14-1])  # (q_0_14,15)
    q_0_cell_2.append(q_0_cell_2[7] - factor_1 * z_booms_ordered_centroid[15-1] - factor_2 * x_booms_ordered_centroid[15-1])  # (q_0_15,16)
    q_0_cell_2.append(q_0_cell_2[8] - factor_1 * z_booms_ordered_centroid[16-1] - factor_2 * x_booms_ordered_centroid[16-1])  # (q_0_16,17)
    q_0_cell_2.append(q_0_cell_2[9] - factor_1 * z_booms_ordered_centroid[17-1] - factor_2 * x_booms_ordered_centroid[17-1])  # (q_0_17,18)
    q_0_cell_2.append(q_0_cell_2[10] - factor_1 * z_booms_ordered_centroid[18-1] - factor_2 * x_booms_ordered_centroid[18-1]) # (q_0_18,19)
    q_0_cell_2.append(q_0_cell_2[11] - factor_1 * z_booms_ordered_centroid[19-1] - factor_2 * x_booms_ordered_centroid[19-1]) # (q_0_19,20)
    q_0_cell_2.append(q_0_cell_2[12] - factor_1 * z_booms_ordered_centroid[20-1] - factor_2 * x_booms_ordered_centroid[20-1]) # (q_0_20,21)
    q_0_cell_2.append(q_0_cell_2[13] - factor_1 * z_booms_ordered_centroid[21-1] - factor_2 * x_booms_ordered_centroid[21-1]) # (q_0_21,22)
    q_0_cell_2.append(q_0_cell_2[14] - factor_1 * z_booms_ordered_centroid[22-1] - factor_2 * x_booms_ordered_centroid[22-1]) # (q_0_22,23)
    q_0_cell_2.append(q_0_cell_2[15] - factor_1 * z_booms_ordered_centroid[23-1] - factor_2 * x_booms_ordered_centroid[23-1]) # (q_0_23,24)
    q_0_cell_2.append(q_0_cell_2[16] - factor_1 * z_booms_ordered_centroid[24-1] - factor_2 * x_booms_ordered_centroid[24-1]) # (q_0_24,25)
    q_0_cell_2.append(q_0_cell_2[17] - factor_1 * z_booms_ordered_centroid[25-1] - factor_2 * x_booms_ordered_centroid[25-1]) # (q_0_25,26)
    q_0_cell_2.append(q_0_cell_2[18] - factor_1 * z_booms_ordered_centroid[26-1] - factor_2 * x_booms_ordered_centroid[26-1]) # (q_0_26,27)
    q_0_cell_2.append(q_0_cell_2[19] - factor_1 * z_booms_ordered_centroid[27-1] - factor_2 * x_booms_ordered_centroid[27-1]) # (q_0_27,28)
    q_0_cell_2.append(q_0_cell_2[20] - factor_1 * z_booms_ordered_centroid[28-1] - factor_2 * x_booms_ordered_centroid[28-1]) # (q_0_28,29)
    q_0_cell_2.append(q_0_cell_2[21] - factor_1 * z_booms_ordered_centroid[29-1] - factor_2 * x_booms_ordered_centroid[29-1]) # (q_0_29,30)
    q_0_cell_2.append(q_0_cell_2[22] - factor_1 * z_booms_ordered_centroid[30-1] - factor_2 * x_booms_ordered_centroid[30-1]) # (q_0_30,31)
    q_0_cell_2.append(q_0_cell_2[23] - factor_1 * z_booms_ordered_centroid[31-1] - factor_2 * x_booms_ordered_centroid[31-1]) # (q_0_31,32)   
    q_0_cell_2.append(-q_spar_1_2) # (q_0_32,7)
    
    # Invert the sign of all elements as I went clockwise for cell 2 
    q_0_cell_2 = [-x for x in q_0_cell_2]
    
    
    # ---- Lengths of the segments and the cells ----
    lengths = dist_2_booms(x_booms_ordered_centroid, z_booms_ordered_centroid)[0]
    intersecting_length = dist_2_pts(x_booms_ordered_centroid[7-1], z_booms_ordered_centroid[7-1], x_booms_ordered_centroid[32-1], z_booms_ordered_centroid[32-1])

    # Cell 1
    lengths_c_1 = dist_2_booms(x_booms_cell_1, z_booms_cell_1)[0]
    l_cell_1 = np.sum(lengths_c_1)
    # print('lengths_c_1 = ', lengths_c_1)
    # print(l_cell_1)
    
    # Cell 2
    lengths_c_2 = dist_2_booms(x_booms_cell_2, z_booms_cell_2)[0]
    l_cell_2 = np.sum(lengths_c_2)
    # print('lengths_c_2 = ', lengths_c_2)
    # print(l_cell_2)
    
    
    # ---- Cell Area ----
    A_h_c_1 = polygon_area(x_booms_cell_1, z_booms_cell_1)
    A_h_c_2 = polygon_area(x_booms_cell_2, z_booms_cell_2)
    #print(f"Cell Areas: A_h_c_1 = {A_h_c_1:.4f} m², A_h_c_2 = {A_h_c_2:.4f} m²")
    
    
    # ---- Open shear flow integration ----
    # Cell 1
    q_0_cell_1 = np.array(q_0_cell_1) # Convert the lists to NumPy arrays
    lengths_c_1 = np.array(lengths_c_1)
    # print('q_0_cell_1 =', q_0_cell_1)
    # print('lengths_c_1 =', lengths_c_1)
    open_shear_flux_1 = np.sum(q_0_cell_1 * lengths_c_1)
    # print('open_shear_flux_1 =', open_shear_flux_1)
    
    # Cell 2
    q_0_cell_2 = np.array(q_0_cell_2) # Convert the lists to NumPy arrays
    lengths_c_2 = np.array(lengths_c_2)
    # print('q_0_cell_2 =', q_0_cell_2)
    # print('lengths_c_2 =', lengths_c_2)
    open_shear_flux_2 = np.sum(q_0_cell_2 * lengths_c_2)
    # print('open_shear_flux_2 =', open_shear_flux_2)
    
    
    # ---- Swept Areas ----
    # Cell 1
    swept_area_cell_1 = swept_area_from_center(x_booms_cell_1, z_booms_cell_1, x_centroid = 0, z_centroid = 0)[1]
    # Cell 2
    swept_area_cell_2 = swept_area_from_center(x_booms_cell_2, z_booms_cell_2, x_centroid = 0, z_centroid = 0)[1]
    # Term in the momentum equilibrium : 
    term_swept_area = 2 * (np.sum(q_0_cell_1 * swept_area_cell_1) + np.sum(q_0_cell_2 * swept_area_cell_2))
    #print('term_swept_area =', term_swept_area)
    
    
    # ---- Correction term computation : solve a 2 eqns syst (because 2 cells) ----
    # Syst parameters : x => q_1_corr and y => q_2_corr
    
    # Twist rate comaptibility equation
    coeff_x_eq1 = l_cell_1 * A_h_c_2 + intersecting_length * A_h_c_1
    coeff_y_eq1 = (-1)* intersecting_length * A_h_c_2 - l_cell_2 * A_h_c_2
    ind_term_eq1 = open_shear_flux_2 * A_h_c_1 - open_shear_flux_1 * A_h_c_2
    
    # Momentum equilibrium
    coeff_x_eq2 = 2 * A_h_c_1
    coeff_y_eq2 = 2 * A_h_c_2
    
    Moment_shear = x_centroid * T_z - z_centroid * T_x 
    term_P_z = B * np.sum(x_booms_ordered_centroid * sigma_yy * delta_z/delta_y)
    term_P_x = B * np.sum(x_booms_ordered_centroid * sigma_yy * delta_x/delta_y)
    ind_term_eq2 = Moment_shear - term_swept_area - term_P_z + term_P_x 
    
    # Solve syst : 
    A = np.array([[coeff_x_eq1, coeff_y_eq1],[coeff_x_eq2, coeff_y_eq2]]) # Coeff matrix
    B = np.array([ind_term_eq1, ind_term_eq2]) # Independent term vector 
    sol = np.linalg.solve(A, B)
    q_1_corr, q_2_corr = sol
    #print(f"q_1_corr = {q_1_corr:.4f}, q_2_corr = {q_2_corr:.4f}")
    
    # ---- Shear flow due to torsion ----
    # Coeff of the eq of the syst :
    a_1 = 2 * A_h_c_1
    b_1 = 2 * A_h_c_2
    c_1 = -M_y
    
    a_2 = (l_cell_1/A_h_c_1) + (intersecting_length/A_h_c_2)
    b_2 = -((l_cell_2/A_h_c_2) + (intersecting_length/A_h_c_1))
    c_2 = 0
    
    # Solve syst : 
    A_torsion = np.array([[a_1, b_1],[a_2, b_2]]) # Coeff matrix
    B_torsion = np.array([c_1, c_2]) # Independent term vector 
    sol_torsion = np.linalg.solve(A_torsion, B_torsion)
    q_1_torsion, q_2_torsion = sol_torsion
    #print(f"q_1_torsion = {q_1_torsion:.4f}, q_2_torsion = {q_2_torsion:.4f}")
    
    # ---- Shear flow (closed) ----
    q_closed_cell_1 = np.array(q_0_cell_1) + q_1_corr + q_1_torsion # Apply the correction and the shear flow due to torsion
    q_closed_cell_2 = np.array(q_0_cell_2) + q_2_corr + q_2_torsion 
    q_closed = np.concatenate([q_closed_cell_1, q_closed_cell_2])
    
    # ---- Thickness computation ----   
    thickness = max(q_closed, key=abs)/tau_max 
    
    return thickness

# -------------------------------------------------------------------

def plotAirfoil(plot_airfoil, n_booms, x_c_max, x_c_cell_1):
    
    # ---- Geometry ----
    wing_semi_span = 10 #[m]
    y_tip = wing_semi_span
    y_root = 0
    y_arrondi = 3.10921985816
    chord_root = 5.875
    chord_tip = 2.35
    chord_length_arrondi = chord_tip + (chord_root - chord_tip)/(y_root - y_tip) * (y_arrondi - y_tip)
    #print('chord_length :',chord_length_arrondi)
    sweep_angle = np.radians(31.599) #[°]
    delta_y = (wing_semi_span - y_arrondi)
    x_dist_root_tip = delta_y *np.tan(sweep_angle)
    #print('delta_y : ', delta_y)
    #print('x_dist_root_tip : ', x_dist_root_tip)
    #print('')
    
    # Load airfoil data
    data_wing = pd.read_csv("sc20710_XYZ.csv")
    x_c_all = data_wing.iloc[:, 0].values
    z_c_all = data_wing.iloc[:, 1].values  # Change y to z (x, z coordinates)

    # Separate and sort upper and lower surfaces
    upper_surface = data_wing[data_wing.iloc[:, 1] >= 0].copy().sort_values(by=data_wing.columns[0])  # z >= 0
    lower_surface = data_wing[data_wing.iloc[:, 1] < 0].copy().sort_values(by=data_wing.columns[0])  # z < 0

    # Add the leading edge boom at (0, 0)
    x_c_booms = [0.0]
    z_c_booms = [0.0]

    # Ensure odd number of booms for symmetry around leading edge
    if (n_booms - 1) % 2 != 0:
        raise ValueError("Number of booms must be odd to have (0,0) + symmetric pairs.")

    n_pairs = (n_booms - 1) // 2


    # Generate x positions from 0 to x_c_max, skipping the leading edge (already added)
    x_c_positions = np.linspace(0, x_c_max, n_pairs + 1)[1:]

    # Interpolate z from upper and lower surfaces
    z_c_upper = np.interp(x_c_positions, upper_surface.iloc[:, 0], upper_surface.iloc[:, 1])
    z_c_lower = np.interp(x_c_positions, lower_surface.iloc[:, 0], lower_surface.iloc[:, 1])

    # Add symmetric booms
    for x, z_u, z_l in zip(x_c_positions, z_c_upper, z_c_lower):
        x_c_booms.extend([x, x])
        z_c_booms.extend([z_u, z_l])


    # ---- Find the points closest to x/c = 0.25 for the dotted line ----
    x_c_dotted_target_1 = x_c_cell_1
    closest_idx_upper_1 = np.argmin(np.abs(x_c_positions - x_c_dotted_target_1))
    closest_x_c_upper_1 = x_c_positions[closest_idx_upper_1]
    z_c_upper_dotted_1 = np.interp(closest_x_c_upper_1, upper_surface.iloc[:, 0], upper_surface.iloc[:, 1])
    z_c_lower_dotted_1 = np.interp(closest_x_c_upper_1, lower_surface.iloc[:, 0], lower_surface.iloc[:, 1])


    # ---- Find the points closest to x/c = x_c_max for the dotted line ----
    x_c_dotted_target_2 = x_c_max
    closest_idx_upper_2 = np.argmin(np.abs(x_c_positions - x_c_dotted_target_2))
    closest_x_c_upper_2 = x_c_positions[closest_idx_upper_2]
    z_c_upper_dotted_2 = np.interp(closest_x_c_upper_2, upper_surface.iloc[:, 0], upper_surface.iloc[:, 1])
    z_c_lower_dotted_2 = np.interp(closest_x_c_upper_2, lower_surface.iloc[:, 0], lower_surface.iloc[:, 1])


    # ---- Centroid location ----
    x_c_centroid = np.sum(x_c_booms) / len(x_c_booms)
    z_c_centroid = np.sum(z_c_booms) / len(z_c_booms)
    # Multiply by the chord
    x_booms = np.array(x_c_booms) * chord_length_arrondi
    z_booms = np.array(z_c_booms) * chord_length_arrondi
    x_centroid = x_c_centroid * chord_length_arrondi
    z_centroid = z_c_centroid * chord_length_arrondi
    print('')
    #print('Centroid location:')
    #print(f'x = {x_centroid} m')
    #print(f'z = {z_centroid} m')
    # Conversion to feet
    print(f'x_centroid = {x_centroid*3.28084} ft')
    print(f'z_centroid = {z_centroid*3.28084} ft')
    
    
    # ---- Verify that the distance between 2 succesive booms is between 0.1 and 0.2 m (Limits suggested by M.Noels) ----
    # Separate booms on upper and lower surfaces (skipping the first LE point)
    upper_booms = [(x, z) for i, (x, z) in enumerate(zip(x_booms[1:], z_booms[1:])) if z > 0]
    lower_booms = [(x, z) for i, (x, z) in enumerate(zip(x_booms[1:], z_booms[1:])) if z < 0]

    # Sort upper by increasing x (already is), lower by decreasing x to follow airfoil profile
    upper_booms.sort(key=lambda pt: pt[0])
    lower_booms.sort(key=lambda pt: -pt[0])

    # Full reordered list, starting at LE, around upper, then lower
    reordered_booms = [(x_booms[0], z_booms[0])] + upper_booms + lower_booms

    # Separate reordered x and z for distance calculation
    x_booms_ordered = [pt[0] for pt in reordered_booms]
    z_booms_ordered = [pt[1] for pt in reordered_booms]
    
    # Compute distances in correct order
    dist, max_dist, max_pair = dist_2_booms(x_booms_ordered, z_booms_ordered)
    
    # Reorder also the x_c_booms and z_c_booms in the same way
    # (normalised coordinates corresponding to reordered_booms)
    reordered_booms_c = [(x / chord_length_arrondi, z / chord_length_arrondi) for x, z in reordered_booms]

    x_c_booms_ordered = [pt[0] for pt in reordered_booms_c]
    z_c_booms_ordered = [pt[1] for pt in reordered_booms_c]
    #print('x_c_booms_ordered  =', x_c_booms_ordered)
    #print('z_c_booms_ordered  =', z_c_booms_ordered)
    
    '''
    print('dist =', dist)
    print("Max distance between successive booms: {:.3f} m".format(max_dist))

    if not all(0.1 <= d <= 0.2 for d in dist):
        print("⚠️ Warning: Some boom distances are outside the allowed range [0.1 m, 0.20 m], except if that is the vertical distance at the TE between cell 2 and 3.")
        print("Booms with max distance are at:")
        print(f"  Point a: x = {max_pair[0][0]:.3f}, z = {max_pair[0][1]:.3f}")
        print(f"  Point b: x = {max_pair[1][0]:.3f}, z = {max_pair[1][1]:.3f}")
    else:
            print("✅ All boom distances are within the required range.")
    '''
    
    # ---- Coordinate transformation to centroid-based system (for plotting) ----
    x_c_booms_centroid = [x - x_c_centroid for x in x_c_booms]
    z_c_booms_centroid = [z - z_c_centroid for z in z_c_booms]
    # Multiply by the chord to remove the normalization
    x_booms_centroid =  np.array(x_c_booms_centroid) * chord_length_arrondi
    z_booms_centroid =  np.array(z_c_booms_centroid) * chord_length_arrondi
    # Multiply by the chord to remove the normalization
    x_all = x_c_all*chord_length_arrondi
    z_all = z_c_all*chord_length_arrondi
    
    # ---- Coordinate transformation to centroid-based system (for the calculations) ----
    x_booms_ordered_centroid = [x - x_centroid for x in x_booms_ordered]
    z_booms_ordered_centroid = [z - z_centroid for z in z_booms_ordered]
    #print('')
    #print('x_booms_ordered_centroid =', x_booms_ordered_centroid)
    #print('')
    #print('z_booms_ordered_centroid =', z_booms_ordered_centroid)
    #print('')
    
    
    # ---- Delimitation of the two cells ----
    
    # -- Cell 1 --
    target_x_1 = closest_x_c_upper_1 * chord_length_arrondi # coord x de la fin de la cell 1   
    target_z_up_1 = z_c_upper_dotted_1 * chord_length_arrondi # coord z de la fin de la cell 1 upper side
    target_z_down_1 = z_c_lower_dotted_1 * chord_length_arrondi # coord z de la fin de la cell 1 lower side
      
    x_indice_cell_1 = np.argmin(np.abs(x_booms_ordered - target_x_1)) # indice which gives the position in the array of the x-coordinate of the boom concerned
    z_indice_cell_up_1 = np.argmin(np.abs(z_booms_ordered - target_z_up_1)) # same for the z-coordinate for the upper point (up)
    z_indice_cell_down_1 = np.argmin(np.abs(z_booms_ordered - target_z_down_1)) # same for the z-coordinate for the lower point (down)
    #print('x_indice_cell_1 = ', x_indice_cell_1)
    #print('z_indice_cell_up_1 = ', z_indice_cell_up_1)
    #print('z_indice_cell_down_1 = ', z_indice_cell_down_1)
    
    x_booms_cell_1 = x_booms_ordered_centroid[:x_indice_cell_1 + 1] + x_booms_ordered_centroid[z_indice_cell_down_1:] # Select all the x-coordinates of the booms in cell 1 and put them in one array
    z_booms_cell_1 = z_booms_ordered_centroid[:x_indice_cell_1 + 1] + z_booms_ordered_centroid[z_indice_cell_down_1:] # Select all the y-coordinates of the booms in cell 1 and put them in one array
    #print('x_booms_cell_1 = ', x_booms_cell_1)
    #print('z_booms_cell_1 = ', z_booms_cell_1)
    
    # -- Cell 2 --
    x_booms_cell_2 = x_booms_ordered_centroid[x_indice_cell_1:z_indice_cell_down_1 + 1] # Select all the x-coordinates of the booms in cell 2 and put them in one array
    z_booms_cell_2 = z_booms_ordered_centroid[x_indice_cell_1:z_indice_cell_down_1 + 1] # Select all the y-coordinates of the booms in cell 2 and put them in one array
    #print('x_booms_cell_2 = ', x_booms_cell_2)
    #print('z_booms_cell_2 = ', z_booms_cell_2)
    
    # -- NB : as there is only two cells, we only need the coordinates the delimiation between the two (if more, we will also need the coordinates of the limitation bertween the cell 2 and 3, etc) --

    
    # ---- Boom area ---- (different pts of the envelope) 
    
    # Point A
    T_x_A = structural_loads_regular_wing(n[0], alpha[0], Lwt[0], Dwt[0], Wwt, Mw[0])[0]
    T_y_A = structural_loads_regular_wing(n[0], alpha[0], Lwt[0], Dwt[0], Wwt, Mw[0])[1]
    T_z_A = structural_loads_regular_wing(n[0], alpha[0], Lwt[0], Dwt[0], Wwt, Mw[0])[2]
    M_x_A = structural_loads_regular_wing(n[0], alpha[0], Lwt[0], Dwt[0], Wwt, Mw[0])[3]
    M_y_A = structural_loads_regular_wing(n[0], alpha[0], Lwt[0], Dwt[0], Wwt, Mw[0])[4]
    M_z_A = structural_loads_regular_wing(n[0], alpha[0], Lwt[0], Dwt[0], Wwt, Mw[0])[5]
    B_pt_A = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_A, M_z_A,sigma_y_0, safety_factor)[0]
    sigma_yy_A = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_A, M_z_A,sigma_y_0, safety_factor)[1]
    
    # Point B
    T_x_B = structural_loads_regular_wing(n[1], alpha[1], Lwt[1], Dwt[1], Wwt, Mw[1])[0]
    T_y_B = structural_loads_regular_wing(n[1], alpha[1], Lwt[1], Dwt[1], Wwt, Mw[1])[1]
    T_z_B = structural_loads_regular_wing(n[1], alpha[1], Lwt[1], Dwt[1], Wwt, Mw[1])[2]
    M_x_B = structural_loads_regular_wing(n[1], alpha[1], Lwt[1], Dwt[1], Wwt, Mw[1])[3]
    M_y_B = structural_loads_regular_wing(n[1], alpha[1], Lwt[1], Dwt[1], Wwt, Mw[1])[4]
    M_z_B = structural_loads_regular_wing(n[1], alpha[1], Lwt[1], Dwt[1], Wwt, Mw[1])[5]
    B_pt_B = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_B, M_z_B,sigma_y_0, safety_factor)[0]
    sigma_yy_B = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_B, M_z_B,sigma_y_0, safety_factor)[1]
    
    # Point C
    T_x_C = structural_loads_regular_wing(n[2], alpha[2], Lwt[2], Dwt[2], Wwt, Mw[2])[0]
    T_y_C = structural_loads_regular_wing(n[2], alpha[2], Lwt[2], Dwt[2], Wwt, Mw[2])[1]
    T_z_C = structural_loads_regular_wing(n[2], alpha[2], Lwt[2], Dwt[2], Wwt, Mw[2])[2]
    M_x_C = structural_loads_regular_wing(n[2], alpha[2], Lwt[2], Dwt[2], Wwt, Mw[2])[3]
    M_y_C = structural_loads_regular_wing(n[2], alpha[2], Lwt[2], Dwt[2], Wwt, Mw[2])[4]
    M_z_C = structural_loads_regular_wing(n[2], alpha[2], Lwt[2], Dwt[2], Wwt, Mw[2])[5]
    B_pt_C = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_C, M_z_C,sigma_y_0, safety_factor)[0]
    sigma_yy_C = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_C, M_z_C,sigma_y_0, safety_factor)[1]
    
    # Point D
    T_x_D = structural_loads_regular_wing(n[3], alpha[3], Lwt[3], Dwt[3], Wwt, Mw[3])[0]
    T_y_D = structural_loads_regular_wing(n[3], alpha[3], Lwt[3], Dwt[3], Wwt, Mw[3])[1]
    T_z_D = structural_loads_regular_wing(n[3], alpha[3], Lwt[3], Dwt[3], Wwt, Mw[3])[2]
    M_x_D = structural_loads_regular_wing(n[3], alpha[3], Lwt[3], Dwt[3], Wwt, Mw[3])[3]
    M_y_D = structural_loads_regular_wing(n[3], alpha[3], Lwt[3], Dwt[3], Wwt, Mw[3])[4]
    M_z_D = structural_loads_regular_wing(n[3], alpha[3], Lwt[3], Dwt[3], Wwt, Mw[3])[5]
    B_pt_D = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_D, M_z_D,sigma_y_0, safety_factor)[0]
    sigma_yy_D = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_D, M_z_D,sigma_y_0, safety_factor)[1]
    
    # Point E
    T_x_E = structural_loads_regular_wing(n[4], alpha[4], Lwt[4], Dwt[4], Wwt, Mw[4])[0]
    T_y_E = structural_loads_regular_wing(n[4], alpha[4], Lwt[4], Dwt[4], Wwt, Mw[4])[1]
    T_z_E = structural_loads_regular_wing(n[4], alpha[4], Lwt[4], Dwt[4], Wwt, Mw[4])[2]
    M_x_E = structural_loads_regular_wing(n[4], alpha[4], Lwt[4], Dwt[4], Wwt, Mw[4])[3]
    M_y_E = structural_loads_regular_wing(n[4], alpha[4], Lwt[4], Dwt[4], Wwt, Mw[4])[4]
    M_z_E = structural_loads_regular_wing(n[4], alpha[4], Lwt[4], Dwt[4], Wwt, Mw[4])[5]
    B_pt_E = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_E, M_z_E,sigma_y_0, safety_factor)[0]
    sigma_yy_E = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_E, M_z_E,sigma_y_0, safety_factor)[1]

    # Max value (most critical)
    B_values = [B_pt_A, B_pt_B, B_pt_C, B_pt_D, B_pt_E]
    sigma_yy_values = [sigma_yy_A, sigma_yy_B, sigma_yy_C, sigma_yy_D, sigma_yy_E]

    B = max(B_values)
    index_of_max = B_values.index(B)
    sigma_yy = sigma_yy_values[index_of_max]
    sigma_yy = np.array(sigma_yy)
      
    #print(f"Index de la contrainte max : {index_of_max}")
    #print(f"sigma_yy associée (à B max) : {sigma_yy}")
    
    # Constants for conversion
    m2_to_mm2 = 1e6           # 1 m² = 1,000,000 mm²
    m2_to_in2 = 1550.0031     # 1 m² ≈ 1550.0031 in²
    print("")
    print(f"Boom area : {B:.9f} m² | {B * m2_to_mm2:.3f} mm² | {B * m2_to_in2:.3f} in²")
    print("")
    # print('The associated sigma_yy are :', sigma_yy) # should be an array for each boom
    
    
    # ---- Skin thickness ---- 
    
    # Taper effect along x
    x_pos_booms_tip = np.array(x_c_booms_ordered) * chord_tip + x_dist_root_tip 
    x_pos_booms_arrondi = np.array(x_c_booms_ordered) * chord_length_arrondi # where chord_length_arrondi is the chord from the tip, just before the 'arrondi'
    delta_x = x_pos_booms_tip - x_pos_booms_arrondi 
    delta_x = np.array(delta_x)
    #print('x_pos_booms_tip :', x_pos_booms_tip)
    #print('')
    #print('x_pos_booms_arrondi :', x_pos_booms_arrondi)
    #print('')
    #print('delta_x :', delta_x) 
    #print('')
    
    # Taper effect along z
    z_pos_booms_tip = np.array(z_c_booms_ordered) * chord_tip 
    z_pos_booms_arrondi = np.array(z_c_booms_ordered) * chord_length_arrondi # where chord_length_arrondi is the chord from the tip, just before the 'arrondi'
    delta_z = z_pos_booms_tip - z_pos_booms_arrondi 
    delta_z = np.array(delta_z)
    #print('z_pos_booms_tip :', z_pos_booms_tip)
    #print('')
    #print('z_pos_booms_arrondi :', z_pos_booms_arrondi)
    #print('')
    #print('delta_z :', delta_z)
    #print('')
    
    # Point A
    thickness_pt_A = skin_thickness(B, sigma_yy, delta_y, delta_x, delta_z, T_x_A, T_z_A, M_y_A, x_booms_ordered_centroid, z_booms_ordered_centroid, x_booms_cell_1, z_booms_cell_1, tau_max, x_booms_cell_2, z_booms_cell_2, x_centroid, z_centroid)
    
    # Point B
    thickness_pt_B = skin_thickness(B, sigma_yy, delta_y, delta_x, delta_z, T_x_B, T_z_B, M_y_B, x_booms_ordered_centroid, z_booms_ordered_centroid, x_booms_cell_1, z_booms_cell_1, tau_max, x_booms_cell_2, z_booms_cell_2, x_centroid, z_centroid)
    
    # Point C
    thickness_pt_C = skin_thickness(B, sigma_yy, delta_y, delta_x, delta_z, T_x_C, T_z_C, M_y_C, x_booms_ordered_centroid, z_booms_ordered_centroid, x_booms_cell_1, z_booms_cell_1, tau_max, x_booms_cell_2, z_booms_cell_2, x_centroid, z_centroid)
    
    # Point D
    thickness_pt_D = skin_thickness(B, sigma_yy, delta_y, delta_x, delta_z, T_x_D, T_z_D, M_y_D, x_booms_ordered_centroid, z_booms_ordered_centroid, x_booms_cell_1, z_booms_cell_1, tau_max, x_booms_cell_2, z_booms_cell_2, x_centroid, z_centroid)
    
    # Point E
    thickness_pt_E = skin_thickness(B, sigma_yy, delta_y, delta_x, delta_z, T_x_E, T_z_E, M_y_E, x_booms_ordered_centroid, z_booms_ordered_centroid, x_booms_cell_1, z_booms_cell_1, tau_max, x_booms_cell_2, z_booms_cell_2, x_centroid, z_centroid)
    
    # Take the maximum value
    thickness_values = [thickness_pt_A, thickness_pt_B, thickness_pt_C, thickness_pt_D, thickness_pt_E]
    thickness = max(thickness_values)

    # Conversion constants
    m_to_mm = 1000           # 1 m = 1000 mm
    m_to_in = 39.3701        # 1 m ≈ 39.3701 in
    print('')
    print(f"Skin thickness : {thickness:.9f} m | {thickness * m_to_mm:.6f} mm | {thickness * m_to_in:.6f} in")
    print('')
    
    
    # ---- Plotting ----
    if not plot_airfoil:
        return
    fig, ax = plt.subplots(figsize=(11, 5),dpi=300)
    airfoil_line, = ax.plot(x_c_all, z_c_all)  # Airfoil plot
    ax.scatter(x_c_booms, z_c_booms, color='red', zorder=3)

    # Get the color of the airfoil plot to use for the dotted line
    airfoil_color = airfoil_line.get_color()    

    # Plot the dotted line between the closest points to x/c = 0.25
    ax.plot([closest_x_c_upper_1, closest_x_c_upper_1], [z_c_upper_dotted_1, z_c_lower_dotted_1], '--', color=airfoil_color)

    # Plot the dotted line between the closest points to x/c = 0.8
    ax.plot([closest_x_c_upper_2, closest_x_c_upper_2], [z_c_upper_dotted_2, z_c_lower_dotted_2], '--', color=airfoil_color)
    
    # Plot the centroid location as a blue point
    ax.scatter(x_c_centroid, z_c_centroid, color='blue', zorder=4)
    
    # Add "C" label next to the centroid
    ax.text(x_c_centroid - 0.05, z_c_centroid - 0.025, 'C', color='blue', fontsize=12)
    
    # Plot coordinate axes at the centroid with arrows
    ax.arrow(x_c_centroid, z_c_centroid, 0.05, 0, head_width=0.02, head_length=0.02, fc='blue', ec='blue')  # X-axis arrow
    ax.text(x_c_centroid + 0.08, z_c_centroid - 0.01, 'x', color='blue', fontsize=12) # X-axis name
    ax.arrow(x_c_centroid, z_c_centroid, 0, 0.05, head_width=0.02, head_length=0.02, fc='blue', ec='blue')  # Z-axis arrow
    ax.text(x_c_centroid - 0.01, z_c_centroid + 0.08, 'z', color='blue', fontsize=12) # Z-axis name
    
    # Plot max distance segment in orange 
    ax.plot([max_pair[0][0]/chord_length_arrondi, max_pair[1][0]/chord_length_arrondi],
            [max_pair[0][1]/chord_length_arrondi, max_pair[1][1]/chord_length_arrondi],
            linestyle='--', color='orange', linewidth=2)
        
    # Set equal scaling for both axes (same units for x and z axes)
    ax.set_aspect('equal')
    
    # Labels and title
    ax.set_xlabel('$x/c$')
    ax.set_ylabel('$z/c$')
    ax.set_title(f'SC(2)-0710 Airfoil with {n_booms} booms')
    ax.grid()
    plt.show()
    #plt.savefig(fname='airfoil_x_c_centroid_location.pdf')
    
    # ---- Plotting (Real Size) ----
    fig2, ax2 = plt.subplots(figsize=(11, 5), dpi=300)
    ax2.plot(x_all, z_all, label='Airfoil')  # Airfoil in meters
    ax2.scatter(x_booms, z_booms, color='red', zorder=3)

    # Plot the dotted lines in meters
    x_dotted_1 = closest_x_c_upper_1 * chord_length_arrondi
    z_dotted_1_top = z_c_upper_dotted_1 * chord_length_arrondi
    z_dotted_1_bot = z_c_lower_dotted_1 * chord_length_arrondi

    x_dotted_2 = closest_x_c_upper_2 * chord_length_arrondi
    z_dotted_2_top = z_c_upper_dotted_2 * chord_length_arrondi
    z_dotted_2_bot = z_c_lower_dotted_2 * chord_length_arrondi

    ax2.plot([x_dotted_1, x_dotted_1], [z_dotted_1_top, z_dotted_1_bot], '--', color=airfoil_color)
    ax2.plot([x_dotted_2, x_dotted_2], [z_dotted_2_top, z_dotted_2_bot], '--', color=airfoil_color)

    # Plot centroid in real coordinates
    ax2.scatter(x_centroid, z_centroid, color='blue', zorder=4)
    ax2.text(x_centroid - 0.25, z_centroid - 0.125, 'C', color='blue', fontsize=12)

    # Plot coordinate axes at the centroid
    ax2.arrow(x_centroid, z_centroid, 0.05 * chord_length_arrondi, 0, head_width=0.02 * chord_length_arrondi, head_length=0.02 * chord_length_arrondi, fc='blue', ec='blue')
    ax2.text(x_centroid + 0.4, z_centroid - 0.05, 'x', color='blue', fontsize=12)
    ax2.arrow(x_centroid, z_centroid, 0, 0.05 * chord_length_arrondi, head_width=0.02 * chord_length_arrondi, head_length=0.02 * chord_length_arrondi, fc='blue', ec='blue')
    ax2.text(x_centroid - 0.05 , z_centroid + 0.4, 'z', color='blue', fontsize=12)

    ax2.set_aspect('equal')
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel('z [m]')
    ax2.set_title(f'SC(2)-0710 Airfoil in Real Scale ({n_booms} booms)')
    ax2.grid()
    #ax2.legend()
    plt.show()

# -------------------------------------------------------------------

plotAirfoil(False, n_booms=37, x_c_max=0.7, x_c_cell_1=0.25) 
