import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import math


# Geometry wings
beta = np.radians(4.3) # Setting angle [rad] 
(dx_ac,dy_ac,dz_ac) = (10.140 - (4.5+3.109), 10.196 - 7.27,0) # Coord AC of the wing ()
(dx_w_cg,dy_w_cg,dz_w_cg) = (2.8471729, 3.1292394, 0.1595797) #lacement of the center of gravity of the wing to the root of the regular wing



# ---- Structural loads ---- (will be imported from another fct, here just expl values)

# Here, import all the n and all the alpha from the envelope : 
n =  [2.5, 2.5, 0, -1, -1] # Load factor
alpha = [ 0.06151096, -0.04670625, -0.10033783, -0.13220202, -0.24113608]   # [rad]
L_overall = [1057331.80668917, 1095922.26154146, 79107.31234472, -336009.15840703, -347684.38494123]   # [N]
D_wing = [ 5313.09423855, 15816.89997096, 15816.89997096, 10122.81598142, 2123.01255461] # [N]
M_wing =  [-303561.349929, -903691.7644038, -903691.7644038, -578362.72921843, -121297.40751009]  # [N.m]

# Trouver un moyen de calculer L_wing via la formule analytique avec C_L = 0.433 via pondération des surfaces 
# Drag : pondérer par les surfaces aussi 

Surf_aile_equiv = 119.85 #[m²]
Surf_demi_aile_equiv = Surf_aile_equiv/2 #[m²]
Surf_demi_wing_studied = 21.86460271713637 # [m²]
Surf_demi_trans_wing = 13.78039728282926 # [m²]

Coeff = (2*Surf_demi_wing_studied)/Surf_aile_equiv
Coeff_2 = (Surf_demi_trans_wing)/ (Surf_demi_wing_studied + Surf_demi_trans_wing)
print(Coeff)
print(Coeff_2)



Lwt = np.array(L_overall)*Coeff # Extract the correct lift
Dwt = np.array(D_wing)*Coeff_2
Wwt = (4735.9488633)*9.81*Coeff_2 # [kg] mass -> to force [N] 
Mw = np.array(M_wing)*Coeff_2


# ---- Material ----

# sigma_x_0 = 702*10**6 # CHOOSE THE MATERIAL
# tau_max = 373.5*10**6 # maximum shear stress
# safety_factor = 1.5

# Alu moyen Granta
sigma_x_0 = (797+683)/2*10**6 # CHOOSE THE MATERIAL
tau_max = (402+345)/2*10**6 # maximum shear stress
safety_factor = 1.5

# CFRP
# sigma_x_0 = 1500*10**6
# tau_max = 94*10**6
# safety_factor = 1.5

print('')
print('Lift 1 wing = ', Lwt/2)
print('')
print('Drag 1 wing =', Dwt/2)
print('')
print('Moment 1 wing =', Mw/2)
print('')
print('Weight (N) =', Wwt/2)
print('')

# -------------------------------------------------------------------

def structural_loads_regular_wing (n, alpha, Lwt, Dwt, Wwt, Mw): # There will be more parameters as the lift also varies, etc
    
    T_x = 0
    T_y = ((n*Wwt/2-Lwt/2)*np.sin(alpha+beta) + Dwt/2*np.cos(alpha+beta))
    T_z = ((-n*Wwt/2+Lwt/2)*np.cos(alpha+beta) + Dwt/2*np.sin(alpha+beta))
    
    M_x = np.cos(alpha+beta)/2 * (-n*Wwt* dy_w_cg + Lwt * dy_ac - Dwt * dz_ac) + np.sin(alpha+beta)/2 * (-n*Wwt * dz_w_cg + Lwt * dz_ac + Dwt * dy_ac) - Mw/2
    M_y = np.cos(alpha+beta)/2 * (n*Wwt * dx_w_cg - Lwt * dx_ac) + np.sin(alpha+beta)/2 * (Dwt * dx_ac) 
    M_z = np.sin(alpha+beta)/2 * (n*Wwt * dx_w_cg - Lwt * dx_ac) + np.cos(alpha+beta)/2 * (Dwt * dx_ac)  

    # T_x = -T_x
    # T_y = -T_y
    # T_z = -T_z

    # M_x = -M_x
    # M_y = -M_y
    # M_z = -M_z

    
    # Old version    
    # T_x = (n*Wwt/2-Lwt/2)*np.sin(alpha+beta) + Dwt/2*np.cos(alpha+beta) #alpha is the angle of attack, beta is the setting angle
    # T_y = 0
    # T_z = (-n*Wwt/2+Lwt/2)*np.cos(alpha+beta) + Dwt/2*np.sin(alpha+beta)
    # M_x = 1/2*(-n*(Wwt*dy_w_cg)+Lwt*dy_ac)*np.cos(alpha+beta) + Dwt/2*dy_ac*np.sin(alpha+beta)#dyw is the y placement of the center of gravity to the root of the regular wing , dya is the distance of the aerodynamic center to the root wing
    # M_y = 1/2*(-n*Wwt*dx_w_cg + Lwt*dx_ac -Dwt*dz_ac)*np.cos(alpha+beta) + 1/2*(-n*Wwt*dz_w_cg +Lwt*dz_ac + Dwt*dx_ac)*np.sin(alpha+beta) + Mw/2 
    # M_z = 1/2*(+n*(Wwt*dy_w_cg)-Lwt*dy_ac)*np.sin(alpha+beta) + Dwt/2*dy_ac*np.cos(alpha+beta)
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

def skin_thickness(B, sigma_xx, delta_x, delta_y, delta_z, T_y, T_z, M_x, y_booms_ordered_centroid, z_booms_ordered_centroid, y_booms_cell_1, z_booms_cell_1, tau_max, y_booms_cell_2, z_booms_cell_2, y_centroid, z_centroid) :
    
    
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

    print('denom :', denom)    
    print('factor_1 :', factor_1)
    print('factor_2 :', factor_2)
    
    # -- Cell 1 --
    q_0_cell_1 = [0] # Cut in cell 1, the open shear flow between booms 7 and 6 is zero                                              (q_0_7,6)
    q_0_cell_1.append(q_0_cell_1[0] - factor_1 * z_booms_ordered_centroid[6-1]*B - factor_2 * y_booms_ordered_centroid[6-1]*B)     # (q_0_6,5)
    q_0_cell_1.append(q_0_cell_1[1] - factor_1 * z_booms_ordered_centroid[5-1]*B - factor_2 * y_booms_ordered_centroid[5-1]*B)     # (q_0_5,4)
    q_0_cell_1.append(q_0_cell_1[2] - factor_1 * z_booms_ordered_centroid[4-1]*B - factor_2 * y_booms_ordered_centroid[4-1]*B)     # (q_0_4,3)
    q_0_cell_1.append(q_0_cell_1[3] - factor_1 * z_booms_ordered_centroid[3-1]*B - factor_2 * y_booms_ordered_centroid[3-1]*B)     # (q_0_3,2)
    q_0_cell_1.append(q_0_cell_1[4] - factor_1 * z_booms_ordered_centroid[2-1]*B - factor_2 * y_booms_ordered_centroid[2-1]*B)     # (q_0_2,1)
    q_0_cell_1.append(q_0_cell_1[5] - factor_1 * z_booms_ordered_centroid[1-1]*B - factor_2 * y_booms_ordered_centroid[1-1]*B)     # (q_0_1,37)
    q_0_cell_1.append(q_0_cell_1[6] - factor_1 * z_booms_ordered_centroid[37-1]*B - factor_2 * y_booms_ordered_centroid[37-1]*B)   # (q_0_37,36)
    q_0_cell_1.append(q_0_cell_1[7] - factor_1 * z_booms_ordered_centroid[36-1]*B - factor_2 * y_booms_ordered_centroid[36-1]*B)   # (q_0_36,35)
    q_0_cell_1.append(q_0_cell_1[8] - factor_1 * z_booms_ordered_centroid[35-1]*B - factor_2 * y_booms_ordered_centroid[35-1]*B)   # (q_0_35,34)
    q_0_cell_1.append(q_0_cell_1[9] - factor_1 * z_booms_ordered_centroid[34-1]*B - factor_2 * y_booms_ordered_centroid[34-1]*B)   # (q_0_34,33)
    q_0_cell_1.append(q_0_cell_1[10] - factor_1 * z_booms_ordered_centroid[33-1]*B - factor_2 * y_booms_ordered_centroid[33-1]*B)  # (q_0_33,32) 
    q_spar_1_2 = - factor_1 * z_booms_ordered_centroid[7-1]*B - factor_2 * y_booms_ordered_centroid[7-1]*B  # (q_0_7,32) et pas (q_0_32,7) !
    q_0_cell_1.append(-q_spar_1_2)  # (q_0_32,7) mettre dans le meme sens que les autres
    
    # Mettre le tableau dans le meme ordre que y_booms_cell_1 et z_booms_cell_1
    n = 6  # nombre d'éléments à déplacer à la fin 
    #print('q_0_cell_1 :', q_0_cell_1)
    q_0_cell_1 = q_0_cell_1[n:] + q_0_cell_1[:n] # Cut the array between case n=6 and 7 and invert the 2 part
    # Invert the order of all elements
    q_0_cell_1 = q_0_cell_1[::-1] 
    # The array "q_0_cell_1_t" is now in the same order than y_booms_cell_1 and z_booms_cell_1
    # print('q_0_cell_1 :', q_0_cell_1)
    # print('max =', np.max(np.abs(q_0_cell_1)))
    
    # -- Cell 2 (t) -- 
    q_0_cell_2 = [0]                                                                                                              # (q_0_7,8)
    q_0_cell_2.append(q_0_cell_2[0] - factor_1 * z_booms_ordered_centroid[8-1]*B - factor_2 * y_booms_ordered_centroid[8-1]*B)    # (q_0_8,9)
    q_0_cell_2.append(q_0_cell_2[1] - factor_1 * z_booms_ordered_centroid[9-1]*B - factor_2 * y_booms_ordered_centroid[9-1]*B)    # (q_0_9,10)
    q_0_cell_2.append(q_0_cell_2[2] - factor_1 * z_booms_ordered_centroid[10-1]*B - factor_2 * y_booms_ordered_centroid[10-1]*B)  # (q_0_10,11)
    q_0_cell_2.append(q_0_cell_2[3] - factor_1 * z_booms_ordered_centroid[11-1]*B - factor_2 * y_booms_ordered_centroid[11-1]*B)  # (q_0_11,12)
    q_0_cell_2.append(q_0_cell_2[4] - factor_1 * z_booms_ordered_centroid[12-1]*B - factor_2 * y_booms_ordered_centroid[12-1]*B)  # (q_0_12,13)
    q_0_cell_2.append(q_0_cell_2[5] - factor_1 * z_booms_ordered_centroid[13-1]*B - factor_2 * y_booms_ordered_centroid[13-1]*B)  # (q_0_13,14)
    q_0_cell_2.append(q_0_cell_2[6] - factor_1 * z_booms_ordered_centroid[14-1]*B - factor_2 * y_booms_ordered_centroid[14-1]*B)  # (q_0_14,15)
    q_0_cell_2.append(q_0_cell_2[7] - factor_1 * z_booms_ordered_centroid[15-1]*B - factor_2 * y_booms_ordered_centroid[15-1]*B)  # (q_0_15,16)
    q_0_cell_2.append(q_0_cell_2[8] - factor_1 * z_booms_ordered_centroid[16-1]*B - factor_2 * y_booms_ordered_centroid[16-1]*B)  # (q_0_16,17)
    q_0_cell_2.append(q_0_cell_2[9] - factor_1 * z_booms_ordered_centroid[17-1]*B - factor_2 * y_booms_ordered_centroid[17-1]*B)  # (q_0_17,18)
    q_0_cell_2.append(q_0_cell_2[10] - factor_1 * z_booms_ordered_centroid[18-1]*B - factor_2 * y_booms_ordered_centroid[18-1]*B) # (q_0_18,19)
    q_0_cell_2.append(q_0_cell_2[11] - factor_1 * z_booms_ordered_centroid[19-1]*B - factor_2 * y_booms_ordered_centroid[19-1]*B) # (q_0_19,20)
    q_0_cell_2.append(q_0_cell_2[12] - factor_1 * z_booms_ordered_centroid[20-1]*B - factor_2 * y_booms_ordered_centroid[20-1]*B) # (q_0_20,21)
    q_0_cell_2.append(q_0_cell_2[13] - factor_1 * z_booms_ordered_centroid[21-1]*B - factor_2 * y_booms_ordered_centroid[21-1]*B) # (q_0_21,22)
    q_0_cell_2.append(q_0_cell_2[14] - factor_1 * z_booms_ordered_centroid[22-1]*B - factor_2 * y_booms_ordered_centroid[22-1]*B) # (q_0_22,23)
    q_0_cell_2.append(q_0_cell_2[15] - factor_1 * z_booms_ordered_centroid[23-1]*B - factor_2 * y_booms_ordered_centroid[23-1]*B) # (q_0_23,24)
    q_0_cell_2.append(q_0_cell_2[16] - factor_1 * z_booms_ordered_centroid[24-1]*B - factor_2 * y_booms_ordered_centroid[24-1]*B) # (q_0_24,25)
    q_0_cell_2.append(q_0_cell_2[17] - factor_1 * z_booms_ordered_centroid[25-1]*B - factor_2 * y_booms_ordered_centroid[25-1]*B) # (q_0_25,26)
    q_0_cell_2.append(q_0_cell_2[18] - factor_1 * z_booms_ordered_centroid[26-1]*B - factor_2 * y_booms_ordered_centroid[26-1]*B) # (q_0_26,27)
    q_0_cell_2.append(q_0_cell_2[19] - factor_1 * z_booms_ordered_centroid[27-1]*B - factor_2 * y_booms_ordered_centroid[27-1]*B) # (q_0_27,28)
    q_0_cell_2.append(q_0_cell_2[20] - factor_1 * z_booms_ordered_centroid[28-1]*B - factor_2 * y_booms_ordered_centroid[28-1]*B) # (q_0_28,29)
    q_0_cell_2.append(q_0_cell_2[21] - factor_1 * z_booms_ordered_centroid[29-1]*B - factor_2 * y_booms_ordered_centroid[29-1]*B) # (q_0_29,30)
    q_0_cell_2.append(q_0_cell_2[22] - factor_1 * z_booms_ordered_centroid[30-1]*B - factor_2 * y_booms_ordered_centroid[30-1]*B) # (q_0_30,31)
    q_0_cell_2.append(q_0_cell_2[23] - factor_1 * z_booms_ordered_centroid[31-1]*B - factor_2 * y_booms_ordered_centroid[31-1]*B) # (q_0_31,32)   
    q_0_cell_2.append(-q_spar_1_2) # (q_0_32,7)
    
    # Invert the sign of all elements as I went clockwise for cell 2 
    q_0_cell_2 = [-y for y in q_0_cell_2]
    # print('q_0_cell_2 :', q_0_cell_2)
    # print('max =', np.max(np.abs(q_0_cell_2)))
    
    
    # ---- Lengths of the segments and the cells ----
    #lengths = dist_2_booms(y_booms_ordered_centroid, z_booms_ordered_centroid)[0]
    intersecting_length = dist_2_pts(y_booms_ordered_centroid[7-1], z_booms_ordered_centroid[7-1], y_booms_ordered_centroid[32-1], z_booms_ordered_centroid[32-1])
    #print('intersecting_length : ', intersecting_length)

    # Cell 1
    lengths_c_1 = dist_2_booms(y_booms_cell_1, z_booms_cell_1)[0]
    l_cell_1 = np.sum(lengths_c_1)
    # print('lengths_c_1 = ', lengths_c_1)
    # print(l_cell_1)
    
    # Cell 2
    lengths_c_2 = dist_2_booms(y_booms_cell_2, z_booms_cell_2)[0]
    l_cell_2 = np.sum(lengths_c_2)
    # print('lengths_c_2 = ', lengths_c_2)
    # print(l_cell_2)
    
    
    # ---- Cell Area ----
    A_h_c_1 = polygon_area(y_booms_cell_1, z_booms_cell_1)
    A_h_c_2 = polygon_area(y_booms_cell_2, z_booms_cell_2)
    #print(f"Cell Areas: A_h_c_1 = {A_h_c_1:.4f} m², A_h_c_2 = {A_h_c_2:.4f} m²")
    
    
    # ---- Open shear flow integration ----
    # Cell 1
    q_0_cell_1 = np.array(q_0_cell_1) # Convert it to NumPy arrays
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
    swept_area_cell_1 = swept_area_from_center(y_booms_cell_1, z_booms_cell_1, y_centroid = 0, z_centroid = 0)[1]
    # Cell 2
    swept_area_cell_2 = swept_area_from_center(y_booms_cell_2, z_booms_cell_2, y_centroid = 0, z_centroid = 0)[1]
    # Term in the momentum equilibrium : 
    term_swept_area = 2 * (np.sum(q_0_cell_1 * swept_area_cell_1) + np.sum(q_0_cell_2 * swept_area_cell_2))
    #print('term_swept_area =', term_swept_area)
    
    
    # ---- Correction term computation : solve a 2 eqns syst (because 2 cells) ----
    # Syst parameters : unknowns : "x" for q_1_corr and "y" for q_2_corr
    
    # Twist rate comaptibility equation
    coeff_x_eq1 = l_cell_1 * A_h_c_2 + intersecting_length * A_h_c_1
    coeff_y_eq1 = (-1)* intersecting_length * A_h_c_2 - l_cell_2 * A_h_c_2
    ind_term_eq1 = open_shear_flux_2 * A_h_c_1 - open_shear_flux_1 * A_h_c_2
    
    # Momentum equilibrium
    coeff_x_eq2 = 2 * A_h_c_1
    coeff_y_eq2 = 2 * A_h_c_2
    
    term_P_z = B * np.sum(y_booms_ordered_centroid * sigma_xx * delta_z/delta_x)
    term_P_y = B * np.sum(y_booms_ordered_centroid * sigma_xx * delta_y/delta_x)
    ind_term_eq2 = (y_centroid * T_z - z_centroid * T_y) - term_swept_area - term_P_z + term_P_y 
    
    # Solve syst : 
    A = np.array([[coeff_x_eq1, coeff_y_eq1],[coeff_x_eq2, coeff_y_eq2]]) # Coeff matrix
    B = np.array([ind_term_eq1, ind_term_eq2]) # Independent term vector 
    sol = np.linalg.solve(A, B)
    q_1_corr, q_2_corr = sol
    print(f"q_1_corr = {q_1_corr:.4f}, q_2_corr = {q_2_corr:.4f}")
    
    # ---- Shear flow due to torsion ----
    # Coeff of the eq of the syst :
    a_1 = 2 * A_h_c_1
    b_1 = 2 * A_h_c_2
    c_1 = M_x
    a_2 = (l_cell_1/A_h_c_1) + (intersecting_length/A_h_c_2)
    b_2 = -((l_cell_2/A_h_c_2) + (intersecting_length/A_h_c_1))
    c_2 = 0
    
    # Solve syst : 
    A_torsion = np.array([[a_1, b_1],[a_2, b_2]]) # Coeff matrix
    B_torsion = np.array([c_1, c_2]) # Independent term vector 
    sol_torsion = np.linalg.solve(A_torsion, B_torsion)
    q_1_torsion, q_2_torsion = sol_torsion
    print(f"q_1_torsion = {q_1_torsion:.4f}, q_2_torsion = {q_2_torsion:.4f}")
    
    # ---- Shear flow (closed) ----
    q_closed_cell_1 = np.array(q_0_cell_1) + q_1_corr + q_1_torsion # Apply the correction and the shear flow due to torsion
    q_closed_cell_2 = np.array(q_0_cell_2) + q_2_corr + q_2_torsion # Apply the correction and the shear flow due to torsion
    
    # Take into account the spar which has the correction added in both sense
    q_closed_spar = q_closed_cell_1[6] - q_2_corr - q_2_torsion # The correction of cell 2 are in the other direction
    # Now put it back
    q_closed_cell_1[6] = q_closed_spar
    q_closed_cell_2[25] = -q_closed_spar
    
    q_closed = np.concatenate([q_closed_cell_1, q_closed_cell_2]) # Note that q_closed contains the first spar shear flow twice but it is not a problem
    #print('q_closed = ',q_closed)
    # ---- Thickness computation ----   
    print("q_max =", np.max(np.abs(q_closed)) )
    print("")
    #print('tau_max', tau_max)
    thickness = np.max(np.abs(q_closed)) / (tau_max/safety_factor) 
    
    return thickness

# -------------------------------------------------------------------

def plotAirfoil(plot_airfoil, n_booms, y_c_max, y_c_cell_1):
    
    # ---- Geometry ----
    wing_semi_span = 10 #[m]
    x_tip = wing_semi_span
    x_root = 0
    x_arrondi = 3.10921985816 #
    chord_root = 4.779 # 
    chord_tip = 2.35
    chord_length_arrondi = chord_tip + (chord_root - chord_tip)/(x_root - x_tip) * (x_arrondi - x_tip)
    #print('chord_length :',chord_length_arrondi)
    sweep_angle = np.radians(31.599) #[°]
    delta_x = (wing_semi_span - x_arrondi)
    y_dist_root_tip = delta_x *np.tan(sweep_angle)
    print('delta_x : ', delta_x)
    #print('y_dist_root_tip : ', y_dist_root_tip)
    #print('')
    
    # Load airfoil data
    data_wing = pd.read_csv("sc20710_XYZ.csv")
    y_c_all = data_wing.iloc[:, 0].values
    z_c_all = data_wing.iloc[:, 1].values  

    # Separate and sort upper and lower surfaces
    upper_surface = data_wing[data_wing.iloc[:, 1] >= 0].copy().sort_values(by=data_wing.columns[0])  # z >= 0
    lower_surface = data_wing[data_wing.iloc[:, 1] < 0].copy().sort_values(by=data_wing.columns[0])  # z < 0

    # Add the leading edge boom at (0, 0)
    y_c_booms = [0.0]
    z_c_booms = [0.0]

    # Ensure odd number of booms for symmetry around leading edge
    if (n_booms - 1) % 2 != 0:
        raise ValueError("Number of booms must be odd to have (0,0) + symmetric pairs.")

    n_pairs = (n_booms - 1) // 2


    # Generate y positions from 0 to y_c_max, skipping the leading edge (already added)
    y_c_positions = np.linspace(0, y_c_max, n_pairs + 1)[1:]

    # Interpolate z from upper and lower surfaces
    z_c_upper = np.interp(y_c_positions, upper_surface.iloc[:, 0], upper_surface.iloc[:, 1])
    z_c_lower = np.interp(y_c_positions, lower_surface.iloc[:, 0], lower_surface.iloc[:, 1])

    # Add symmetric booms
    for y, z_u, z_l in zip(y_c_positions, z_c_upper, z_c_lower):
        y_c_booms.extend([y, y])
        z_c_booms.extend([z_u, z_l])


    # ---- Find the points closest to y/c = 0.25 for the dotted line ----
    y_c_dotted_target_1 = y_c_cell_1
    closest_idy_upper_1 = np.argmin(np.abs(y_c_positions - y_c_dotted_target_1))
    closest_y_c_upper_1 = y_c_positions[closest_idy_upper_1]
    z_c_upper_dotted_1 = np.interp(closest_y_c_upper_1, upper_surface.iloc[:, 0], upper_surface.iloc[:, 1])
    z_c_lower_dotted_1 = np.interp(closest_y_c_upper_1, lower_surface.iloc[:, 0], lower_surface.iloc[:, 1])


    # ---- Find the points closest to y/c = y_c_max for the dotted line ----
    y_c_dotted_target_2 = y_c_max
    closest_idy_upper_2 = np.argmin(np.abs(y_c_positions - y_c_dotted_target_2))
    closest_y_c_upper_2 = y_c_positions[closest_idy_upper_2]
    z_c_upper_dotted_2 = np.interp(closest_y_c_upper_2, upper_surface.iloc[:, 0], upper_surface.iloc[:, 1])
    z_c_lower_dotted_2 = np.interp(closest_y_c_upper_2, lower_surface.iloc[:, 0], lower_surface.iloc[:, 1])


    # ---- Centroid location ----
    y_c_centroid = np.sum(y_c_booms) / len(y_c_booms)
    z_c_centroid = np.sum(z_c_booms) / len(z_c_booms)
    # Multiply by the chord
    y_booms = np.array(y_c_booms) * chord_length_arrondi
    z_booms = np.array(z_c_booms) * chord_length_arrondi
    y_centroid = y_c_centroid * chord_length_arrondi
    z_centroid = z_c_centroid * chord_length_arrondi
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
    
    # Compute distances in correct order
    dist, max_dist, max_pair = dist_2_booms(y_booms_ordered, z_booms_ordered)
    
    # Reorder also the y_c_booms and z_c_booms in the same way
    # (normalised coordinates corresponding to reordered_booms)
    reordered_booms_c = [(y / chord_length_arrondi, z / chord_length_arrondi) for y, z in reordered_booms]

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
    # Multiply by the chord to remove the normalization
    y_booms_centroid =  np.array(y_c_booms_centroid) * chord_length_arrondi
    z_booms_centroid =  np.array(z_c_booms_centroid) * chord_length_arrondi
    # Multiply by the chord to remove the normalization
    y_all = y_c_all*chord_length_arrondi
    z_all = z_c_all*chord_length_arrondi
    
    # ---- Coordinate transformation to centroid-based system (for the calculations) ----
    y_booms_ordered_centroid = [y - y_centroid for y in y_booms_ordered]
    z_booms_ordered_centroid = [z - z_centroid for z in z_booms_ordered]
    #print('')
    #print('y_booms_ordered_centroid =', y_booms_ordered_centroid)
    #print('')
    #print('z_booms_ordered_centroid =', z_booms_ordered_centroid)
    #print('')
    
    
    # ---- Delimitation of the two cells ----
    
    # -- Cell 1 --
    target_y_1 = closest_y_c_upper_1 * chord_length_arrondi # coord y de la fin de la cell 1   
    target_z_up_1 = z_c_upper_dotted_1 * chord_length_arrondi # coord z de la fin de la cell 1 upper side
    target_z_down_1 = z_c_lower_dotted_1 * chord_length_arrondi # coord z de la fin de la cell 1 lower side
      
    y_indice_cell_1 = np.argmin(np.abs(y_booms_ordered - target_y_1)) # indice which gives the position in the array of the y-coordinate of the boom concerned
    z_indice_cell_up_1 = np.argmin(np.abs(z_booms_ordered - target_z_up_1)) # same for the z-coordinate for the upper point (up)
    z_indice_cell_down_1 = np.argmin(np.abs(z_booms_ordered - target_z_down_1)) # same for the z-coordinate for the lower point (down)
    #print('y_indice_cell_1 = ', y_indice_cell_1)
    #print('z_indice_cell_up_1 = ', z_indice_cell_up_1)
    #print('z_indice_cell_down_1 = ', z_indice_cell_down_1)
    
    y_booms_cell_1 = y_booms_ordered_centroid[:y_indice_cell_1 + 1] + y_booms_ordered_centroid[z_indice_cell_down_1:] # Select all the y-coordinates of the booms in cell 1 and put them in one array
    z_booms_cell_1 = z_booms_ordered_centroid[:y_indice_cell_1 + 1] + z_booms_ordered_centroid[z_indice_cell_down_1:] # Select all the y-coordinates of the booms in cell 1 and put them in one array
    #print('y_booms_cell_1 = ', y_booms_cell_1)
    #print('z_booms_cell_1 = ', z_booms_cell_1)
    
    # -- Cell 2 --
    y_booms_cell_2 = y_booms_ordered_centroid[y_indice_cell_1:z_indice_cell_down_1 + 1] # Select all the y-coordinates of the booms in cell 2 and put them in one array
    z_booms_cell_2 = z_booms_ordered_centroid[y_indice_cell_1:z_indice_cell_down_1 + 1] # Select all the y-coordinates of the booms in cell 2 and put them in one array
    #print('y_booms_cell_2 = ', y_booms_cell_2)
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
    #print('B_values =', B_values)
    sigma_xx_values = [sigma_xx_A, sigma_xx_B, sigma_xx_C, sigma_xx_D, sigma_xx_E]

    B = max(B_values)
    index_of_max = B_values.index(B)
    sigma_xx = sigma_xx_values[index_of_max]
    sigma_xx = np.array(sigma_xx)
      
    print(f"Index de la contrainte max : {index_of_max}")
    print(f"sigma_xx associée (à B max) : {sigma_xx}")
    
    # Constants for conversion
    m2_to_mm2 = 1e6           # 1 m² = 1,000,000 mm²
    m2_to_in2 = 1550.0031     # 1 m² ≈ 1550.0031 in²
    print("")
    print(f"Boom area : {B:.9f} m² | {B * m2_to_mm2:.3f} mm² | {B * m2_to_in2:.3f} in²")
    print("")
    # print('The associated sigma_xx are :', sigma_xx) # should be an array for each boom
    
    
    # ---- Skin thickness ---- 
    
    # Taper effect along y
    y_pos_booms_tip = np.array(y_c_booms_ordered) * chord_tip + y_dist_root_tip 
    y_pos_booms_arrondi = np.array(y_c_booms_ordered) * chord_length_arrondi # where chord_length_arrondi is the chord from the tip, just before the 'arrondi'
    delta_y = y_pos_booms_tip - y_pos_booms_arrondi 
    delta_y = np.array(delta_y)
    # print('y_pos_booms_tip :', y_pos_booms_tip)
    # print('')
    # print('y_pos_booms_arrondi :', y_pos_booms_arrondi)
    # print('')
    # print('delta_y :', delta_y) 
    # print('')
    
    # Taper effect along z
    z_pos_booms_tip = np.array(z_c_booms_ordered) * chord_tip 
    z_pos_booms_arrondi = np.array(z_c_booms_ordered) * chord_length_arrondi # where chord_length_arrondi is the chord from the tip, just before the 'arrondi'
    delta_z = z_pos_booms_tip - z_pos_booms_arrondi 
    delta_z = np.array(delta_z)
    # print('z_pos_booms_tip :', z_pos_booms_tip)
    # print('')
    # print('z_pos_booms_arrondi :', z_pos_booms_arrondi)
    # print('')
    # print('delta_z :', delta_z)
    # print('')
    
    # Point A
    thickness_pt_A = skin_thickness(B, sigma_xx, delta_x, delta_y, delta_z, T_y_A, T_z_A, M_x_A, y_booms_ordered_centroid, z_booms_ordered_centroid, y_booms_cell_1, z_booms_cell_1, tau_max, y_booms_cell_2, z_booms_cell_2, y_centroid, z_centroid)
    
    # Point B
    thickness_pt_B = skin_thickness(B, sigma_xx, delta_x, delta_y, delta_z, T_y_B, T_z_B, M_x_B, y_booms_ordered_centroid, z_booms_ordered_centroid, y_booms_cell_1, z_booms_cell_1, tau_max, y_booms_cell_2, z_booms_cell_2, y_centroid, z_centroid)
    
    # Point C
    thickness_pt_C = skin_thickness(B, sigma_xx, delta_x, delta_y, delta_z, T_y_C, T_z_C, M_x_C, y_booms_ordered_centroid, z_booms_ordered_centroid, y_booms_cell_1, z_booms_cell_1, tau_max, y_booms_cell_2, z_booms_cell_2, y_centroid, z_centroid)
    
    # Point D
    thickness_pt_D = skin_thickness(B, sigma_xx, delta_x, delta_y, delta_z, T_y_D, T_z_D, M_x_D, y_booms_ordered_centroid, z_booms_ordered_centroid, y_booms_cell_1, z_booms_cell_1, tau_max, y_booms_cell_2, z_booms_cell_2, y_centroid, z_centroid)
    
    # Point E
    thickness_pt_E = skin_thickness(B, sigma_xx, delta_x, delta_y, delta_z, T_y_E, T_z_E, M_x_E, y_booms_ordered_centroid, z_booms_ordered_centroid, y_booms_cell_1, z_booms_cell_1, tau_max, y_booms_cell_2, z_booms_cell_2, y_centroid, z_centroid)
    
    # Take the maximum value
    thickness_values = [thickness_pt_A, thickness_pt_B, thickness_pt_C, thickness_pt_D, thickness_pt_E]
    print('Thickness_values =', thickness_values)
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
    airfoil_line, = ax.plot(y_c_all, z_c_all)  # Airfoil plot
    ax.scatter(y_c_booms, z_c_booms, color='red', zorder=3)

    # Get the color of the airfoil plot to use for the dotted line
    airfoil_color = airfoil_line.get_color()    

    # Plot the dotted line between the closest points to y/c = 0.25
    ax.plot([closest_y_c_upper_1, closest_y_c_upper_1], [z_c_upper_dotted_1, z_c_lower_dotted_1], '--', color=airfoil_color)

    # Plot the dotted line between the closest points to y/c = 0.8
    ax.plot([closest_y_c_upper_2, closest_y_c_upper_2], [z_c_upper_dotted_2, z_c_lower_dotted_2], '--', color=airfoil_color)
    
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
    ax.plot([max_pair[0][0]/chord_length_arrondi, max_pair[1][0]/chord_length_arrondi],
            [max_pair[0][1]/chord_length_arrondi, max_pair[1][1]/chord_length_arrondi],
            linestyle='--', color='orange', linewidth=2)
        
    # Set equal scaling for both axes (same units for y and z axes)
    ax.set_aspect('equal')
    
    # Labels and title
    ax.set_xlabel('$y/c$')
    ax.set_ylabel('$z/c$')
    ax.set_title(f'SC(2)-0710 Airfoil with {n_booms} booms')
    ax.grid()
    #fig.savefig(fname='airfoil_y_c_centroid_location.pdf')
    plt.show()

    
    # ---- Plotting (Real Size) ----
    fig2, ax2 = plt.subplots(figsize=(11, 5), dpi=300)
    ax2.plot(y_all, z_all, label='Airfoil')  # Airfoil in meters
    ax2.scatter(y_booms, z_booms, color='red', zorder=3)

    # Plot the dotted lines in meters
    y_dotted_1 = closest_y_c_upper_1 * chord_length_arrondi
    z_dotted_1_top = z_c_upper_dotted_1 * chord_length_arrondi
    z_dotted_1_bot = z_c_lower_dotted_1 * chord_length_arrondi

    y_dotted_2 = closest_y_c_upper_2 * chord_length_arrondi
    z_dotted_2_top = z_c_upper_dotted_2 * chord_length_arrondi
    z_dotted_2_bot = z_c_lower_dotted_2 * chord_length_arrondi

    ax2.plot([y_dotted_1, y_dotted_1], [z_dotted_1_top, z_dotted_1_bot], '--', color=airfoil_color)
    ax2.plot([y_dotted_2, y_dotted_2], [z_dotted_2_top, z_dotted_2_bot], '--', color=airfoil_color)

    # Plot centroid in real coordinates
    ax2.scatter(y_centroid, z_centroid, color='blue', zorder=4)
    ax2.text(y_centroid - 0.25, z_centroid - 0.125, 'C', color='blue', fontsize=12)

    # Plot coordinate axes at the centroid
    ax2.arrow(y_centroid, z_centroid, 0.05 * chord_length_arrondi, 0, head_width=0.02 * chord_length_arrondi, head_length=0.02 * chord_length_arrondi, fc='blue', ec='blue')
    ax2.text(y_centroid + 0.4, z_centroid - 0.05, 'y', color='blue', fontsize=12)
    ax2.arrow(y_centroid, z_centroid, 0, 0.05 * chord_length_arrondi, head_width=0.02 * chord_length_arrondi, head_length=0.02 * chord_length_arrondi, fc='blue', ec='blue')
    ax2.text(y_centroid - 0.05 , z_centroid + 0.4, 'z', color='blue', fontsize=12)

    ax2.set_aspect('equal')
    ax2.set_xlabel('y [m]')
    ax2.set_ylabel('z [m]')
    ax2.set_title(f'SC(2)-0710 Airfoil in Real Scale ({n_booms} booms)')
    ax2.grid()
    #ax2.legend()
    #fig2.savefig(fname='airfoil_y_centroid_location.pdf')
    plt.show()


    # ---- Third Plot: Minimal Clean Version ----
    fig3, ax3 = plt.subplots(figsize=(11, 5), dpi=300)
    ax3.plot(y_all, z_all, color='#f7944d')  # Airfoil outline in black
    ax3.scatter(y_booms, z_booms, color='black', s=11, zorder=3)  # Booms as small black dots

    # Plot the cell division lines
    ax3.plot([y_dotted_1, y_dotted_1], [z_dotted_1_top, z_dotted_1_bot], '--', color='#f7944d', linewidth=1)
    ax3.plot([y_dotted_2, y_dotted_2], [z_dotted_2_top, z_dotted_2_bot], '--', color='#f7944d', linewidth=1)

    # Plot the centroid as a small point only
    ax3.scatter(y_centroid, z_centroid, color='black', s=8, zorder=2)

    # Remove extras
    ax3.set_aspect('equal')
    ax3.axis('off')  # No axis, ticks, grid, labels, or title

    # Save the clean plot
    fig3.savefig('airfoil_minimal_clean.pdf', bbox_inches='tight', pad_inches=0)
    plt.show()

# -------------------------------------------------------------------

plotAirfoil(True, n_booms=37, y_c_max=0.7, y_c_cell_1=0.25) 
