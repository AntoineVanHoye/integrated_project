
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


#geometry wings
beta = 1 #setting angle (degree ou radians ?)
(dx_ac,dy_ac,dz_ac) = (1,1,1) #distance of the aerodynamic center to the root wing
(dx_w_cg,dy_w_cg,dz_w_cg) = (2,2,2) #y placement of the center of gravity of the wing to the root of the regular wing
#(dxa_emp,dya_emp,dza_emp) = (0,0,0) #distance of the aerodynamic center to the root wing of the empenage

#Forces
Lwt = 1  #lift of the wings tappered ( talk to the professor wether we need to separate the differents kind of lifts)
Lw = 1   #lift of the wing overall
Dwt = 1   #drag of the wings tappered ( talk to the professor wether we need to separate the differents kind of drags)
Dw = 1   #drag of the wing overall
Lfus = 1 #lift of the fuselage
Dfus = 1 #drag of the fuselage

P =  1   #lift empenage (add drag ?)
Ffin = 1 #lift fin (formulas of the course ?)
T = 1   #thrust (data from Amos)

#Weight
Wwt = 1  #weight of the wing tappered
Ww = 1    #weight of the wing overall
Wfus = 1  #weight of the fuselage
#Moments
Mwt = 1    #pitch down moment wings tappered
Mw = 1     #pitch down moment wing overall
Mfus = 1  #pitch down moment fuselage
Memp = 1  #pitch down moment empenage

# -------------------------------------------------------------------

def structural_loads_regular_wing (n, alpha): # There will be more parameters as the lift also varies, etc
    
    T_x = (n*Wwt/2-Lwt/2)*np.sin(alpha+beta) + Dwt/2*np.cos(alpha+beta) #alpha is the angle of attack, beta is the setting angle
    T_y = 0
    T_z = (-n*Wwt/2+Lwt/2)*np.cos(alpha+beta) + Dwt/2*np.sin(alpha+beta)
    M_x = 1/2*(-n*(Wwt*dy_w_cg)+Lwt*dy_ac)*np.cos(alpha+beta) - Dwt/2*dy_ac*np.sin(alpha+beta)#dyw is the y placement of the center of gravity to the root of the regular wing , dya is the distance of the aerodynamic center to the root wing
    M_y = 1/2*(-n*Wwt*dx_w_cg + Lwt*dx_ac -Dwt*dz_ac)*np.cos(alpha+beta) + 1/2*(-n*Wwt*dx_w_cg +Lwt*dx_ac + Dwt*dz_ac)*np.sin(alpha+beta) + Mw/2 
    M_z = 1/2*(+n*(Wwt*dy_w_cg)-Lwt*dy_ac)*np.sin(alpha+beta) + Dwt/2*dy_ac*np.cos(alpha+beta)
    return (T_x,T_y,T_z,M_x,M_y,M_z)

# -------------------------------------------------------------------

def dist_2_booms(x_coord, y_coord): 
    dist = [] 
    max_dist = 0
    max_pair = (None, None)  # Will hold the indices (or points) of the max distance

    for i in range(1, len(x_coord)):
        dx = x_coord[i] - x_coord[i-1]
        dy = y_coord[i] - y_coord[i-1]
        d = np.sqrt(dx**2 + dy**2)
        dist.append(d)

        if d > max_dist:
            max_dist = d
            max_pair = ((x_coord[i-1], y_coord[i-1]), (x_coord[i], y_coord[i]))

    return dist, max_dist, max_pair 

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

def skin_thickness(B, sigma_yy, delta_x, delta_z, T_x, T_z, x_booms_ordered_centroid, z_booms_ordered_centroid, x_booms_cell_1, z_booms_cell_1) :
    
    # Taper effect (suite) :
    T_x_web = T_x - B* np.sum(sigma_yy * delta_x)
    T_z_web = T_z - B* np.sum(sigma_yy * delta_z)
    
    # Inertia per unit area (the boom area 'B' is considered to be the same everywhere)
    I_xx_over_B = np.sum(np.array(z_booms_ordered_centroid)**2) 
    I_zz_over_B = np.sum(np.array(x_booms_ordered_centroid)**2) 
    I_xz_over_B = np.sum(np.array(x_booms_ordered_centroid) * np.array(z_booms_ordered_centroid)) 
    
    # ---- Open shear flow ---- 
    
    # -- Cell 1 --
    q_0_cell_1 = [0] # Cut in cell 1, the open shear flow between booms 1 and 2 is zero
    
    #for x in x_booms_cell_1 : 
    
    # -- Cell 2 -- 
    
    
    # ---- Correction term ----
    
    
    # ---- Shear flow (closed) ----
    
    
    # ---- Thickness computation ----
     
    thickness = 1
    
    return thickness
# -------------------------------------------------------------------

def plotAirfoil(plot_airfoil, n_booms, x_c_max, x_c_cell_1):
    if not plot_airfoil:
        return
    
    # ---- Geometry ----
    wing_semi_span = 10 #[m]
    x_tip = wing_semi_span
    x_root = 0
    x_arrondi = 3.10921985816
    chord_root = 5.875
    chord_tip = 2.35
    chord_length_arrondi = chord_tip + (chord_root - chord_tip)/(x_root - x_tip) * (x_arrondi -x_tip)
    # print('chord_length :',chord_length_arrondi)
    sweep_angle = np.radians(31.599) #[°]
    x_dist_root_tip = (wing_semi_span - x_arrondi)*np.tan(sweep_angle)
    print('x_dist_root_tip : ', x_dist_root_tip)
    print('')
    
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
    print('Centroid location:')
    print(f'x = {x_centroid}')
    print(f'z = {z_centroid}')
    
    
    # ---- Verify that the distance between 2 succesive booms is between 0.1 and 0.2 m (Limits imposed by M.Noels) ----
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
    print('x_booms_ordered_centroid =', x_booms_ordered_centroid)
    print('z_booms_ordered_centroid =', z_booms_ordered_centroid)
    
    
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

    
    # ---- Structural loads ---- (will be imported from another fct, here just expl values)
    
    # Here, import all the n and all the alpha from the envelope : 
    # Import from Louis' code
    n = [1,1,1,1,1]
    alpha = [1,1,1,1,1] # degree ou rad ??
    
    
    # ---- Material ----
    
    sigma_y_0 = 1 # CHOOSE THE MATERIAL
    safety_factor = 1.5
    
    
    # ---- Boom area ---- (different pts of the envelope) 
    
    # Point A
    T_x_A = structural_loads_regular_wing(n[0], alpha[0])[0]
    T_y_A = structural_loads_regular_wing(n[0], alpha[0])[1]
    T_z_A = structural_loads_regular_wing(n[0], alpha[0])[2]
    M_x_A = structural_loads_regular_wing(n[0], alpha[0])[3]
    M_y_A = structural_loads_regular_wing(n[0], alpha[0])[4]
    M_z_A = structural_loads_regular_wing(n[0], alpha[0])[5]
    B_pt_A = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_A, M_z_A,sigma_y_0, safety_factor)[0]
    sigma_yy_A = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_A, M_z_A,sigma_y_0, safety_factor)[1]
    
    # Point B
    T_x_B = structural_loads_regular_wing(n[1], alpha[1])[0]
    T_y_B = structural_loads_regular_wing(n[1], alpha[1])[1]
    T_z_B = structural_loads_regular_wing(n[1], alpha[1])[2]
    M_x_B = structural_loads_regular_wing(n[1], alpha[1])[3]
    M_y_B = structural_loads_regular_wing(n[1], alpha[1])[4]
    M_z_B = structural_loads_regular_wing(n[1], alpha[1])[5]
    B_pt_B = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_B, M_z_B,sigma_y_0, safety_factor)[0]
    sigma_yy_B = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_B, M_z_B,sigma_y_0, safety_factor)[1]
    
    # Point C
    T_x_C = structural_loads_regular_wing(n[2], alpha[2])[0]
    T_y_C = structural_loads_regular_wing(n[2], alpha[2])[1]
    T_z_C = structural_loads_regular_wing(n[2], alpha[2])[2]
    M_x_C = structural_loads_regular_wing(n[2], alpha[2])[3]
    M_y_C = structural_loads_regular_wing(n[2], alpha[2])[4]
    M_z_C = structural_loads_regular_wing(n[2], alpha[2])[5]
    B_pt_C = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_C, M_z_C,sigma_y_0, safety_factor)[0]
    sigma_yy_C = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_C, M_z_C,sigma_y_0, safety_factor)[1]
    
    # Point D
    T_x_D = structural_loads_regular_wing(n[3], alpha[3])[0]
    T_y_D = structural_loads_regular_wing(n[3], alpha[3])[1]
    T_z_D = structural_loads_regular_wing(n[3], alpha[3])[2]
    M_x_D = structural_loads_regular_wing(n[3], alpha[3])[3]
    M_y_D = structural_loads_regular_wing(n[3], alpha[3])[4]
    M_z_D = structural_loads_regular_wing(n[3], alpha[3])[5]
    B_pt_D = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_D, M_z_D,sigma_y_0, safety_factor)[0]
    sigma_yy_D = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_D, M_z_D,sigma_y_0, safety_factor)[1]
    
    # Point E
    T_x_E = structural_loads_regular_wing(n[4], alpha[4])[0]
    T_y_E = structural_loads_regular_wing(n[4], alpha[4])[1]
    T_z_E = structural_loads_regular_wing(n[4], alpha[4])[2]
    M_x_E = structural_loads_regular_wing(n[4], alpha[4])[3]
    M_y_E = structural_loads_regular_wing(n[4], alpha[4])[4]
    M_z_E = structural_loads_regular_wing(n[4], alpha[4])[5]
    B_pt_E = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_E, M_z_E,sigma_y_0, safety_factor)[0]
    sigma_yy_E = boom_area(z_booms_ordered_centroid, x_booms_ordered_centroid, M_x_E, M_z_E,sigma_y_0, safety_factor)[1]

    # Max value (most critical)
    B_values = [B_pt_A, B_pt_B, B_pt_C, B_pt_D, B_pt_E]
    sigma_yy_values = [sigma_yy_A, sigma_yy_B, sigma_yy_C, sigma_yy_D, sigma_yy_E]

    B = max(B_values)
    index_of_max = B_values.index(B)
    sigma_yy = sigma_yy_values[index_of_max]
      
    print(f"Index de la contrainte max : {index_of_max}")
    print(f"sigma_yy associée (à B max) : {sigma_yy}")

    print('The boom area is :', B)
    # print('The associated sigma_yy are :', sigma_yy) # should be an array for each boom
    
    
    # ---- Skin thickness ---- 
    
    # Taper effect along x
    x_pos_booms_tip = np.array(x_c_booms_ordered) * chord_tip + x_dist_root_tip 
    x_pos_booms_arrondi = np.array(x_c_booms_ordered) * chord_length_arrondi # where chord_length_arrondi is the chord from the tip, just before the 'arrondi'
    delta_x = x_pos_booms_tip - x_pos_booms_arrondi 
    print('x_pos_booms_tip :', x_pos_booms_tip)
    print('')
    print('x_pos_booms_arrondi :', x_pos_booms_arrondi)
    print('')
    print('delta_x :', delta_x) 
    print('')
    
    # Taper effect along z
    z_pos_booms_tip = np.array(z_c_booms_ordered) * chord_tip 
    z_pos_booms_arrondi = np.array(z_c_booms_ordered) * chord_length_arrondi # where chord_length_arrondi is the chord from the tip, just before the 'arrondi'
    delta_z = z_pos_booms_tip - z_pos_booms_arrondi 
    print('z_pos_booms_tip :', z_pos_booms_tip)
    print('')
    print('z_pos_booms_arrondi :', z_pos_booms_arrondi)
    print('')
    print('delta_z :', delta_z)
    print('')
    
    # Point A
    thickness_pt_A = skin_thickness(B, sigma_yy, delta_x, delta_z, T_x_A, T_z_A, x_booms_ordered_centroid, z_booms_ordered_centroid, x_booms_cell_1, z_booms_cell_1)
    
    # Point B
    thickness_pt_B = skin_thickness(B, sigma_yy, delta_x, delta_z, T_x_B, T_z_B, x_booms_ordered_centroid, z_booms_ordered_centroid, x_booms_cell_1, z_booms_cell_1)
    
    # Point C
    thickness_pt_C = skin_thickness(B, sigma_yy, delta_x, delta_z, T_x_C, T_z_C, x_booms_ordered_centroid, z_booms_ordered_centroid, x_booms_cell_1, z_booms_cell_1)
    
    # Point D
    thickness_pt_D = skin_thickness(B, sigma_yy, delta_x, delta_z, T_x_D, T_z_D, x_booms_ordered_centroid, z_booms_ordered_centroid, x_booms_cell_1, z_booms_cell_1)
    
    # Point E
    thickness_pt_E = skin_thickness(B, sigma_yy, delta_x, delta_z, T_x_E, T_z_E, x_booms_ordered_centroid, z_booms_ordered_centroid, x_booms_cell_1, z_booms_cell_1)
    
    # Take the maximum value
    thickness_values = [thickness_pt_A, thickness_pt_B, thickness_pt_C, thickness_pt_D, thickness_pt_E]
    thickness = max(thickness_values)
    print('')
    print('The skin thickness is ', thickness)
    print('')
    
    
    # ---- Plotting ----
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

# Example usage: 27 booms from x/c = 0 to 0.8
plotAirfoil(True, n_booms=37, x_c_max=0.7, x_c_cell_1=0.25)  # For now, it only works with odd numbers of booms
