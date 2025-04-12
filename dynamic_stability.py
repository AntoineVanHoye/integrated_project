import numpy as np 

##################################################################
######CONFIGURATION SETTING
##################################################################

config = 3
fuel = 2

##################################################################
######GENERAL PARAMETERS     
##################################################################

M = 0.9
beta = np.sqrt(1-M**2)
e = 0.85 # Oswald efficiency factor
alpha_e = 0*np.pi/180
delta = 0.005
rho = 0.28724050871107903
V0 = 265.5380870986307
g = 9.81
CL = 0.45
CD = 0

##################################################################
######WING PARAMETERS      
##################################################################

sweep_quarter_wing = 29.0*np.pi/180
surface_wing_ideal = 119.05
sweep_LE_wing = 31.599 *np.pi/180
surf_wing = 71.11
cl_alpha_wing = 6.680687891225398
AR_wing = 7.04
taper_ratio_wing = 0.492
sweep_half_chord_wing = np.arctan(np.tan(sweep_quarter_wing)+4/AR_wing*(1-taper_ratio_wing)/(1+taper_ratio_wing)*(0.25-0.5))
MAC_wing = 3.706

##################################################################
######FUSELAGE PARAMETERS      
##################################################################

surf_fus = 127.07
sweep_LE_fus = 50.0*np.pi/180
cl_alpha_fus = 5.641342450858086
b_fus = 9
AR_fus = 0.637
taper_ratio_fus = 0.681
sweep_quarter_fus = 41.791 *np.pi/180
sweep_half_chord_fus = np.arctan(np.tan(sweep_quarter_fus)+4/AR_fus*(1-taper_ratio_fus)/(1+taper_ratio_fus)*(0.25-0.5))

##################################################################
######GENERAL GEOMETRICAL PARAMETERS      
##################################################################

b = 28.95*3.28084 # Wingspan in ft
b_wing = b - b_fus
l_cabin = 10.1
cabin_width = 9       
l_cockpit = 2.01
AR = 7.04
surf_tot = 198.17
MAC_tot = 11.098
x_AC_tot = 8.556
cl_alpha = (cl_alpha_wing * b_wing + cl_alpha_fus * b_fus) / b
kappa = cl_alpha/(2*np.pi)
sweep_half_chord = (sweep_half_chord_fus*surf_fus + sweep_half_chord_wing*surf_wing)/surf_tot
CL_alpha = 0

##################################################################
######TAIL PARAMETERS
##################################################################

z_f = 1.51 #distance between AC of the tail and AC of the aircraft
surf_vert_tail = 19
c_root_tail = 3.5
MAC_tail = 2.5120665203177994
x_AC_tail_local = 2.660991632261686
x_AC_tail = l_cabin + l_cockpit + x_AC_tail_local + 1
eta_tail = 0
hor_tail_surf = 35
AR_tail = 5.079670344831602
CL_alpha_tail = 2.745145450926034

##################################################################
######INERTIAS
##################################################################

I_y = 165669.0  / 2.20462
I_x = 0
I_xz = 0
I_z = 0

##################################################################
######PARAMETERS TO CHANGE
##################################################################

Kn = 0
Le = 0
m = 0
x_CG_tot = 8.5

l_f = x_AC_tail - x_CG_tot
V_tail = hor_tail_surf * (x_AC_tail - x_CG_tot)/(surf_tot* MAC_tot)

##################################################################
######COMPUTATION OF METHOD PARAMETERS
##################################################################

U_e = V0 * np.cos(alpha_e)
W_e = V0 * np.sin(alpha_e)
theta_e = 0

##################################################################
######LONGITUDINAL DYNAMIC DERIVATIVES
##################################################################

##################################################################
######ALPHA DERIVATIVES 
##################################################################

def alpha_der(i,d):

    CL_alpha = 2*np.pi*AR/(2+np.sqrt(4+(AR**2 *beta**2 * kappa**2)*(1 + np.tan(sweep_half_chord)**2/beta**2)))
    CD_alpha = 2*CL*CL_alpha/(np.pi*AR*e)
    CM_alpha = -Kn * CL_alpha

    return CL_alpha, CD_alpha, CM_alpha


##################################################################
######U DERIVATIVES 
##################################################################

def u_der1(i,d): 

    CL_alpha = alpha_der(i,d)[0]
    CM_alpha = alpha_der(i,d)[2]
    CL_M = 0
    CD_alpha = alpha_der(i,d)[1]
    CD_M = 0
    CM_alpha = alpha_der(i,d)[2]
    x_AC_M = 0

    CL_u = 2*CL - alpha_e*CL_alpha + M*CL_M
    CD_u = 2*CD - alpha_e*CD_alpha + M*CD_M
    CM_u = -alpha_e * CM_alpha + CL * x_AC_M

    return CL_u, CD_u, CM_u

def u_der2(i,d): 

    CL_u = u_der1(i,d)[0]
    CD_u = u_der1(i,d)[1]
    CM_u = u_der1(i,d)[2]

    Z_u = -CL_u *np.cos(alpha_e) - CD_u *np.sin(alpha_e)
    X_u = CL_u *np.sin(alpha_e) - CD_u* np.cos(alpha_e)
    M_u = CM_u

    X_u = X_u * (1/2*rho*V0*surface_wing_ideal)
    Z_u = Z_u * (1/2*rho*V0*surface_wing_ideal)
    M_u = M_u * (1/2*rho*V0*surface_wing_ideal*MAC_tot)

    return Z_u, X_u, M_u

##################################################################
######Q DERIVATIVES 
##################################################################

def q_der1(i,d):

    CL_alpha = alpha_der(i,d)[0]   
    l_T = x_AC_tail - x_CG_tot
    X_W = np.abs(x_AC_tot - x_CG_tot)
    b = 28.95*3.28084
    mean_chord = AR_tail/b

    CL_q_wing_M0 = (1/2 + 2*X_W/mean_chord)*CL_alpha
    B = np.sqrt(1 - M**2 * np.cos(sweep_quarter_wing)**2)
    CL_q_wing = (AR_tail + 2*np.cos(sweep_quarter_wing))/(AR_tail*B + 2*np.cos(sweep_quarter_wing))*CL_q_wing_M0
    CL_q_tail = 2*CL_alpha_tail*eta_tail * V_tail

    K = 0.7
    CM_q_wing_M0 = -K * cl_alpha_wing * np.cos(sweep_quarter_wing) * ((AR_tail*(2 * (X_W/MAC_wing)**2 + 1/2 *(X_W/MAC_wing)))(AR + 2*np.cos(sweep_quarter_wing)) + 1/24 * (AR_tail**3 * np.tan(sweep_quarter_wing)**2)/(AR_tail + 6*np.cos(sweep_quarter_wing) + 1/8))
    CM_q_wing = (((AR_tail**3 * np.tan(sweep_quarter_wing)**2)/(AR_tail*B + 6*np.cos(sweep_quarter_wing)) + 3/B)/(((AR_tail**3 * np.tan(sweep_quarter_wing)**2)/(AR_tail*B + 6*np.cos(sweep_quarter_wing))) + 3)) * CM_q_wing_M0
    CM_q_tail = -2*CL_alpha_tail*eta_tail*V_tail*l_T/MAC_tail

    CD_q = 0 #approximation
    CL_q = CL_q_wing + CL_q_tail #only horizontal tail
    CM_q = CM_q_wing + CM_q_tail

    return CL_q, CD_q, CM_q

def q_der2(i,d):

    CL_q = q_der1(i,d)[0]
    CD_q = q_der1(i,d)[1]
    CM_q = q_der1(i,d)[2]

    Z_q = -1/2*(CL_q*np.cos(alpha_e) + CD_q*np.sin(alpha_e))
    X_q = 1/2*(CL_q*np.sin(alpha_e) - CD_q*np.cos(alpha_e))
    M_q = 1/2*CM_q

    Z_q = Z_q * (1/2*rho*V0*surface_wing_ideal*MAC_tot)
    X_q = X_q * (1/2*rho*V0*surface_wing_ideal*MAC_tot)
    M_q = M_q * (1/2*rho*V0*surface_wing_ideal*MAC_tot**2)

    return Z_q, X_q, M_q

##################################################################
######ALPHA DOT DERIVATIVES 
##################################################################

def alpha_dot_der(i,d):

    deps_dalpha = 0 #approximation
    l_T = x_AC_tail - x_CG_tot
    CD_alpha_dot = 0

    CL_alpha_dot_W = 0
    CL_alpha_dot_H = 2*CL_alpha_tail*eta_tail*V_tail*deps_dalpha
    CL_alpha_dot = CL_alpha_dot_W + CL_alpha_dot_H

    CM_alpha_dot_W = 0
    CM_alpha_dot_H = -2*CL_alpha_tail*eta_tail*V_tail*deps_dalpha*l_T/MAC_tail
    CM_alpha_dot = CM_alpha_dot_W + CM_alpha_dot_H

    return CL_alpha_dot, CD_alpha_dot, CM_alpha_dot

##################################################################
######W DOT DERIVATIVES 
##################################################################

def w_dot_der(i,d): 

    CL_alpha_dot = alpha_dot_der(i,d)[0]
    CD_alpha_dot = alpha_dot_der(i,d)[1]
    CM_alpha_dot = alpha_dot_der(i,d)[2]

    X_w_dot = 1/(2*np.cos(alpha_e)) * (CL_alpha_dot*np.sin(alpha_e) - CD_alpha_dot*np.cos(alpha_e))
    Z_w_dot = 1/(2*np.cos(alpha_e)) * (-CL_alpha_dot*np.cos(alpha_e) - CD_alpha_dot*np.sin(alpha_e))
    M_w_dot = 1/(2*np.cos(alpha_e)) * CM_alpha_dot

    X_w_dot = X_w_dot * (1/2*rho*surf_tot*MAC_tot)
    Z_w_dot = Z_w_dot * (1/2*rho*surf_tot*MAC_tot)
    M_w_dot = M_w_dot * (1/2*rho*surf_tot*MAC_tot**2)

    return X_w_dot, Z_w_dot, M_w_dot

##################################################################
######W DERIVATIVES 
##################################################################

def w_der(i,d):

    CL_alpha = alpha_der(i,d)[0]
    CD_alpha = alpha_der(i,d)[1]
    CM_alpha = alpha_der(i,d)[2]
   
    De = CD*1/2*rho*V0**2*surface_wing_ideal
    Xe = Le*np.sin(alpha_e) - De*np.cos(alpha_e)
    Ze = -Le*np.cos(alpha_e) - De*np.sin(alpha_e)

    C_Xe = Xe/(1/2*rho*V0**2*surf_tot)
    C_Ze = Ze/(1/2*rho*V0**2*surf_tot)

    Z_w = 1/np.cos(alpha_e) * (C_Xe -CL_alpha*np.cos(alpha_e) - CD_alpha*np.sin(alpha_e))
    X_w = 1/np.cos(alpha_e) * (-C_Ze + CL_alpha*np.sin(alpha_e) - CD_alpha*np.cos(alpha_e))
    M_w = 1/np.cos(alpha_e) * CM_alpha

    Z_w = Z_w * (1/2*rho*V0*surface_wing_ideal*MAC_tot)
    X_w = X_w * (1/2*rho*V0*surface_wing_ideal*MAC_tot)
    M_w = M_w * (1/2*rho*V0*surface_wing_ideal*MAC_tot**2)

    return Z_w, X_w, M_w

##################################################################
######CONSTRUCTION OF THE A MATRIX FOR THE LONGITUDINAL MOTION
##################################################################

def long_dyn_stab(i,d,Cl): 

    X_w_dot = w_dot_der(i,d)[0]
    Z_w_dot = w_dot_der(i,d)[1]
    M_w_dot = w_dot_der(i,d)[2]

    X_u = u_der2(i,d)[1]
    Z_u = u_der2(i,d)[0]
    M_u = u_der2(i,d)[2]

    X_w = w_der(i,d)[1]
    Z_w = w_der(i,d)[0]
    M_w = w_der(i,d)[2]

    X_q = q_der2(i,d)[1]
    Z_q = q_der2(i,d)[0]
    M_q = q_der2(i,d)[2]
    
    M_matrix = np.array([
        [m, -X_w_dot, 0, 0],
        [0, (m - Z_w_dot), 0, 0],
        [0, -M_w_dot, I_y, 0],
        [0, 0, 0, 1]
    ])

    B_matrix = np.array([
        [-X_u, -X_w, -(X_q - m * W_e), m * g * np.cos(theta_e)],
        [-Z_u, -Z_w, -(Z_q + m * U_e), m * g * np.sin(theta_e)],
        [-M_u, -M_w, -M_q, 0],
        [0, 0, -1, 0]
    ])

    A_matrix = -np.linalg.inv(M_matrix) @ B_matrix
    eigenvalues = np.linalg.eigvals(A_matrix)
    real_parts = np.real(eigenvalues)

    conjugates = np.conj(eigenvalues)  

    all_conjugates_exist = np.all(np.isin(conjugates, eigenvalues))

    if all_conjugates_exist:
        test = 1
    else:
        test = 0
    
    return real_parts, test

##################################################################
######LATERAL DYNAMIC DERIVATIVES
##################################################################

def dsigma_dbeta(sweep_quarter_wing,AR): 

    average_fus_sect = 0
    d = np.sqrt(average_fus_sect/0.7854)
    Z_w = 0#vertical distance from the wing root quarter-chord to the fuselage centerline  
    d_sigma = -0.276 + 3.06*surf_vert_tail/surf_tot /(1 + np.cos(sweep_quarter_wing)) + 0.4*Z_w + 0.009*AR
    return d_sigma

##################################################################
######V DERIVATIVES 
##################################################################

def v_der(sweep_quarter_wing,AR):

    d_sigma = dsigma_dbeta(sweep_quarter_wing,AR)
    c_1 = 0#fin lift curve slope
    Y_v = 1/2*rho*V0**2 *surf_vert_tail * c_1 * (1-d_sigma)
    L_v = -1/2 * rho*V0**2 * surf_vert_tail * z_f * c_1 * (1-d_sigma)
    N_v = 1/2 * rho*V0*surf_vert_tail * l_f * c_1 * (1-d_sigma)

    Y_v = Y_v * (1/2*rho*V0*surface_wing_ideal)
    L_v = L_v * (1/2*rho*V0*surface_wing_ideal*b)
    N_v = N_v * (1/2*rho*V0*surface_wing_ideal*b)

    return Y_v, L_v, N_v

##################################################################
######R DERIVATIVES 
##################################################################

def r_der(sweep_quarter_wing,AR):

    Y_v = v_der(sweep_quarter_wing,AR)[0]
    Y_r = -2*Y_v * (z_f * np.sin(alpha_e) - l_f * np.cos(alpha_e)/b)
    L_r = -2*Y_v*((z_f*np.sin(alpha_e) + l_f*np.cos(alpha_e)*(z_f*np.cos(alpha_e) - l_f*np.sin(alpha_e)))/b**2)
    N_r = 2*Y_v * ((z_f*np.sin(alpha_e + l_f*np.cos(alpha_e))**2)/b**2)

    Y_r = Y_r * (1/2*rho*V0*surface_wing_ideal*b)
    L_r = L_r * (1/2*rho*V0*surface_wing_ideal*b**2)
    N_r = N_r * (1/2*rho*V0*surface_wing_ideal*b**2)

    return Y_r, L_r, N_r

##################################################################
######P DERIVATIVES 
##################################################################

def p_der(sweep_quarter_wing,AR):
    
    Y_v = v_der(sweep_quarter_wing,AR)[0]
    L_r = r_der(sweep_quarter_wing,AR)[1]
    Y_p = 2*Y_v * (z_f * np.cos(alpha_e) - l_f * np.sin(alpha_e)/b)
    L_p = 2 * Y_v * (z_f/b)**2
    N_p = L_r

    Y_p = Y_p * (1/2*rho*V0*surface_wing_ideal*b)   
    L_p = L_p * (1/2*rho*V0*surface_wing_ideal*b**2)
    N_p = N_p * (1/2*rho*V0*surface_wing_ideal*b**2)

    return Y_p, L_p, N_p


##################################################################
######CONSTRUCTION OF THE A MATRIX FOR THE LATERAL MOTION
##################################################################


def lat_dyn_stab(i,d,Cl): 

    Y_v = v_der(sweep_quarter_wing,AR)[0]
    Y_p = p_der(sweep_quarter_wing,AR)[0]
    Y_r = r_der(sweep_quarter_wing,AR)[0]

    L_v = v_der(sweep_quarter_wing,AR)[1]
    L_p = p_der(sweep_quarter_wing,AR)[1]
    L_r = r_der(sweep_quarter_wing,AR)[1]

    N_v = v_der(sweep_quarter_wing,AR)[2]
    N_p = p_der(sweep_quarter_wing,AR)[2]
    N_r = r_der(sweep_quarter_wing,AR)[2]

    M_matrix = np.array([
    [m, 0, 0, 0, 0],
    [0, I_x, -I_xz, 0, 0],
    [0, -I_xz, I_z, 0, 0],
    [0, 0, 0, 1, 0],
    [0, 0, 0, 0, 1]
    ])


    B_matrix = np.array([
        [-Y_v, -(Y_p + m * W_e), -(Y_r - m * U_e), -m * g * np.cos(theta_e), -m * g * np.sin(theta_e)],
        [-L_v, -L_p, -L_r, 0, 0],
        [-N_v, -N_p, -N_r, 0, 0],
        [0, -1, 0, 0, 0],
        [0, 0, -1, 0, 0]
    ])

    # Calcul de la matrice A
    A_matrix = -np.linalg.inv(M_matrix) @ B_matrix
    eigenvalues = np.linalg.eigvals(A_matrix)
    real_parts = np.real(eigenvalues)

    conjugates = np.conj(eigenvalues)  

    all_conjugates_exist = np.all(np.isin(conjugates, eigenvalues))

    if all_conjugates_exist:
        test = 1
    else:
        test = 0
    
    return real_parts, test

##################################################################
######ANALYSIS OF THE PHUGOID MODE
##################################################################

def phugoid(): 

    omega_p = g*np.sqrt(2)/V0
    xi_p = g*CD/(CL*V0)/omega_p

    return omega_p, xi_p

##################################################################
######SHORT PERIOD OSCIELLATIONS MODE
##################################################################

def short_period(i,d,Cl):

    M_q = q_der2(config,fuel)[2]
    M_w = w_der(config,fuel)[2]
    Z_w = w_der(config,fuel)[0]
    M_w_dot = w_dot_der(config,fuel)[2]

    omega_s = np.sqrt((M_q/I_y)*(Z_w/m)-(M_w/I_y)*Ue)
    xi_s = -(M_q/I_y + Z_w/m + M_w_dot/I_y*Ue)/(2*omega_s)

    return omega_s, xi_s
##################################################################
######DEFINITION OF THE PRINT FUNCTION
##################################################################

def main(): 

    real_parts_long, test_long = long_dyn_stab(config,fuel,0.45)
    real_parts_lat, test_lat = lat_dyn_stab(config,fuel,0.45)
    omega_p, xi_p = phugoid()
    omega_s, xi_s = short_period(config,fuel,0.45)

    print("--------------------------CHECK OF THE SYSTEM--------------------------------------------")
    if test_long == 1:
        print("The complex conjugate of all the eigenvalues are also eigenvalues for the longitudinal motion.")
    else:
        print("There exists at least one complex conjugate of an eigenvalue which is not an eigenvalue for the longitudinal motion.")
    
    if test_lat == 1:
        print("The complex conjugate of all the eigenvalues are also eigenvalues for the lateral motion.")
    else: 
        print("There exists at least one complex conjugate of an eigenvalue which is not an eigenvalue for the lateral motion.")
    
    if test_lat == 0 or test_long == 0:
        print("Thus, there is a problem with the definition of the system. The system is not well defined.")
    else: 
        print("Thus, the system is well defined.")

    print("----------------------------------------------------------------------")

    print("--------------------------STABILITY--------------------------------------------")
    if np.all(real_parts_long < 0):
        print("The aircraft is dynamically stable for the longitudinal motion.")
    elif np.any(real_parts_long > 0):
        print("The system is dynamiccaly unstable for the longitudinal motion.")
    elif np.any(real_parts_long == 0):
        print("The system is dynamically neutrally stable for the longitudinal motion.")
    
    if np.all(real_parts_long < 0):
        print("The aircraft is dynamically stable for the lateral motion.")
    elif np.any(real_parts_long > 0):
        print("The system is dynamiccaly unstable for the lateral motion.")
    elif np.any(real_parts_long == 0):
        print("The system is dynamically neutrally stable for the lateral motion.")
    print("----------------------------------------------------------------------")

    print("--------------------------PHUGOID MODE--------------------------------------------")
    print("The frequency of the phugoid mode is: ", omega_p, "rad/s")
    print("The damping ratio of the phugoid mode is: ", xi_p)
    print("----------------------------------------------------------------------")

    print("--------------------------SHORT PERIOD OSCILLATIONS MODE--------------------------------------------")
    print("The frequency of the short period oscillations mode is: ", omega_s, "rad/s")
    print("The damping ratio of the short period oscillations mode is: ", xi_s)
    print("----------------------------------------------------------------------")

if __name__ == "__main__":
    main()
