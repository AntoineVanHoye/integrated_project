import numpy as np

##################################################################
######GENERAL PARAMETERS     
##################################################################

M = 0.9
beta = np.sqrt(1-M**2)
e = 0.65 # Oswald efficiency factor
alpha_e = 0*np.pi/180
delta = 0.005
rho = 0.28724050871107903*(0.0685218/35.3147) #air density in slugs/ft^3
rho_adim = 0.28724050871107903 *(0.0685218/35.3147)
V0 = 265.5380870986307*3.28084 #velocity in ft/s
V0_adim = 265.5380870986307*3.28084
g = 9.81 *3.28084 #gravity in ft/s^2
CL = 0.45
CD = 0.026
CD0 = 0.02
beta_M0 = 1

##################################################################
######WING PARAMETERS      
##################################################################

sweep_quarter_wing = 29.0*np.pi/180
surface_wing_ideal = 119.05*10.7639 #surface of the wing in ft^2
surface_wing_ideal_adim = 119.05*10.7639
sweep_LE_wing = 31.599 *np.pi/180
surf_wing = 71.11*10.7639
cl_alpha_wing = 6.680687891225398
AR_wing = 7.04
taper_ratio_wing = 0.492
sweep_half_chord_wing = np.arctan(np.tan(sweep_quarter_wing)+4/AR_wing*(1-taper_ratio_wing)/(1+taper_ratio_wing)*(0.25-0.5))
MAC_wing = 3.706*3.28084 #mean aerodynamic chord in ft
B = np.sqrt(1 - M**2 * np.cos(sweep_quarter_wing)**2)
wing_aoa = 2.92*np.pi/180
twist_angle = -2.63*np.pi/180

##################################################################
######FUSELAGE PARAMETERS      
##################################################################

surf_fus = 127.07*10.7639
sweep_LE_fus = 50.0*np.pi/180
cl_alpha_fus = 5.641342450858086
b_fus = 9*3.28084
AR_fus = 0.637
taper_ratio_fus = 0.681
sweep_quarter_fus = 41.791 *np.pi/180
sweep_half_chord_fus = np.arctan(np.tan(sweep_quarter_fus)+4/AR_fus*(1-taper_ratio_fus)/(1+taper_ratio_fus)*(0.25-0.5))
av_section_fus = 15.71563982*10.7639 
length_fus = 16.8*3.28084 #length of the fuselage in ft
x1 = 0.4*length_fus
x0 = (0.378 + 0.527*x1/length_fus)*length_fus
section_x0 = 19.77543273327567*10.7639

##################################################################
######GENERAL GEOMETRICAL PARAMETERS      
##################################################################

b = 28.95*3.28084 # Wingspan in ft
b_adim = 28.95*3.28084
b_wing = b - b_fus
l_cabin = 10.1*3.28084
cabin_width = 9*3.28084    
l_cockpit = 2.01*3.28084
AR = 7.04
surf_tot = 198.17*10.7639
MAC_tot = 11.098*3.28084
MAC_tot_adim = 11.098*3.28084
x_AC_tot = 8.556*3.28084
cl_alpha = (cl_alpha_wing * b_wing + cl_alpha_fus * b_fus) / b
kappa = cl_alpha/(2*np.pi)
sweep_half_chord = (sweep_half_chord_fus*b_fus + sweep_half_chord_wing*b_wing)/b
CL_alpha = 2*np.pi*AR/(2+np.sqrt(4+(AR**2 *beta**2 / kappa**2)*(1 + np.tan(sweep_half_chord)**2/beta**2)))
mean_chord = b/AR
gamma = 0 #dihedral angle 

##################################################################
######TAIL PARAMETERS
##################################################################

surf_vert_tail = 19*10.7639
c_root_tail = 3.5*3.28084
MAC_tail = 2.5120665203177994*3.28084
x_AC_tail_local = 2.660991632261686*3.28084
x_AC_tail = l_cabin + l_cockpit + x_AC_tail_local + 1*3.28084 
eta_tail =  0.9
hor_tail_surf = 35*10.7639
span_hor_tail = 12.5*3.28084
AR_tail = 5.079670344831602
AR_hor_tail = 4.464285714285714
AR_vert_tail = 2.4234693877551026
y_AC_tail = 1.8645301924973645*3.28084
sweep_quarter_tail = 34.89785663924308 * np.pi/180
z_w = 0.7*3.28084 #distance from body centerline to quarter chord point of exposed wing root chord, positive for the quarter chord point below the body centerline
z_f = x_AC_tail - x_AC_tot
cl_alpha_tail = (0.00606 - 0.00593)/(2.5-2.25)
kappa_tail = cl_alpha_tail/(2*np.pi)
CL_alpha_tail = 2*np.pi*AR_hor_tail/(2 + np.sqrt(4 + (AR_hor_tail**2*beta**2/kappa_tail**2)*(1 + np.tan(sweep_quarter_tail)**2/beta**2)))
CL_alpha_tail_vert = 2*np.pi*AR_vert_tail/(2 + np.sqrt(4 + (AR_vert_tail**2*beta**2/kappa_tail**2)*(1 + np.tan(sweep_quarter_tail)**2/beta**2)))
CL_tail = 0.20545613522169995

##################################################################
######INERTIAS
##################################################################

I_y = 398469.35 * 0.737562149
I_x = 1070353.72  * 0.737562149
I_xz = 75798.46919250* 0.737562149
I_z = 1344583.40 * 0.737562149

##################################################################
######PARAMETERS TO CHANGE
##################################################################

Kn = 0.1
Le = 584000*0.224809 #lift force in lb
m = 105682.71* 0.031080950037834#*0.45359237 #mass in kg
x_CG_tot = 8.5*3.28084

x = x_AC_tot - x_CG_tot
l_V = x_AC_tail - x_CG_tot
V_tail = hor_tail_surf * (x_AC_tail - x_CG_tot)/(surface_wing_ideal* MAC_tot)
X_W = np.abs(x_AC_tot - x_CG_tot)

##################################################################
######COMPUTATION OF METHOD PARAMETERS
##################################################################

U_e = V0 * np.cos(alpha_e)
W_e = V0 * np.sin(alpha_e)
theta_e = 0

##################################################################
######PARAMETERS TO BE OBTAINED FROM THE GRAPHS 
##################################################################

K = 0.7 #for q derivatives
#for v derivatives
k = 0.9 
Ki = 1.25
term1 = -0.3
term2 = -0.4 #for r derivatives
#for u derivatives 
K1 = 1.3
xAC_cR = 0.95
K2 = 0.8
K1_bis = 1.3
xAC_cR_bis = 0.8
K2_bis = 0.8

##################################################################
######LONGITUDINAL DYNAMIC DERIVATIVES
##################################################################

##################################################################
######ALPHA DERIVATIVES 
##################################################################

def alpha_der():

    CD_alpha = 2*CL*CL_alpha/(np.pi*AR*e)
    CM_alpha = -Kn * CL_alpha

    return CL_alpha, CD_alpha, CM_alpha

##################################################################
######U DERIVATIVES 
##################################################################

def u_der1(): 

    CL_alpha = alpha_der()[0]
    CM_alpha = alpha_der()[2]
    CL_M = 0
    CD_alpha = alpha_der()[1]
    CD_M = 0
    CM_alpha = alpha_der()[2]
    x_AC_M = 0

    CL_u = 2*CL - alpha_e*CL_alpha + M*CL_M
    CD_u = 2*CD - alpha_e*CD_alpha + M*CD_M
    CM_u = -alpha_e * CM_alpha - CL * x_AC_M

    return CL_u, CD_u, CM_u

def u_der2(): 

    CL_u = u_der1()[0]
    CD_u = u_der1()[1]
    CM_u = u_der1()[2]

    Z_u = -CL_u *np.cos(alpha_e) - CD_u *np.sin(alpha_e)
    X_u = CL_u *np.sin(alpha_e) - CD_u* np.cos(alpha_e)
    M_u = CM_u

    X_u = X_u * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim)
    Z_u = Z_u * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim)
    M_u = M_u * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*MAC_tot_adim)

    return Z_u, X_u, M_u

##################################################################
######Q DERIVATIVES 
##################################################################

def q_der1():

    CL_alpha_M0 = 2*np.pi*AR/(2+np.sqrt(4+(AR**2 *beta_M0**2 * kappa**2)*(1 + np.tan(sweep_half_chord)**2/beta_M0**2)))

    CL_q_wing_M0 = (1/2 + 2*X_W/mean_chord)*CL_alpha_M0
    CL_q_wing = (AR + 2*np.cos(sweep_quarter_wing))/(AR*B + 2*np.cos(sweep_quarter_wing))*CL_q_wing_M0
    CL_q_tail = 2*CL_alpha_tail*eta_tail * V_tail

    CM_q_wing_M0 = -K * cl_alpha_wing * np.cos(sweep_quarter_wing) * (((AR*(2 * (X_W/MAC_wing)**2 + 1/2 *(X_W/MAC_wing)))/(AR + 2*np.cos(sweep_quarter_wing))) + (1/24 * (AR**3 * np.tan(sweep_quarter_wing)**2)/(AR + 6*np.cos(sweep_quarter_wing)) + 1/8))
    CM_q_wing = ((((AR**3 * np.tan(sweep_quarter_wing)**2)/(AR*B + 6*np.cos(sweep_quarter_wing))) + 3/B)/(((AR**3 * np.tan(sweep_quarter_wing)**2)/(AR*B + 6*np.cos(sweep_quarter_wing))) + 3)) * CM_q_wing_M0
    CM_q_tail = -2*CL_alpha_tail*eta_tail*V_tail*l_V/MAC_tail

    CD_q = 0 #approximation, /!\ check for transonic 
    CL_q = CL_q_wing + CL_q_tail #only horizontal tail
    CM_q = CM_q_wing + CM_q_tail

    return CL_q, CD_q, CM_q

def q_der2():

    CL_q = q_der1()[0]
    CD_q = q_der1()[1]
    CM_q = q_der1()[2]

    Z_q = -1/2*(CL_q*np.cos(alpha_e) + CD_q*np.sin(alpha_e))
    X_q = 1/2*(CL_q*np.sin(alpha_e) - CD_q*np.cos(alpha_e))
    M_q = 1/2*CM_q

    Z_q = Z_q * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*MAC_tot_adim)
    X_q = X_q * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*MAC_tot_adim)
    M_q = M_q * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*MAC_tot_adim**2)

    return Z_q, X_q, M_q

##################################################################
######ALPHA DOT DERIVATIVES 
##################################################################

def alpha_dot_der():

    deps_dalpha = 0 #approximation
    CD_alpha_dot = 0

    CL_alpha_dot_W = 0
    CL_alpha_dot_H = 2*CL_alpha_tail*eta_tail*V_tail*deps_dalpha
    CL_alpha_dot = CL_alpha_dot_W + CL_alpha_dot_H

    CM_alpha_dot_W = 0
    CM_alpha_dot_H = -2*CL_alpha_tail*eta_tail*V_tail*deps_dalpha*l_V/MAC_tail
    CM_alpha_dot = CM_alpha_dot_W + CM_alpha_dot_H

    return CL_alpha_dot, CD_alpha_dot, CM_alpha_dot

##################################################################
######W DOT DERIVATIVES 
##################################################################

def w_dot_der(): 

    CL_alpha_dot = alpha_dot_der()[0]
    CD_alpha_dot = alpha_dot_der()[1]
    CM_alpha_dot = alpha_dot_der()[2]

    X_w_dot = 1/(2*np.cos(alpha_e)) * (CL_alpha_dot*np.sin(alpha_e) - CD_alpha_dot*np.cos(alpha_e))
    Z_w_dot = 1/(2*np.cos(alpha_e)) * (-CL_alpha_dot*np.cos(alpha_e) - CD_alpha_dot*np.sin(alpha_e))
    M_w_dot = 1/(2*np.cos(alpha_e)) * CM_alpha_dot

    X_w_dot = X_w_dot * (1/2*rho_adim*surface_wing_ideal_adim*MAC_tot_adim)
    Z_w_dot = Z_w_dot * (1/2*rho_adim*surface_wing_ideal_adim*MAC_tot_adim)
    M_w_dot = M_w_dot * (1/2*rho_adim*surface_wing_ideal_adim*MAC_tot_adim**2)

    return X_w_dot, Z_w_dot, M_w_dot

##################################################################
######W DERIVATIVES 
##################################################################

def w_der():

    CL_alpha = alpha_der()[0]
    CD_alpha = alpha_der()[1]
    CM_alpha = alpha_der()[2]
   
    De = CD*1/2*rho*V0**2*surface_wing_ideal
    Xe = Le*np.sin(alpha_e) - De*np.cos(alpha_e)
    Ze = -Le*np.cos(alpha_e) - De*np.sin(alpha_e)

    C_Xe = Xe/(1/2*rho*V0**2*surface_wing_ideal)
    C_Ze = Ze/(1/2*rho*V0**2*surface_wing_ideal)

    Z_w = 1/np.cos(alpha_e) * (C_Xe -CL_alpha*np.cos(alpha_e) - CD_alpha*np.sin(alpha_e))
    X_w = 1/np.cos(alpha_e) * (-C_Ze + CL_alpha*np.sin(alpha_e) - CD_alpha*np.cos(alpha_e))
    M_w = 1/np.cos(alpha_e) * CM_alpha

    Z_w = Z_w * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim)
    X_w = X_w * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim)
    M_w = M_w * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*MAC_tot_adim)

    return Z_w, X_w, M_w

##################################################################
######CONSTRUCTION OF THE A MATRIX FOR THE LONGITUDINAL MOTION
##################################################################

def long_dyn_stab(): 

    X_w_dot = w_dot_der()[0]
    Z_w_dot = w_dot_der()[1]
    M_w_dot = w_dot_der()[2]

    X_u = u_der2()[1]
    Z_u = u_der2()[0]
    M_u = u_der2()[2]

    X_w = w_der()[1]
    Z_w = w_der()[0]
    M_w = w_der()[2]

    X_q = q_der2()[1]
    Z_q = q_der2()[0]
    M_q = q_der2()[2]
    
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

    for i in [0, -1]:  
        eig = eigenvalues[i]
        sigma = eig.real
        omega_d = eig.imag
        omega_n = np.sqrt(sigma**2 + omega_d**2)
        zeta = -sigma / omega_n if omega_n != 0 else 0
        if i == 0 :
            print(f"Short period oscillations mode :")
            print(f"  ω_n = {omega_n:.4f} rad/s")
            print(f"  ζ = {zeta:.4f}")
            if zeta > 0.3 and zeta < 2:
                print("You are level 1\n")
            if zeta > 0.2 and zeta < 0.3:
                print("You are level 2\n")
            if zeta > 0.15 and zeta < 0.2:
                print("You are level 3\n")

        if i == -1 :
            print(f"Phugoid mode :")
            print(f"  ω_n = {omega_n:.4f} rad/s")
            print(f"  ζ = {zeta:.4f}")
            if zeta > 0.04: 
                print("You are level 1\n")
            if zeta >0 and zeta<0.04:
                print("You are level 2\n")

    return real_parts, test, eigenvalues

##################################################################
######LATERAL DYNAMIC DERIVATIVES
##################################################################

##################################################################
######V DERIVATIVES 
##################################################################

def beta_der():
   
    d = np.sqrt(av_section_fus/0.7854)
    term = 0.724 * 3.06*(surf_vert_tail/surface_wing_ideal)/(1+np.cos(sweep_quarter_wing))+0.4*z_w/d + 0.009*AR
    term3 = -0.0025
    KM1 = 1.4
    KF1 = 0.95
    term4 = -0.001
    term5 = -0.0023
    KM2 = 1.23
    KF2 = 0.97
    term6 =-0.001
    KN = 0.0012
    KRl = 1.1
    fus_side_area = 35.08483365482052*10.7639

    Y_beta_W = -0.0001*np.abs(gamma)*57.3
    Y_beta_B = -2*Ki * section_x0/surf_fus
    Y_beta_V = -k*CL_alpha_tail_vert * term * surf_vert_tail/surface_wing_ideal
    
    L_beta_WB = 57.3*(CL*(term3 * KM1 * KF1 + term4))
    L_beta_H = 57.3*(CL_tail*(term5 * KM2 * KF2 + term6))*hor_tail_surf*span_hor_tail/(surface_wing_ideal*b)
    L_beta_V = Y_beta_V * (y_AC_tail*np.cos(alpha_e) - l_V*np.sin(alpha_e))/b

    N_beta_W = 0
    N_beta_B = -57.3*KN*KRl*fus_side_area/surface_wing_ideal *length_fus/b
    N_beta_V = -Y_beta_V * (l_V*np.cos(alpha_e) + y_AC_tail*np.sin(alpha_e))/b
    
    L_beta = L_beta_WB + L_beta_H + L_beta_V
    Y_beta = Y_beta_W + Y_beta_B + Y_beta_V
    N_beta = N_beta_W + N_beta_B + N_beta_V

    Y_beta_normed = Y_beta * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim)
    L_beta = L_beta * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*b_adim)
    N_beta = N_beta * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*b_adim)

    return L_beta, N_beta, Y_beta_normed, Y_beta_V

##################################################################
######R DERIVATIVES 
##################################################################

def r_der():

    Y_beta_V = beta_der()[3]
    Clr_CL_M0 = 0.35
    term_twist = 0.0002
    Lr_CL = ((1 + ((AR*(1-B**2)/(2*B*(AR*B+2*np.cos(sweep_quarter_wing))))) + ((AR*B+2*np.cos(sweep_quarter_wing)/(AR*B+4*np.cos(sweep_quarter_wing))))*np.tan(sweep_quarter_wing)**2/8)/(1 + ((AR + 2*np.cos(sweep_quarter_wing))/(AR + 4*np.cos(sweep_quarter_wing)))*np.tan(sweep_quarter_wing)**2/8))*Clr_CL_M0

    Y_r = -2*Y_beta_V * (y_AC_tail * np.sin(alpha_e) + l_V * np.cos(alpha_e))/b
    
    L_rW = CL *Lr_CL +  term_twist*twist_angle 
    L_rV = -2/b**2*(l_V * np.cos(alpha_e)+y_AC_tail*np.sin(alpha_e))*(y_AC_tail*np.cos(alpha_e) - l_V*np.sin(alpha_e))*Y_beta_V
    
    N_rV = 2/b**2*(l_V*np.cos(alpha_e)+y_AC_tail*np.sin(alpha_e))**2*Y_beta_V
    N_rW = term1 * CL**2 + term2*CD0

    L_r = L_rW + L_rV
    N_r = N_rW + N_rV

    Y_r = Y_r * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*b_adim)
    L_r = L_r * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*b_adim**2)
    N_r = N_r * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*b_adim**2)

    return Y_r, L_r, N_r

##################################################################
######P DERIVATIVES 
##################################################################

def p_der():

    term_twist = 0.0002
    Y_beta_V = beta_der()[3]
    factor = -0.23
    Np_CL_M0 = -1/6 * ((AR + 6*(AR+np.cos(sweep_quarter_wing))*(x/MAC_tot*np.tan(sweep_quarter_wing)/AR + np.tan(sweep_quarter_wing)**2/12))/(AR + 4*np.cos(sweep_quarter_wing)))
    Np_CL = ((AR + 4*np.cos(sweep_quarter_wing))/(AR*B+4*np.cos(sweep_quarter_wing))*((AR*B + 1/2*(AR*B+np.cos(sweep_quarter_wing))*np.tan(sweep_quarter_wing)**2))/(AR+1/2*(AR+np.cos(sweep_quarter_wing))*np.tan(sweep_quarter_wing)**2))*Np_CL_M0

    Y_p = 2*Y_beta_V * (y_AC_tail * np.cos(alpha_e) - l_V * np.sin(alpha_e))/b

    L_p_WB = factor *kappa/beta
    L_p_H = 0.5*factor*kappa/beta*hor_tail_surf/surface_wing_ideal*(span_hor_tail/b)**2
    L_p_V = 2 * Y_beta_V * (y_AC_tail/b)**2

    L_p = L_p_WB + L_p_H + L_p_V

    N_p_V = -2/b * (l_V*np.cos(alpha_e)+y_AC_tail*np.sin(alpha_e))*Y_beta_V*(y_AC_tail*np.cos(alpha_e)-l_V*np.sin(alpha_e))/b
    N_p_W = - L_p_WB*np.tan(wing_aoa) - (-L_p*np.tan(wing_aoa) - Np_CL * CL) + term_twist*twist_angle 
    N_p = N_p_W + N_p_V

    Y_p = Y_p * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*b_adim)   
    L_p = L_p * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*b_adim**2)
    N_p = N_p * (1/2*rho_adim*V0_adim*surface_wing_ideal_adim*b_adim**2)

    return Y_p, L_p, N_p

##################################################################
######CONSTRUCTION OF THE A MATRIX FOR THE LATERAL MOTION
##################################################################

def lat_dyn_stab(): 

    Y_beta = beta_der()[2]
    Y_p = p_der()[0]
    Y_r = r_der()[0]

    L_beta = beta_der()[0]
    L_p = p_der()[1]
    L_r = r_der()[1]

    N_beta = beta_der()[2]
    N_p = p_der()[2]
    N_r = r_der()[2]

    M_matrix = np.array([
    [m, 0, 0, 0, 0],
    [0, I_x, -I_xz, 0, 0],
    [0, -I_xz, I_z, 0, 0],
    [0, 0, 0, 1, 0],
    [0, 0, 0, 0, 1]
    ])

    B_matrix = np.array([
        [-Y_beta, -(Y_p + m * W_e), -(Y_r - m * U_e), -m * g * np.cos(theta_e), -m * g * np.sin(theta_e)],
        [-L_beta, -L_p, -L_r, 0, 0],
        [-N_beta, -N_p, -N_r, 0, 0],
        [0, -1, 0, 0, 0],
        [0, 0, -1, 0, 0]
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
    
    eig = eigenvalues[2]
    sigma = eig.real
    omega_d = eig.imag
    omega_n = np.sqrt(sigma**2 + omega_d**2)
    zeta = -sigma / omega_n if omega_n != 0 else 0
    print(f"Dutch roll mode :")
    print(f"  ω_n = {omega_n:.4f} rad/s")
    print(f"  ζ = {zeta:.4f}\n")
    if zeta > 0.08 and zeta*omega_n > 0.15 and omega_n > 0.4: 
        print("You are level 1\n")
    if zeta > 0.02 and zeta*omega_n > 0.05 and omega_n > 0.4:
        print("You are level 2\n")
    time_constant_2 = -1/eigenvalues[1].real
    time_constant_1 = -1/eigenvalues[4].real
    print(f"Roll subsidence mode:")
    print(f"The time constant for the roll subsidence mode = {time_constant_1:.4f} s")
    if time_constant_1 > 20:
        print("You are level 1\n")
    if time_constant_1 > 8 and time_constant_1 < 20:
        print("You are level 2\n")
    if time_constant_1 > 4 and time_constant_1 < 8:
        print("You are level 3\n")
    print(f"Spiral mode:")
    print(f"The time constant for the spiral mode = {time_constant_2:.4f} s")
    if time_constant_2 < 1.4: 
        print("You are level 1\n")
    if time_constant_2 > 1.4 and time_constant_2 < 3:
        print("You are level 2\n")
    if time_constant_2 > 3 and time_constant_2 < 10:
        print("You are level 3\n")

    return real_parts, test, eigenvalues

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

def short_period():

    M_q = q_der2()[2]
    M_w = w_der()[2]
    Z_w = w_der()[0]
    M_w_dot = w_dot_der()[2]

    omega_s = np.sqrt((M_q/I_y)*(Z_w/m)-(M_w/I_y)*U_e)
    xi_s = -(M_q/I_y + Z_w/m + M_w_dot/I_y*U_e)/(2*omega_s)

    return omega_s, xi_s

##################################################################
######DEFINITION OF THE PRINT FUNCTION
##################################################################

def main(): 

    real_parts_long, test_long, eigenvalues_long = long_dyn_stab()
    real_parts_lat, test_lat, eigenvalues_lat = lat_dyn_stab()
    omega_p, xi_p = phugoid()
    omega_s, xi_s = short_period()

    print("--------------------------EIGENVALUES--------------------------------------------")
    print("The eigenvales for the longitudinal motion are: ", eigenvalues_long)
    print("The eigenvales for the lateral motion are: ", eigenvalues_lat)
    print("----------------------------------------------------------------------")

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
    
    if np.all(real_parts_lat < 0):
        print("The aircraft is dynamically stable for the lateral motion.")
    elif np.any(real_parts_lat > 0):
        print("The system is dynamiccaly unstable for the lateral motion.")
    elif np.any(real_parts_lat == 0):
        print("The system is dynamically neutrally stable for the lateral motion.")
    print("----------------------------------------------------------------------")

    print("--------------------------PHUGOID MODE--------------------------------------------")
    print("The frequency of the phugoid mode is: ", omega_p, "rad/s")
    print("The damping ratio of the phugoid mode is: ", xi_p)

    print("--------------------------SHORT PERIOD MODE--------------------------------------------")
    print("The frequency of the short period mode is: ", omega_s, "rad/s")
    print("The damping ratio of the short period mode is: ", xi_s)
    print("----------------------------------------------------------------------")

    print("--------------------------LONGITUDINAL DERIVATIVES--------------------------------------------")
    print("The longitudinal derivatives are:")
    print("CL_alpha: ", alpha_der()[0])
    print("CD_alpha: ", alpha_der()[1])
    print("CM_alpha: ", alpha_der()[2])
    print("CL_u: ", u_der1()[0])
    print("CD_u: ", u_der1()[1])
    print("CM_u: ", u_der1()[2])
    print("CL_q: ", q_der1()[0])
    print("CD_q: ", q_der1()[1])
    print("CM_q: ", q_der1()[2])
    print("CL_alpha_dot: ", alpha_dot_der()[0])
    print("CD_alpha_dot: ", alpha_dot_der()[1])
    print("CM_alpha_dot: ", alpha_dot_der()[2])
    print("")
    print("And, in the matrices, we have:")
    print("Z_u: ", u_der2()[0])
    print("X_u: ", u_der2()[1])
    print("M_u: ", u_der2()[2])
    print("Z_q: ", q_der2()[0])
    print("X_q: ", q_der2()[1])
    print("M_q: ", q_der2()[2])
    print("X_w_dot: ", w_dot_der()[0])
    print("Z_w_dot: ", w_dot_der()[1])
    print("M_w_dot: ", w_dot_der()[2])
    print("X_w: ", w_der()[0])
    print("Z_w: ", w_der()[1])
    print("M_w: ", w_der()[2])
    print("")
    print("----------------------------------------------------------------------")
    print("--------------------------LATERAL DERIVATIVES--------------------------------------------")
    print("The lateral derivatives are:")
    print("L_beta: ", beta_der()[0])
    print("N_beta: ", beta_der()[1])
    print("Y_beta: ", beta_der()[2])
    print("Y_r: ", r_der()[0])
    print("L_r: ", r_der()[1])
    print("N_r: ", r_der()[2])
    print("Y_p: ", p_der()[0])
    print("L_p: ", p_der()[1])
    print("N_p: ", p_der()[2])
    print("----------------------------------------------------------------------")


if __name__ == "__main__":
    main()
