import numpy as np 
from static_stability import long_stat_stab_cruise
from wings import getClAlfa
from wings import getCl
from wings import get_Lift_and_drag
from tail import LiftCurveSlope
from tail import setting_angle
from static_stability import tail_eff
from tail import geomtail
from static_stability import CG_position
from wings import air_density
from wings import detSurfac
from wings import true_airspeed_at_altitude
from tail import surf_tail
from wings import getSurface_And_AR
from wings import wingGeometryIDEAL
from wings import getMAC
from wings import getAirfoilWing
from wings import getAirfoilFus
from static_stability import boucleForce

##################################################################
######CONFIGURATION SETTING
##################################################################

config = 3
fuel = 2

##################################################################
######INPUT PARAMETERS
##################################################################

sweep_LE_fus = 50.0
sweep_quarter_wing = 29.0
CL = 0.45
Cm0_airfoil_wing = getAirfoilWing()[4]
Cm0_airfoil_fus = getAirfoilFus()[4]
force = boucleForce(config,fuel,Cm0_airfoil_fus,Cm0_airfoil_wing, CL, sweep_LE_fus, sweep_quarter_wing)
e = 0.85 # Oswald efficiency factor
alpha_e = 0*np.pi/180
delta = 0.005
z_f = 1.51 #distance between AC of the tail and AC of the aircraft
l_cabin = 10.1
l_cockpit = 2.01

##################################################################
######GENERAL PARAMETERS
##################################################################

surface_wing_ideal, AR = getSurface_And_AR(CL, force)
M = 0.9
rho = air_density(12500)[0]
V0 = true_airspeed_at_altitude(12500,0.9)
surface_wing_ideal, AR, taper_wing, croot, c_tip_wing, chord_middle, sweep_LE_wing, sweep_beta = wingGeometryIDEAL(CL, force, sweep_quarter_wing)
surface_wing_ideal, surf_fus, surf_wing, surf_tot = detSurfac(CL, sweep_LE_fus, sweep_quarter_wing, force)
Cl_tot0, Cd_tot0, cl_max, AoA_L0, Cl_tot, CD, AoA, cd0, CL_alfa = get_Lift_and_drag(AR, delta, sweep_LE_fus, sweep_LE_wing, force)
MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(CL, sweep_LE_fus, sweep_quarter_wing, force)
surf_vert_tail = surf_tail()[1]
c_root_tail,span_hor_tail,span_vert_tail,AR_h_tail, AR_tail,surf_vert_tail, surf_tot_tail, MAC_tail,y_AC_tail,x_AC_tail_local = geomtail()
x_AC_tail = l_cabin + l_cockpit + x_AC_tail_local + 1
x_CG_tot = CG_position(config,fuel, AR, sweep_LE_fus, sweep_LE_wing)
l_f = x_AC_tail - x_CG_tot

##################################################################
######COMPUTATION OF GENERAL PARAMETERS
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

    Kn = long_stat_stab_cruise(i,d,AR, sweep_LE_fus, sweep_LE_wing)[0]
    CL_alpha = getClAlfa(AR, sweep_LE_fus, sweep_LE_wing)
    CD_alpha = 2*CL*CL_alpha/(np.pi*AR*e)
    CM_alpha = -Kn * CL_alpha

    return CL_alpha, CD_alpha, CM_alpha

##################################################################
######U DERIVATIVES 
##################################################################

def u_der1(i,d): 

    CL_alpha = alpha_der(i,d)[0]
    CM_alpha = alpha_der(i,d)[2]
    CL_M = 
    CD_alpha = alpha_der(i,d)[1]
    CD_M = 
    CM_alpha = alpha_der(i,d)[2]
    x_AC_M = 

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

    X_u = X_u * (1/2*rho*V0*surf_tot)

    return Z_u, X_u, M_u

##################################################################
######Q DERIVATIVES 
##################################################################

def q_der1(i,d):

    CL_alpha = alpha_der(i,d)[0]
    CL_alpha_tail = LiftCurveSlope()
    cl_alpha_wing = getAirfoilWing()[0]
    eta_tail = setting_angle(force)
    V_tail = tail_eff(i,d, AR, sweep_LE_fus, sweep_LE_wing)
    c_root_tail,span_hor,span_vert,AR_h, AR,gamma_h, surf_tot_tail, MAC_tail,yac,xac = geomtail()
    tail_pos = 1
    x_AC_tail = tail_pos+xac
    l_T = x_AC_tail - x_CG_tot
    X_W = np.abs(x_AC_tot - x_CG_tot)
    b = 28.95*3.28084
    mean_chord = AR/b

    CL_q_wing_M0 = (1/2 + 2*X_W/mean_chord)*CL_alpha
    B = np.sqrt(1 - M**2 * np.cos(sweep_quarter_wing)**2)
    CL_q_wing = (AR + 2*np.cos(sweep_quarter_wing))/(AR*B + 2*np.cos(sweep_quarter_wing))*CL_q_wing_M0
    CL_q_tail = 2*CL_alpha_tail*eta_tail * V_tail

    K = 0.7
    CM_q_wing_M0 = -K * cl_alpha_wing * np.cos(sweep_quarter_wing) * ((AR*(2 * (X_W/MAC_wing)**2 + 1/2 *(X_W/MAC_wing)))(AR + 2*np.cos(sweep_quarter_wing)) + 1/24 * (AR**3 * np.tan(sweep_quarter_wing)**2)/(AR + 6*np.cos(sweep_quarter_wing) + 1/8))
    CM_q_wing = (((AR**3 * np.tan(sweep_quarter_wing)**2)/(AR*B + 6*np.cos(sweep_quarter_wing)) + 3/B)/(((AR**3 * np.tan(sweep_quarter_wing)**2)/(AR*B + 6*np.cos(sweep_quarter_wing))) + 3)) * CM_q_wing_M0
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

    return Z_q, X_q, M_q

##################################################################
######ALPHA DOT DERIVATIVES 
##################################################################

def alpha_dot_der(i,d):

    CL_alpha_H = LiftCurveSlope()
    deps_dalpha = 0 #approximation
    eta_H = setting_angle(force)
    V_tail = tail_eff(i,d, AR, sweep_LE_fus, sweep_LE_wing)
    c_root_tail,span_hor,span_vert,AR_h, AR,gamma_h, surf_tot_tail, MAC_tail,yac,xac = geomtail()
    tail_pos = 1
    x_AC_tail = tail_pos+xac
    x_CG_tot = CG_position(i,d, AR, sweep_LE_fus, sweep_LE_wing)
    l_T = x_AC_tail - x_CG_tot
    CD_alpha_dot = 0

    CL_alpha_dot_W = 0
    CL_alpha_dot_H = 2*CL_alpha_H*eta_H*V_tail*deps_dalpha
    CL_alpha_dot = CL_alpha_dot_W + CL_alpha_dot_H

    CM_alpha_dot_W = 0
    CM_alpha_dot_H = -2*CL_alpha_H*eta_H*V_tail*deps_dalpha*l_T/MAC_tail
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

    M_w_dot = M_w_dot * (1/2*rho*surf_tot*MAC_tot**2)

    return X_w_dot, Z_w_dot, M_w_dot

##################################################################
######W DERIVATIVES 
##################################################################

def w_der(i,d):

    CL_alpha = alpha_der(i,d)[0]
    CD_alpha = alpha_der(i,d)[1]
    CM_alpha = alpha_der(i,d)[2]
    Le = force
    lift_coef, CD, CL_max, AoA_L0, cl, _, aoa, Cd_tot0, CL_alfa = get_Lift_and_drag(CL, delta, sweep_LE_fus, sweep_quarter_wing, force)
    De = CD*1/2*rho*V0**2*surface_wing_ideal
    Xe = Le*np.sin(alpha_e) - De*np.cos(alpha_e)
    Ze = -Le*np.cos(alpha_e) - De*np.sin(alpha_e)

    C_Xe = Xe/(1/2*rho*V0**2*surf_tot)
    C_Ze = Ze/(1/2*rho*V0**2*surf_tot)

    Z_w = 1/np.cos(alpha_e) * (C_Xe -CL_alpha*np.cos(alpha_e) - CD_alpha*np.sin(alpha_e))
    X_w = 1/np.cos(alpha_e) * (-C_Ze + CL_alpha*np.sin(alpha_e) - CD_alpha*np.cos(alpha_e))
    M_w = 1/np.cos(alpha_e) * CM_alpha

    return Z_w, X_w, M_w

##################################################################
######CONSTRUCTION OF THE A MATRIX FOR THE LONGITUDINAL MOTION
##################################################################

def long_dyn_stab(i,d,Cl): 

    position, pourc_wings, motors_pos, total_weight = CG_position(i,d, Cl, sweep_LE_fus, sweep_quarter_wing, force)
    m = 
    g = 9.81
    I_y = 

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

    surf_vert_tail = surf_tail()[1]
    average_fus_sect = 
    d = np.sqrt(average_fus_sect/0.7854)
    Z_w = #vertical distance from the wing root quarter-chord to the fuselage centerline 
    d_sigma = -0.276 + 3.06*surf_vert_tail/surf_tot /(1 + np.cos(sweep_quarter_wing)) + 0.4*Z_w + 0.009*AR
    return d_sigma

##################################################################
######V DERIVATIVES 
##################################################################

def v_der(sweep_quarter_wing,AR):

    d_sigma = dsigma_dbeta(sweep_quarter_wing,AR)
    c_1 = #fin lift curve slope
    Y_v = 1/2*rho*V0**2 *surf_vert_tail * c_1 * (1-d_sigma)
    L_v = -1/2 * rho*V0**2 * surf_vert_tail * z_f * c_1 * (1-d_sigma)
    N_v = 1/2 * rho*V0*surf_vert_tail * l_f * c_1 * (1-d_sigma)

    return Y_v, L_v, N_v

##################################################################
######R DERIVATIVES 
##################################################################

def r_der(sweep_quarter_wing,AR):

    b = 28.95*3.28084
    Y_v = v_der(sweep_quarter_wing,AR)[0]
    Y_r = -2*Y_v * (z_f * np.sin(alpha_e) - l_f * np.cos(alpha_e)/b)
    L_r = -2*Y_v*((z_f*np.sin(alpha_e) + l_f*np.cos(alpha_e)*(z_f*np.cos(alpha_e) - l_f*np.sin(alpha_e)))/b**2)
    N_r = 2*Y_v * ((z_f*np.sin(alpha_e + l_f*np.cos(alpha_e))**2)/b**2)

    return Y_r, L_r, N_r

##################################################################
######P DERIVATIVES 
##################################################################

def p_der(sweep_quarter_wing,AR):
    
    b = 28.95*3.28084
    Y_v = v_der(sweep_quarter_wing,AR)[0]
    L_r = r_der(sweep_quarter_wing,AR)[1]
    Y_p = 2*Y_v * (z_f * np.cos(alpha_e) - l_f * np.sin(alpha_e)/b)
    L_p = 2 * Y_v * (z_f/b)**2
    N_p = L_r

    return Y_p, L_p, N_p


##################################################################
######CONSTRUCTION OF THE A MATRIX FOR THE LATERAL MOTION
##################################################################


def lat_dyn_stab(sweep_quarter_wing, AR): 

    m = 
    g = 9.81
    I_x = 
    I_xz = 
    I_z = 

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
######DEFINITION OF THE PRINT FUNCTION
##################################################################

def main(): 

    real_parts_long, test_long = long_dyn_stab(config,fuel,Cl_tot)
    real_parts_lat, test_lat = lat_dyn_stab(sweep_quarter_wing, AR)

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

if __name__ == "__main__":
    main()
