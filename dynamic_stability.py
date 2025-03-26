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

##################################################################
######CONFIGURATION SETTING
##################################################################

config = 3
fuel = 2

##################################################################
######GENERAL PARAMETERS
##################################################################

AR = 
sweep_LE_fus =
sweep_LE_wing =
force =
CL = getCl(AR, sweep_LE_fus, force)
e = 0.85 # Oswald efficiency factor
alpha_e = 2*np.pi/180
M = 
Cl_tot = 
delta = 0.005
rho = air_density(12500)[0]
V0 = true_airspeed_at_altitude(12500,0.9)
surface_wing_ideal, surf_fus, surf_wing, surf_tot = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, force)
Cl_tot0, Cd_tot0, cl_max, AoA_L0, Cl_tot, CD, AoA, cd0, CL_alfa = get_Lift_and_drag(AR, delta, sweep_LE_fus, sweep_LE_wing, force)

##################################################################
######COMPUTATION OF GENERAK PARAMETERS
##################################################################

U_e = V0 * np.cos(alpha_e)
W_e = V0 * np.sin(alpha_e)

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
    CL_M = CM_alpha/CL_alpha
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
    X_u = CL_u *np.sin(alpha_e) + CD_u* np.cos(alpha_e)
    M_u = CM_u

    return Z_u, X_u, M_u

##################################################################
######Q DERIVATIVES 
##################################################################

def q_der1(i,d):

    CL_alpha = alpha_der(i,d)[0]
    CL_alpha_tail = LiftCurveSlope()
    eta_tail = setting_angle(force)
    V_tail = tail_eff(i,d, AR, sweep_LE_fus, sweep_LE_wing)
    c_root_tail,span_hor,span_vert,AR_h, AR,gamma_h, surf_tot_tail, MAC_tail,yac,xac = geomtail()
    tail_pos = 1
    x_AC_tail = tail_pos+xac
    x_CG_tot = CG_position(i,d, AR, sweep_LE_fus, sweep_LE_wing)
    l_T = x_AC_tail - x_CG_tot

    CL_q_wing_M0 = (1/2 + 2*X_W/mean_chord)*CL_alpha
    B = np.sqrt(1 - M**2 * np.cos(sweep_quarter)**2)
    CL_q_wing = (AR + 2*np.cos(sweep_quarter))/(AR*B + 2*np.cos(sweep_quarter))*CL_q_wing_M0
    CL_q_tail = 2*CL_alpha_tail*eta_tail * V_tail

    K = 0.7
    CM_q_wing_M0 = -K * cl_alpha_wing * np.cos(sweep_quarter) * ((AR*(2 * (X_W/MAC_wing)**2 + 1/2 *(X_W/MAC_wing)))(A + 2*np.cos(sweep_quarter)) + 1/24 * (AR**3 * np.tan(sweep_quarter)**2)/(AR + 26*np.cos(sweep_quarter) + 1/8))
    CM_q_wing = (((AR**3 * np.tan(sweep_quarter)**2)/(AR*B + 6*np.cos(sweep_quarter)) + 3/B)/(((AR**3 * np.tan(sweep_quarter)**2)/(AR*B + 6*np.cos(sweep_quarter)) + 3/B) + 3)) * CM_q_wing_M0
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

    return X_w_dot, Z_w_dot, M_w_dot

##################################################################
######W DERIVATIVES 
##################################################################

def w_der(i,d):

    CL_alpha = alpha_der(i,d)[0]
    CD_alpha = alpha_der(i,d)[1]
    CM_alpha = alpha_der(i,d)[2]

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
######DEFINITION OF THE PRINT FUNCTION
##################################################################

def main(): 

    real_parts, test = long_dyn_stab(config,fuel,Cl_tot)

    print("--------------------------CHECK OF THE SYSTEM--------------------------------------------")
    if test == 1:
        print("The complex conjugate of all the eigenvalues are also eigenvalues. Thus, the system is well defined.")
    else:
        print("There is a problem with the definition of the system. There exists at least one complex conjugate of an eigenvalue which is not an eigenvalue.")
    print("----------------------------------------------------------------------")

    print("--------------------------STABILITY--------------------------------------------")
    if np.all(real_parts < 0):
        print("The aircraft is dynamically stable.")
    elif np.any(real_parts > 0):
        print("The system is dynamiccaly unstable.")
    elif np.any(real_parts == 0):
        print("The system is dynamically neutrally stable.")
    print("----------------------------------------------------------------------")

if __name__ == "__main__":
    main()

##################################################################
######LATERAL DYNAMIC DERIVATIVES
##################################################################

def dsigma_dbeta(sweep_quarter,AR): 

    surf_vert_tail = surf_tail()[1]
    average_fus_sect = 
    d = np.sqrt(average_fus_sect/0.7854)
    Z_w = #vertical distance from the wing root quarter-chord to the fuselage centerline 
    d_sigma = -0.276 + 3.06*surf_vert_tail/surf_tot /(1 + np.cos(sweep_quarter)) + 0.4*Z_w + 0.009*AR
    return d_sigma

##################################################################
######V DERIVATIVES 
##################################################################

def v_der(sweep_quarter,AR):
    d_sigma = dsigma_dbeta(sweep_quarter,AR)
    c_1 = #fin lift curve slope
    Y_v = 1/2*rho*speed**2 *surf_vert_tail * c_1 * (1-d_sigma)

    return Y_v

##################################################################
######P DERIVATIVES 
##################################################################

def p_der(sweep_quarter,AR):
    Y_v = v_der(sweep_quarter,AR)[0]
    Y_p = 2*Y_v * (z_f * np.cos(alpha_e) - l_f * np.sin(alpha_e)/b)
    return

##################################################################
######R DERIVATIVES 
##################################################################

def r_der(sweep_quarter,AR):
    Y_v = v_der(sweep_quarter,AR)[0]
    Y_r = -2*Y_v * (z_f * np.sin(alpha_e) - l_f * np.cos(alpha_e)/b)
    return

##################################################################
######CONSTRUCTION OF THE A MATRIX FOR THE LATERAL MOTION
##################################################################


def lat_dyn_stab(): 

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

    return 