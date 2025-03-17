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
alpha_e = 
M = 
delta = 0.005
Cl_tot0, Cd_tot0, cl_max, AoA_L0, Cl_tot, CD, AoA, cd0, CL_alfa = get_Lift_and_drag(AR, delta, sweep_LE_fus, sweep_LE_wing, force)

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

def u_der(i,d): 

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

##################################################################
######Q DERIVATIVES 
##################################################################

def q_der(i,d):

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