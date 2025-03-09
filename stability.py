import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import trapz
from wings import fusGeometry as fusgeom
from wings import wingGeometry as winggeom
from sympy import symbols, Eq, solve
from wings import air_density
from wings import true_airspeed_at_altitude
from weight import get_weight
from wings import wingFuelvolume
from wings import getMAC
from wings import fusGeometry
from wings import detSurfac
from wings import wingCL
from wings import fuselageCL
from wings import getAirfoilFus
from wings import getAirfoilWing
from wings import getClAlfa
from wings import wingCL
from wings import wingGeometry
from tail import LiftCurveSlope
from tail import geomtail
from tail import surfhor_tail
from tail import LiftCurveSlope


##################################################################
######GENERAL QUANTTITES
##################################################################

rho = air_density(12500)[0]
speed = true_airspeed_at_altitude(12500)

delta = 0.005
b = 28.95
l_fus = 16.8
l_cabin = 10.1
l_cockpit = 2.01 
span_wings = 20 
l_aft = l_fus - l_cabin - l_cockpit
wings_taper_ratio = 0.4

##################################################################
######COMPUTATION OF THE CM0_TOT
##################################################################

Cm0_wing = getAirfoilWing()[4]
Cm0_fus = getAirfoilFus()[4]

def Cm0(Cm0_airfoil_fus,Cm0_airfoil_wing, AR ,sweep_LE_fus, sweep_LE_wing):
    _, _, _, _, _, _, _, _, y_wing, leading_wing, trailing_wing, quarter_wing,_,_ = winggeom(AR,sweep_LE_fus, sweep_LE_wing) 
    _, _, _, _, _, _, _, y_fus, leading_fus, trailing_fus, quarter_fus = fusgeom(AR,sweep_LE_fus) 
    surf_tot,surf_fus,surf_wing = detSurfac(AR, sweep_LE_fus)
    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(AR, sweep_LE_fus, sweep_LE_wing)

    c_fus = trailing_fus - leading_fus
    c_wing = trailing_wing - leading_wing
    Cm0_wing = (2/(surf_wing*MAC_wing)) * trapz(Cm0_airfoil_wing*c_wing**2, y_wing)
    Cm0_fus = (2/(surf_fus*MAC_fus)) * trapz(Cm0_airfoil_fus*c_fus**2, y_fus)
    Cm0_tot = (Cm0_wing*surf_wing + Cm0_fus*surf_fus)/surf_tot
    return Cm0_tot,Cm0_fus,Cm0_wing


##################################################################
######GEOMETRICAL POSITIONS OF IMPORTANT POINTS
##################################################################

z_AC_tot = 0
z_CG_tot = 0
z_CG_motors = 2

##################################################################
######TAIL QUANTITIES
##################################################################

c_root_tail,span_hor_tail,span_vert_tail,AR_h_tail, AR_tail,surf_vert_tail, surf_tot_tail, MAC_tail,y_AC_tail,x_AC_tail_local = geomtail()

hor_tail_surf = surfhor_tail()
a1 = LiftCurveSlope()
x_AC_tail = l_cabin + l_cockpit + x_AC_tail_local +2
l_tail = MAC_tail

z_AC_tail = 1.51

##################################################################
######CONFIGURATION SETTING
##################################################################

config = 3
fuel = 1

##################################################################
######CG POSITION
##################################################################

def passengers(i): 

    passengers_weight = 8*215
    if i == 1 : #if no passengers 
       passengers_pos = 1 #arbitrary
       passengers_weight = 0
    if i == 2 : #if passengers are as close as possible to the nose
       passengers_pos = l_cockpit
    if i == 3 : #if passengers are as far as possible from the nose 
       passengers_pos = l_cockpit + l_cabin
    
    return passengers_weight, passengers_pos

def CG_position(i,d, AR, sweep_LE_fus, sweep_LE_wing): 
    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(AR, sweep_LE_fus, sweep_LE_wing)
    fus_weight, aft_weight, wing_weight, land_gear_weight,motors_weight,nacelle_weight,APU_weight,enginst_weight,instr_weight,hydr_syst_weight,furn_weight,air_cond_weight,payload_weight, ops_weight,elec_syst_weight,surf_cont_weight,_,_,_= get_weight()
    _,_,_,_,_,_,chord_tip_fus,_,_,_,_ = fusGeometry(AR, sweep_LE_fus)
    #force = CL(config,fuel,Cm0_fus,Cm0_wing, AR, sweep_LE_fus, sweep_LE_wing)[0] ## Boucle infinie en utilisant ceci

    sweep_angle_wing = sweep_LE_wing*np.pi/180
    sweep_angle_fus = sweep_LE_fus*np.pi/180
    pilots_weight = 220.462*2
    pilots_pos = l_cockpit*0.35

    wing_pos = (l_fus - chord_tip_fus) + MAC_wing*0.24 + y_AC_wing*np.tan(sweep_angle_wing)
    
    fus_pos = (l_cabin + l_cockpit)*0.465 #0.2*MAC_fus + y_AC_fus*np.tan(sweep_angle_fus)

    aft_pos = l_cockpit + l_cabin + l_aft*0.3

    APU_pos = l_cockpit + l_cabin

    hydr_pos = 0.75*wing_pos + 0.25 *x_AC_tail

    payload_pos = 0.35*(l_cockpit + l_cabin) + 0.65*l_cockpit

    ops_pos = l_cockpit 

    land_gear_pos = 9.703*0.7 + 0.3*0.673
    surf_cont_pos = (l_fus - chord_tip_fus) + MAC_wing*0.35 + y_AC_wing*np.tan(sweep_angle_wing)

    instr_pos = 0.35*l_cockpit

    furn_pos = 0.4*l_fus

    air_cond_pos = l_cockpit-0.2
    
    motors_pos = l_cockpit + l_cabin + 0.6

    nacelle_pos = motors_pos

    enginst_pos = motors_pos

    elec_syst_pos = 0.75*l_fus/2 + 0.25*motors_pos
  
    available_fuel_vol = 5 * 1000 #wingFuelvolume(AR, sweep_LE_fus, sweep_LE_wing, force)*1000 
    if d == 1 : #no fuel
        fuel_weight = 0
        fuel_pos = 1
        pourc_wings = 0
    if d == 2 : #full of fuel
        vol_fuel = 27300
        fuel_weight = vol_fuel*0.8*2.20462
        pourc_wings = available_fuel_vol/vol_fuel
        fuel_pos = wing_pos*pourc_wings + (1-pourc_wings)*l_fus*0.44

    passengers_weight = passengers(i)[0]
    passengers_pos = passengers(i)[1] 

    """
    print("fuel:",fuel_pos,"m and",fuel_pos*3.28084,"ft ->",fuel_pos/l_fus*100)
    print("wing:",wing_pos, "m and",wing_pos*3.28084,"ft ->", wing_pos/l_fus *100)
    print("fus:",fus_pos, "m and",fus_pos*3.28084,"ft ->", fus_pos/l_fus *100)
    print("aft:",aft_pos, "m and",aft_pos*3.28084,"ft ->", aft_pos/l_fus *100)
    print("APU:",APU_pos, "m and", APU_pos*3.28084,"ft ->", APU_pos/l_fus *100)
    print("hydr:",hydr_pos, "m and", hydr_pos*3.28084,"ft ->", hydr_pos/l_fus *100)
    print("payload:",payload_pos, "m and",payload_pos*3.28084,"ft ->", payload_pos/l_fus *100) 
    print("ops:",ops_pos, "m and",ops_pos*3.28084,"ft ->", ops_pos/l_fus *100)
    print("land_gear:",land_gear_pos, "m and", land_gear_pos*3.28084,"ft ->", land_gear_pos/l_fus *100)
    print("surf_cont:",surf_cont_pos, "m and", surf_cont_pos*3.28084,"ft ->", surf_cont_pos/l_fus *100)
    print("instr:",instr_pos, "m and", instr_pos*3.28084,"ft ->", instr_pos/l_fus *100)
    print("furn:",furn_pos, "m and",furn_pos*3.28084,"ft ->", furn_pos/l_fus *100)
    print("air_cond:",air_cond_pos, "m and",air_cond_pos*3.28084,"ft ->", air_cond_pos/l_fus *100)
    print("motors:",motors_pos, "m and",motors_pos*3.28084,"ft ->", motors_pos/l_fus *100)
    print("nacelle:",nacelle_pos, "m and", nacelle_pos*3.28084,"ft ->", nacelle_pos/l_fus *100)
    print("enginst:",enginst_pos, "m and",enginst_pos*3.28084,"ft ->", enginst_pos/l_fus *100)
    print("elec_syst:",elec_syst_pos, "m and",elec_syst_pos*3.28084,"ft ->", elec_syst_pos/l_fus *100)
    """
    
    total_mom = (wing_weight*wing_pos) + (fus_weight*fus_pos) + (land_gear_weight*land_gear_pos) + (surf_cont_weight*surf_cont_pos) + (instr_weight*instr_pos) + (elec_syst_weight*elec_syst_pos) + (furn_weight*furn_pos) + (air_cond_weight*air_cond_pos) + (passengers_weight*passengers_pos) + (motors_weight*motors_pos) + (fuel_pos*fuel_weight) + (aft_pos*aft_weight) + (nacelle_pos*nacelle_weight) + (APU_pos*APU_weight) + (enginst_pos*enginst_weight) + (hydr_pos*hydr_syst_weight) + (payload_pos*payload_weight) + (ops_pos*ops_weight) + (pilots_weight*pilots_pos)
    total_weight = wing_weight + fus_weight + land_gear_weight + surf_cont_weight + instr_weight + elec_syst_weight + furn_weight + air_cond_weight + passengers_weight + motors_weight + fuel_weight + aft_weight + nacelle_weight + APU_weight + enginst_weight + hydr_syst_weight + payload_weight + ops_weight + pilots_weight
    position = total_mom/total_weight

    return position, pourc_wings, motors_pos/MAC_tot, total_weight 


##################################################################
######EQUILIBRIUM IN PITCH
##################################################################


#tail volume ratio effectivness 
def tail_eff(i,d, AR, sweep_LE_fus, sweep_LE_wing):
    surf_tot,surf_fus,surf_wing = detSurfac(AR, sweep_LE_fus)
    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(AR, sweep_LE_fus, sweep_LE_wing)

    x_CG_tot = CG_position(i,d, AR, sweep_LE_fus, sweep_LE_wing)[0]
    V_T = hor_tail_surf * (x_AC_tail - x_CG_tot)/(surf_tot* MAC_tot)
    return V_T

def prop_force(): 
    A = 2*1.32**2 *np.pi
    m_dot = rho*A*speed*0.453592/32.2
    Fp = m_dot*speed
    return Fp

def wings_lift(AR ,sweep_LE_fus, sweep_LE_wing):
    surf_tot,surf_fus,surf_wing = detSurfac(AR, sweep_LE_fus)
    force = CL(config,fuel,Cm0_fus,Cm0_wing, AR, sweep_LE_fus, sweep_LE_wing)[0]
    _,CL_w0,_,_,_,_ = wingCL(AR ,sweep_LE_fus, sweep_LE_wing, force)

    L_w = CL_w0 * (1/2)*rho*surf_wing*speed**2
    return L_w

def fus_lift(AR ,sweep_LE_fus, sweep_LE_wing):
    surf_tot,surf_fus,surf_wing = detSurfac(AR, sweep_LE_fus) 
    _,CL_f0,_,_,_ = fuselageCL(AR, sweep_LE_fus)

    L_f = CL_f0*(1/2)*rho*surf_fus*speed**2
    return L_f


#i : configuration numerotation
def CL(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing, AR, sweep_LE_fus, sweep_LE_wing): 
    surf_tot,surf_fus,surf_wing = detSurfac(AR, sweep_LE_fus)
    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(AR, sweep_LE_fus, sweep_LE_wing)
    L_tot, L_T = symbols('L_tot L_T')
    V_T = tail_eff(i,d, AR, sweep_LE_fus, sweep_LE_wing)
    T = 49734.78
    x_CG_tot = CG_position(i,d, AR, sweep_LE_fus, sweep_LE_wing)[0]
    weight = CG_position(i,d, AR, sweep_LE_fus, sweep_LE_wing)[3]*9.81*0.453592 + passengers(i)[0]*9.81*0.453592
    Cm0_tot = Cm0(Cm0_airfoil_fus,Cm0_airfoil_wing, AR ,sweep_LE_fus, sweep_LE_wing)[0]
    M0 = Cm0_tot*((1/2)*rho*(speed**2)*MAC_tot*surf_tot)
    x_CG_motors = CG_position(i,d, AR, sweep_LE_fus, sweep_LE_wing)[2]*MAC_tot
    Fp = prop_force()

    #translation equilibrium 
    eq1 = Eq(L_tot + L_T - weight + Fp,0)
    
    #rotation equilibrium 
    eq2 = Eq(M0 + L_tot*(x_CG_tot-x_AC_tot)-L_T*(x_AC_tail - x_CG_tot)-T*(z_CG_motors - z_CG_tot) + Fp*(x_CG_tot-x_CG_motors),0) 

    solution = solve((eq1,eq2),(L_tot,L_T))
    L_tot = solution[L_tot]
    L_T = solution[L_T]

    return L_tot,L_T


##################################################################
######LONGITUDINAL STABILITY
##################################################################


def long_stat_stab_cruise(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing, AR, sweep_LE_fus, sweep_LE_wing): #in the pitching plane
    surf_tot,surf_fus,surf_wing = detSurfac(AR, sweep_LE_fus)
    a = getClAlfa(AR, sweep_LE_fus, sweep_LE_wing)
    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(AR, sweep_LE_fus, sweep_LE_wing)
    #check the stability
    #neutral point : position of the cg in order to have the derivative equals 0
    x_CG_tot = CG_position(i,d, AR, sweep_LE_fus, sweep_LE_wing)[0]
    engines_pos = CG_position(i,d, AR, sweep_LE_fus, sweep_LE_wing)[2]
    deps = 0
    Fp = prop_force()
    eta = 0.9
    #hn = x_AC_tot/MAC_tot + V_T*a1_over_a*(1- deps) + dalpha_prop * Fp * (engines_pos-x_CG_tot)/MAC_tot#- 0.5*fus_width**2 * fus_length/(S_wing*a*MAC_wing)  #position of the neutral point  

    #page 416 Raymer
    q = 1/2*rho*speed**2
    hn = (a * x_AC_tot/MAC_tot +eta*a1*hor_tail_surf/surf_tot*(1-deps)*x_AC_tail/MAC_tot+Fp/(q*surf_tot)*(1-deps)*engines_pos)/(a + eta*hor_tail_surf/surf_tot*a1*(1-deps)+Fp/(q*surf_tot))

    Kn = hn - x_CG_tot/MAC_tot 
    return Kn, hn


##################################################################
######CG POSITIONS RANGE
##################################################################


def get_CG(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing,Kn, AR, sweep_LE_fus, sweep_LE_wing): 
    surf_tot,surf_fus,surf_wing = detSurfac(AR, sweep_LE_fus)
    a = getClAlfa(AR, sweep_LE_fus, sweep_LE_wing)
    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(AR, sweep_LE_fus, sweep_LE_wing)
    engines_pos = CG_position(i,d, AR, sweep_LE_fus, sweep_LE_wing)[2]
    deps = 0
    V_T = tail_eff(i,d, AR, sweep_LE_fus, sweep_LE_wing)
    dalpha_prop = 1 - deps
    Fp = prop_force()
    eta = 0.9
    q = 1/2*rho*speed**2
    hn = (a * x_AC_tot/MAC_tot + eta*a1*hor_tail_surf/surf_tot*(1-deps)*x_AC_tail/MAC_tot+Fp/(q*surf_tot)*(1-deps)*engines_pos)/(a + eta*hor_tail_surf/surf_tot*a1*(1-deps)+Fp/(q*surf_tot))

    #x_CG = (hn*MAC_tot) - (Kn*MAC_tot)
    x_CG1 = (hn*MAC_tot) - (0.05*MAC_tot)
    x_CG2 = (hn*MAC_tot) - (0.15*MAC_tot)
    return x_CG1, x_CG2


##################################################################
######DIRECTIONAL STABILITY
##################################################################


def interpolation(x1, y1, x2, y2, x3):
    t=(x3-x1)/(x2-x1) 
    y3 = y1 + t * (y2 - y1)
    
    return y3

def dir_stat_stab_cruise(CG_position,AR, sweep_LE_fus):  
    surf_tot,surf_fus,surf_wing = detSurfac(AR, sweep_LE_fus)
    hf1 = (interpolation(0.2481847, 0.129909, 0.2632646, 0.1303305, 0.25)+interpolation(0.2491559,0.04703807 , 0.2613637,0.0475382 , 0.25))*l_fus  # forward fuselage height
    hf2 = (interpolation(0.7477641, 0.04391509, 0.7610847, 0.04077155, 0.75)+interpolation(0.7403715, 0.05064632, 0.7543387, 0.04952601, 0.75))*l_fus # rear fuselage height
    bf1 = 5.881743321 # forward fuselage width
    bf2 = 9  # rear fuselage width
    hf_max = 0.179*16.8  # maximum fuselage height
    x_CG = CG_position 
    L_f=6  # distance between center of gravity and aerodynamic center of the tail
     
    K_beta = 0.3 * (x_CG / l_fus) + 0.75 * (hf_max / l_fus) - 0.105
    CN_beta_fuselage = -K_beta * (surf_fus*l_fus/(surf_wing*span_wings))*((hf1/hf2)**0.5)*((bf2/bf1)**(1/3))
    
    CN_beta_w=0.012 #{High, mid, low}-mounted wing effect = {-0.017,0.012,0.024}
    
    a = LiftCurveSlope()
    c_root_tail,span_hor,span_vert,AR_h, AR,surf_vert_tail, surf_tot_tail, MAC_tail,yac_wing,xac_wing = geomtail()
    CN_beta_fin=a*surf_vert_tail*L_f/(surf_wing*span_wings)
    
    CN_beta_tot=CN_beta_fin+CN_beta_w+CN_beta_fuselage
    
    return CN_beta_fin, CN_beta_fuselage, CN_beta_w, CN_beta_tot


##################################################################
######LATERAL STABILITY
##################################################################

def lat_stat_stab_cruise(dihedral_angle,AR,sweep_LE_fus, sweep_LE_wing): 

    gamma = dihedral_angle*np.pi/180
    weight = CL(config,fuel,Cm0_fus,Cm0_wing, AR, sweep_LE_fus, sweep_LE_wing)[0]
    CL_wings, CL_w0, CD_wing, CL_max, alpha_L0, CL_alpha_wings = wingCL(AR,sweep_LE_fus, sweep_LE_wing, weight)

    CL_beta_dihedral = -0.25 *CL_alpha_wings * gamma* (2 * (1 + 2*wings_taper_ratio)/(3*(1+wings_taper_ratio)))

    graph_value = 0.35 #value from Nicolai's book page 590
    CL_beta_wing_sweep = graph_value * CL_wings
    CL_beta_wings = CL_beta_wing_sweep + CL_beta_dihedral

    CL_beta_wings_fus = 0 #because middle mounted wing

    CL_beta_tot = CL_beta_wings_fus + CL_beta_wings

    return CL_beta_tot, CL_beta_wings, CL_beta_wings_fus

##################################################################
######PRINT OF STABILITY
##################################################################


def printFunction(AR, sweep_LE_fus, sweep_LE_wing, dihedral_angle):

    MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC(AR, sweep_LE_fus, sweep_LE_wing)
    CN_beta_fin, CN_beta_fuselage, CN_beta_w, CN_beta_tot = dir_stat_stab_cruise(CG_position(config,fuel, AR, sweep_LE_fus, sweep_LE_wing)[0],AR, sweep_LE_fus)
    CL_beta_tot, CL_beta_wings, CL_beta_wings_fus = lat_stat_stab_cruise(dihedral_angle,AR,sweep_LE_fus, sweep_LE_wing)

    print("--------------------------AERODYNAMIC CENTER--------------------------------------------")
    print("The aerodynamic center is positioned at",x_AC_tot,"from the nose of the airplane, which represents",x_AC_tot*100/l_fus,"% of the total length.")
    print("----------------------------------------------------------------------")

    print("--------------------------CM0--------------------------------------------")
    print("The blended wing body has a Cm0_tot of",Cm0(Cm0_fus,Cm0_wing, AR ,sweep_LE_fus, sweep_LE_wing)[0])
    print("The wings have a pitching moment (up) of",Cm0(Cm0_fus,Cm0_wing, AR ,sweep_LE_fus, sweep_LE_wing)[2],". It represents",Cm0(Cm0_fus,Cm0_wing, AR ,sweep_LE_fus, sweep_LE_wing)[2]/Cm0(Cm0_fus,Cm0_wing, AR ,sweep_LE_fus, sweep_LE_wing)[0]*100,"% of the pitching moment of the airplane.")
    print("The fuselage has a pitching moment (up) of",Cm0(Cm0_fus,Cm0_wing, AR ,sweep_LE_fus, sweep_LE_wing)[1],". It represents",Cm0(Cm0_fus,Cm0_wing, AR ,sweep_LE_fus, sweep_LE_wing)[1]/Cm0(Cm0_fus,Cm0_wing, AR ,sweep_LE_fus, sweep_LE_wing)[0]*100,"% of the pitching moment of the airplane.")
    print("----------------------------------------------------------------------")

    print("--------------------------TAIL--------------------------------------------")
    print("The position of the aerodynamic center of the tail is at",x_AC_tail," m and the end of the tail is at",x_AC_tail- x_AC_tail_local + l_tail,"m from the nose and the distance between x_AC_tail and x_CG_tot is",x_AC_tail - CG_position(config,fuel, AR, sweep_LE_fus, sweep_LE_wing)[0],"m.")
    print("----------------------------------------------------------------------")

    print("--------------------------FUEL STORAGE--------------------------------------------")
    print(CG_position(config,fuel, AR, sweep_LE_fus, sweep_LE_wing)[1]*100,"% of the total fuel volume is stored in the wings and we need",(1-CG_position(config,fuel, AR, sweep_LE_fus, sweep_LE_wing)[1])*26854.56/1000,"m³ in the fuselage.")
    print("----------------------------------------------------------------------")

    print("--------------------------CENTER OF GRAVITY--------------------------------------------")
    print("The center of gravity is positioned at",CG_position(config,fuel, AR, sweep_LE_fus, sweep_LE_wing)[0],"m (in feet :",CG_position(config,fuel, AR, sweep_LE_fus, sweep_LE_wing)[0]*3.28," (",CG_position(config,fuel, AR, sweep_LE_fus, sweep_LE_wing)[0]*100/l_fus,"%) from the nose of the airplane.")
    print("----------------------------------------------------------------------")

    print("--------------------------CRUISE--------------------------------------------")
    print("The new lift force generated by the body is",CL(config,fuel,Cm0_fus,Cm0_wing, AR, sweep_LE_fus, sweep_LE_wing)[0],"[N]")
    print("The new lift force generated by the tail is",CL(config,fuel,Cm0_fus,Cm0_wing, AR, sweep_LE_fus, sweep_LE_wing)[1],"[N]")
    print("----------------------------------------------------------------------")

    print("--------------------------STATIC MARGIN AND NEUTRAL POINT--------------------------------------------")
    Kn, hn = long_stat_stab_cruise(config,fuel,Cm0_fus,Cm0_wing, AR, sweep_LE_fus, sweep_LE_wing)
    if Kn >= 0.05 and Kn < 0.3 : 
        print("The static margin has a correct value and is equal to : ", (Kn*100), "%  and the neutral point is positioned at",hn*MAC_tot,"from the nose, which represents",hn*MAC_tot*100/l_fus,"% of the total length.")
    else : 
        print("The static margin has to be changed and is equal to : ", (Kn*100), "%.  and the neutral point is positioned at",hn*MAC_tot,"m from the nose, which represents",hn*MAC_tot*100/l_fus,"% of the total length.")
    print("----------------------------------------------------------------------")

    print("--------------------------ACCEPTABLE POSITIONS FOR THE CENTER OF GRAVITY--------------------------------------------")
    print("The center of gravity is positioned at",get_CG(config,fuel,Cm0_fus,Cm0_wing,0.05, AR, sweep_LE_fus, sweep_LE_wing)[0],"m from the nose when the static margin is equal to",0.05*100,"%. It is the maximal value of the range.")
    print("The center of gravity is positioned at",get_CG(config,fuel,Cm0_fus,Cm0_wing,0.05, AR, sweep_LE_fus, sweep_LE_wing)[1],"m from the nose when the static margin is equal to",0.15*100,"%. It is the minimal value of the range.")
    print("----------------------------------------------------------------------")

    print("--------------------------DIRECTIONAL STABILITY--------------------------------------------")
    if (CN_beta_tot<0.1):
        print ("The aircraft is not directionally stable because CN_beta_tot is negative and equals",CN_beta_tot)
        print("The contribution of the fuselage is",CN_beta_fuselage)
        print("The contribution of the wings is",CN_beta_w)
        print("The contribution of the fin is",CN_beta_fin)
    else:
        print ("The aircraft is directionally stable because CN_beta_tot is positive and equals",CN_beta_tot)
        print("The contribution of the fuselage is",CN_beta_fuselage)
        print("The contribution of the wings is",CN_beta_w)
        print("The contribution of the fin is",CN_beta_fin)
    print("----------------------------------------------------------------------")

    print("--------------------------LATERAL STABILITY--------------------------------------------")
    print("The CL_beta coefficient is equal to",CL_beta_tot)
    print("The contribution of the wings is",CL_beta_wings)
    print("----------------------------------------------------------------------")
    return

printFunction(3.4, 42, 25,3)