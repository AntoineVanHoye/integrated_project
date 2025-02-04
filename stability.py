import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import trapz
from wings import fusGeometry as fusgeom
from wings import wingGeometry as winggeom
from sympy import symbols, Eq, solve
from wings import air_density
from wings import true_airspeed_at_altitude

#values to calculate the coefficients 
rho = air_density(12500)[0]
speed = true_airspeed_at_altitude(12500)

#CL for equilibrium 
Cl_tot = 0
CL_T = 1.76

#Important general values 
MAC_tot = 9.949 
surf_tot = 182.56
AR_tot = 1.5
b = 29
MAC_tail = 1.483
hor_tail_surf = 22.5
a = 1.754
a1 = 1.5
a1_over_a = a1/a
weight = 471511.49122
l_fus = 16.8
l_cabin = 10.1
l_cockpit = 2.01
MAC_fus = 13.305
MAC_wing = 4.141

chord_tip_fus = 4.622
sweep_angle_wing = 40*np.pi/180
sweep_angle_fus = 60*np.pi/180

y_AC_wing = 3.796
x_AC_wing = 3.733
y_AC_fus = 2.018


Cm0_fus = -0.1158
Cm0_wing = -0.1533

#Calculation of the Cm0_tot ==> separate the integral in to parts
def Cm0(Cm0_airfoil_fus,Cm0_airfoil_wing):
    _, _, _, _, _, _, _, _, y_wing, leading_wing, trailing_wing, quarter_wing = winggeom() 
    _, _, _, _, _, _, _, y_fus, leading_fus, trailing_fus, quarter_fus = fusgeom() 
    c_fus = trailing_fus - leading_fus
    c_wing = trailing_wing - leading_wing

    Cm0_wing = (2/(surf_tot*MAC_tot)) * trapz(Cm0_airfoil_wing*c_wing**2, y_wing)
    Cm0_fus = (2/(surf_tot*MAC_tot)) * trapz(Cm0_airfoil_fus*c_fus**2, y_fus)
    Cm0_tot = Cm0_wing + Cm0_fus
    return Cm0_tot
Cm0_plane = Cm0(Cm0_fus,Cm0_wing)

print("--------------------------CM0--------------------------------------------")
print("The blended wing body has a Cm0_tot of",Cm0_plane)
print("----------------------------------------------------------------------")

#Position of the important points

x_CG_tot = 9.5
x_AC_tot = 8.789
z_AC_tot = 0

x_AC_tail = 16.5
z_AC_tail = 1.51

#tail volume ratio effectivness 
V_T = hor_tail_surf * (x_AC_tail - x_CG_tot)/(surf_tot* MAC_tot)


config = 2
##################################################################
######CG POSITION
##################################################################

def CG_position(i): 
    wing_weight = 
    wing_pos = (l_fus - chord_tip_fus) + MAC_wing*0.2 + y_AC_wing*np.tan(sweep_angle_wing)

    hor_tail_weight = 
    hor_tail_pos = l_fus - l_cabin - l_cabin
    vert_tail_weight = 
    vert_tail_pos = l_fus - l_cabin - l_cabin

    fus_weight = 
    fus_pos = 0.2*MAC_fus + y_AC_fus*np.tan(sweep_angle_fus)

    land_gear_weight = 
    land_gear_pos = 0.6*MAC_fus + y_AC_fus*np.tan(sweep_angle_fus)

    surf_cont_weight = 
    surf_cont_pos = (l_fus - chord_tip_fus) + MAC_wing*0.4 + y_AC_wing*np.tan(sweep_angle_wing)

    instr_weight = 
    instr_pos = 0.4*l_cockpit

    furn_weight = 
    furn_pos = 0.5*l_fus

    air_cond_weight = 
    air_cond_pos = 

    avio_weight = 
    avio_pos = 
    
    motors_weight = 
    motors_pos = l_cockpit + l_cabin + 1

    elec_syst_weight = 
    elec_syst_pos = 0.75*l_fus/2 + 0.25*motors_pos

    passengers_weight = 8*80
    if i == 1 : 
       passengers_pos = 1 #arbitrary
       passengers_weight = 0
    if i == 2 : #if passengers are as close as possible to the nose
       passengers_pos = l_cockpit
    if i == 3 : #if passengers are as far as possible from the nose 
       passengers_pos = l_cockpit + l_cabin

    total_weight = wing_weight + hor_tail_weight + vert_tail_weight + fus_weight + land_gear_weight + surf_cont_weight + instr_weight + elec_syst_weight + furn_weight + air_cond_weight + avio_weight + passengers_weight + motors_weight
    
    total_mom = wing_weight*wing_pos + hor_tail_weight*hor_tail_pos + vert_tail_weight*vert_tail_pos + fus_weight*fus_pos + land_gear_weight*land_gear_pos + surf_cont_weight*surf_cont_pos + instr_weight*instr_pos + elec_syst_weight*elec_syst_pos + furn_weight*furn_pos + air_cond_weight*air_cond_pos + avio_weight*avio_pos + passengers_weight*passengers_pos + motors_weight*motors_pos
    position = total_mom/total_weight

    return position

print("--------------------------CENTER OF GRAVITY--------------------------------------------")
print("The center of gravity is positionned at",CG_position(config)," meters from the nose of the airplane.")
print("----------------------------------------------------------------------")


##################################################################
######EQUILIBRIUM IN PITCH
##################################################################

#i : configuration numerotation
def CL(i,Cm0_airfoil_fus,Cm0_airfoil_wing): 
    L_tot, L_T = symbols('L_tot L_T')
    Cm_T = V_T * CL_T

    CG_pos = [2,9.5,8]
    if i == 1: 
        x_CG_tot = CG_pos[0]
    if i == 2: 
        x_CG_tot = CG_pos[1]
    if i == 3: 
        x_CG_tot = CG_pos[2]

    Cm0_tot = Cm0(Cm0_airfoil_fus,Cm0_airfoil_wing)
    Cm_tot = Cm0_tot + Cl_tot* (x_CG_tot - x_AC_tot)/MAC_tot - Cm_T
    #translation equilibrium 
    eq1 = Eq(L_tot + L_T - weight,0)
    
    #rotation equilibrium 
    eq2 = Eq(Cm0_tot + L_tot*(x_CG_tot-x_AC_tot)/(1/2*rho*speed**2 * MAC_tot*surf_tot)-L_T*(x_AC_tail - x_CG_tot)/(1/2*rho*speed**2 * MAC_tail*hor_tail_surf)*V_T,0)

    solution = solve((eq1,eq2),(L_tot,L_T))
    L_tot = solution[L_tot]
    L_T = solution[L_T]

    CL_tot = L_tot/(1/2*rho*speed**2 * surf_tot)
    CL_tail = L_T/(1/2*rho*speed**2 * hor_tail_surf)
   
    return L_tot,L_T,CL_tot, CL_tail

print("--------------------------CRUISE--------------------------------------------")
print("The new lift force generated by the body is",CL(config,Cm0_fus,Cm0_wing)[0],"[N]")
print("The new lift force generated by the tail is",CL(config,Cm0_fus,Cm0_wing)[1],"[N]")
print("The new lift coefficient generated by the body is",CL(config,Cm0_fus,Cm0_wing)[2],"[-]")
print("The new lift coefficient generated by the tail is",CL(config,Cm0_fus,Cm0_wing)[3],"[-]")
print("----------------------------------------------------------------------")

def downwash(i,Cm0_airfoil_fus,Cm0_airfoil_wing):
    Cl_tot = CL(i,Cm0_airfoil_fus,Cm0_airfoil_wing)[2]
    deps = 2*a/(np.pi*AR_tot)
    eps = 2*Cl_tot/(np.pi*AR_tot)
    return eps, deps

eps = downwash(config,Cm0_fus,Cm0_wing)

print("--------------------------DOWNWASH--------------------------------------------")
print("The downwash effect eps equals", eps[0]," and deps/dalpha is equal to", eps[1])
print("----------------------------------------------------------------------")

##################################################################
######LONGITUDINAL STATIC STABILITY
##################################################################

def long_stat_stab_cruise(i,Cm0_airfoil_fus,Cm0_airfoil_wing): #in the pitching plane
    #check the stability
    #neutral point : position of the cg in order to have the derivative equals 0
    deps = downwash(i,Cm0_airfoil_fus,Cm0_airfoil_wing)[1]
    
    hn = x_AC_tot/MAC_tot + V_T*a1_over_a*(1- deps) #- 0.5*fus_width**2 * fus_length/(S_wing*a*MAC_wing)#position of the neutral point  

    derivative = hn - x_CG_tot/MAC_tot 
    Kn = - derivative #static margin

    if Kn >= 0.05 and Kn < 0.2 : 
        print("The static margin has a correct value and is equal to : ", (Kn*100), "%")

    else : 
        print("The static margin has to be changed and is equal to : ", (Kn*100), "%")
    
    return 
print("--------------------------STATIC MARGIN--------------------------------------------")
long_stat_stab_cruise(config,Cm0_fus,Cm0_wing)
print("----------------------------------------------------------------------")
