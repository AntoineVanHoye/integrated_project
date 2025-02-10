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
from wings import get_Lift_and_drag
from wings import getAR
from wings import getSweep
from wings import getAirfoilFus
from wings import getAirfoilWing
from wings import getClAlfa

#values to calculate the coefficients 
rho = air_density(12500)[0]
speed = true_airspeed_at_altitude(12500)

delta = 0.005
AR_tot = getAR()
#Important general values 
MAC_fus, y_AC_fus,x_AC_fus,MAC_wing,y_AC_wing,x_AC_wing,MAC_tot,y_AC_tot,x_AC_tot = getMAC()


surf_tot,surf_fus,surf_wing = detSurfac()
_,_,_,_,_,_,_,_,_= get_Lift_and_drag(AR_tot,delta)
a = getClAlfa()
b = 29
#a = 8
mean_chord = surf_tot/b

hor_tail_surf = 47#34.75
a1 = 3.16
a1_over_a = a1/a
l_fus = 16.8
l_cabin = 10.1
l_cockpit = 2.01
l_aft = l_fus - l_cabin - l_cockpit

_,_,_,_,_,_,chord_tip_fus,_,_,_,_ = fusGeometry()
sweep_angle_wing = getSweep()[1]*np.pi/180
sweep_angle_fus = getSweep()[0]*np.pi/180

Cm0_wing = getAirfoilWing()[4]
Cm0_fus = getAirfoilFus()[4]

#Calculation of the Cm0_tot ==> separate the integral in to parts
def Cm0(Cm0_airfoil_fus,Cm0_airfoil_wing):
    _, _, _, _, _, _, _, _, y_wing, leading_wing, trailing_wing, quarter_wing,_,_ = winggeom() 
    _, _, _, _, _, _, _, y_fus, leading_fus, trailing_fus, quarter_fus = fusgeom() 
    c_fus = trailing_fus - leading_fus
    c_wing = trailing_wing - leading_wing
    Cm0_wing = (2/(surf_wing*MAC_wing)) * trapz(Cm0_airfoil_wing*c_wing**2, y_wing)
    #Cm0_wing = 2*trapz(Cm0_airfoil_wing*c_wing**2, y_wing)
    Cm0_fus = (2/(surf_fus*MAC_fus)) * trapz(Cm0_airfoil_fus*c_fus**2, y_fus)
    #Cm0_fus = 2 * trapz(Cm0_airfoil_fus*c_fus**2, y_fus)
    #Cm0_tot = (Cm0_wing + Cm0_fus)/(MAC_tot*surf_tot)
    Cm0_tot = (Cm0_wing*surf_wing + Cm0_fus*surf_fus)/surf_tot
    return Cm0_tot,Cm0_fus,Cm0_wing


print("--------------------------CM0--------------------------------------------")
print("The blended wing body has a Cm0_tot of",Cm0(Cm0_fus,Cm0_wing)[0])
print("The wings have a pitching moment (up) of",Cm0(Cm0_fus,Cm0_wing)[2],". It represents",Cm0(Cm0_fus,Cm0_wing)[2]/Cm0(Cm0_fus,Cm0_wing)[0]*100,"% of the pitching moment of the airplane.")
print("The fuselage has a pitching moment (up) of",Cm0(Cm0_fus,Cm0_wing)[1],". It represents",Cm0(Cm0_fus,Cm0_wing)[1]/Cm0(Cm0_fus,Cm0_wing)[0]*100,"% of the pitching moment of the airplane.")
print("----------------------------------------------------------------------")

#Position of the important points
z_AC_tot = 0
z_CG_tot = 0
MAC_tail = 4.56
z_CG_motors = 1

x_AC_tail_local = 0.91
x_AC_tail = l_cabin + l_cockpit + x_AC_tail_local
l_tail = 4.26


z_AC_tail = 1.51

config = 1
fuel = 2
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

def CG_position(i,d): 

    fus_weight, aft_weight, wing_weight, land_gear_weight,motors_weight,nacelle_weight,APU_weight,enginst_weight,instr_weight,hydr_syst_weight,furn_weight,air_cond_weight,payload_weight, ops_weight,elec_syst_weight,surf_cont_weight,_,_,_= get_weight()

    wing_pos = (l_fus - chord_tip_fus) + MAC_wing*0.2 + y_AC_wing*np.tan(sweep_angle_wing)
    
    fus_pos = (l_cabin + l_cockpit)*0.6 #0.2*MAC_fus + y_AC_fus*np.tan(sweep_angle_fus)

    aft_pos = l_cockpit + l_cabin + l_aft*0.4

    APU_pos = l_cockpit + l_cabin

    hydr_pos = 0.75*((l_fus - chord_tip_fus) + MAC_wing*0.2 + y_AC_wing*np.tan(sweep_angle_wing)) + 0.25 *(l_cockpit + l_cabin + 2)

    payload_pos = 0.75*(l_cockpit + l_cabin) + 0.25*l_cabin

    ops_pos = l_cockpit + 1

    land_gear_pos = (0.6*MAC_fus + y_AC_fus*np.tan(sweep_angle_fus))*0.8 + 0.2*0.22*l_fus

    surf_cont_pos = (l_fus - chord_tip_fus) + MAC_wing*0.4 + y_AC_wing*np.tan(sweep_angle_wing)

    instr_pos = 0.4*l_cockpit

    furn_pos = 0.5*l_fus

    air_cond_pos = l_cockpit + l_cabin + 0.5
    
    motors_pos = l_cockpit + l_cabin + 1

    nacelle_pos = motors_pos

    enginst_pos = motors_pos

    elec_syst_pos = 0.75*l_fus/2 + 0.25*motors_pos

    available_fuel_vol = wingFuelvolume()*1000
    if d == 1 : #no fuel
        fuel_weight = 0
        fuel_pos = 1
        pourc_wings = 0
    if d == 2 : #full of fuel
        vol_fuel = 26854.56
        fuel_weight = vol_fuel*0.8*2.20462
        pourc_wings = available_fuel_vol/vol_fuel
        fuel_pos = wing_pos*pourc_wings + (1-pourc_wings)*l_fus*0.44

    passengers_weight = passengers(i)[0]
    passengers_pos = passengers(i)[1] 

    """
    print("fuel:",fuel_pos,"->",fuel_pos/l_fus*100)
    print("wing:",wing_pos, "->", wing_pos/l_fus *100)
    print("fus:",fus_pos, "->", fus_pos/l_fus *100)
    print("aft:",aft_pos, "->", aft_pos/l_fus *100)
    print("APU:",APU_pos, "->", APU_pos/l_fus *100)
    print("hydr:",hydr_pos, "->", hydr_pos/l_fus *100)
    print("payload:",payload_pos, "->", payload_pos/l_fus *100) 
    print("ops:",ops_pos, "->", ops_pos/l_fus *100)
    print("land_gear:",land_gear_pos, "->", land_gear_pos/l_fus *100)
    print("surf_cont:",surf_cont_pos, "->", surf_cont_pos/l_fus *100)
    print("instr:",instr_pos, "->", instr_pos/l_fus *100)
    print("furn:",furn_pos, "->", furn_pos/l_fus *100)
    print("air_cond:",air_cond_pos, "->", air_cond_pos/l_fus *100)
    print("motors:",motors_pos, "->", motors_pos/l_fus *100)
    print("nacelle:",nacelle_pos, "->", nacelle_pos/l_fus *100)
    print("enginst:",enginst_pos, "->", enginst_pos/l_fus *100)
    print("elec_syst:",elec_syst_pos, "->", elec_syst_pos/l_fus *100)
    """
    #total_weight = wing_weight + fus_weight + land_gear_weight + surf_cont_weight + instr_weight + elec_syst_weight + furn_weight + air_cond_weight + passengers_weight + motors_weight + fuel_weight + aft_weight + nacelle_weight + APU_weight + enginst_weight + hydr_syst + payload_weight + ops_weight 
    total_mom = (wing_weight*wing_pos) + (fus_weight*fus_pos) + (land_gear_weight*land_gear_pos) + (surf_cont_weight*surf_cont_pos) + (instr_weight*instr_pos) + (elec_syst_weight*elec_syst_pos) + (furn_weight*furn_pos) + (air_cond_weight*air_cond_pos) + (passengers_weight*passengers_pos) + (motors_weight*motors_pos) + (fuel_pos*fuel_weight) + (aft_pos*aft_weight) + (nacelle_pos*nacelle_weight) + (APU_pos*APU_weight) + (enginst_pos*enginst_weight) + (hydr_pos*hydr_syst_weight) + (payload_pos*payload_weight) + (ops_pos*ops_weight) 
    total_weight = wing_weight + fus_weight + land_gear_weight + surf_cont_weight + instr_weight + elec_syst_weight + furn_weight + air_cond_weight + passengers_weight + motors_weight + fuel_weight + aft_weight + nacelle_weight + APU_weight + enginst_weight + hydr_syst_weight + payload_weight + ops_weight
    position = total_mom/total_weight

    return position, pourc_wings, motors_pos/MAC_tot, total_weight 

print("--------------------------TAIL--------------------------------------------")
print("The position of the aerodynamic center of the tail is at",x_AC_tail," m and the end of the tail is at",x_AC_tail- x_AC_tail_local + l_tail,"m from the nose and the distance between x_AC_tail and x_CG_tot is",x_AC_tail - CG_position(config,fuel)[0],"m.")
print("----------------------------------------------------------------------")

print("--------------------------FUEL STORAGE--------------------------------------------")
print(CG_position(config,fuel)[1]*100,"% of the total fuel volume is stored in the wings and we need",(1-CG_position(config,fuel)[1])*26854.56/1000,"m³ in the fuselage.")
print("----------------------------------------------------------------------")

print("--------------------------CENTER OF GRAVITY--------------------------------------------")
print("The center of gravity is positioned at",CG_position(config,fuel)[0],"m (",CG_position(config,fuel)[0]*100/l_fus,"%) from the nose of the airplane.")
print("----------------------------------------------------------------------")

#tail volume ratio effectivness 
def tail_eff(i,d):
    x_CG_tot = CG_position(i,d)[0]
    V_T = hor_tail_surf * (x_AC_tail - x_CG_tot)/(surf_tot* MAC_tot)
    return V_T

def prop_force(): 
    A = 2*1.32**2 *np.pi
    m_dot = rho*A*speed*0.453592/32.2
    Fp = m_dot*speed
    return Fp

def wings_lift():
    _,CL_w0,_,_,_,_ = wingCL()
    L_w = CL_w0 * (1/2)*rho*surf_wing*speed**2
    return L_w

def fus_lift(): 
    _,CL_f0,_,_,_ = fuselageCL()
    L_f = CL_f0*(1/2)*rho*surf_fus*speed**2
    return L_f

##################################################################
######EQUILIBRIUM IN PITCH
##################################################################

#i : configuration numerotation
def CL(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing): 
    L_tot, L_T = symbols('L_tot L_T')
    V_T = tail_eff(i,d)
    T = 42676.35
    #T = 160000
    #drag_tail = 
    x_CG_tot = CG_position(i,d)[0]
    weight = CG_position(i,d)[3]*9.81*0.453592 + passengers(i)[0]*9.81*0.453592
    Cm0_tot = Cm0(Cm0_airfoil_fus,Cm0_airfoil_wing)[0]
    M0 = Cm0_tot*((1/2)*rho*(speed**2)*MAC_tot*surf_tot)
    x_CG_motors = CG_position(i,d)[2]
    Fp = prop_force()

    #Cm_tot = Cm0_tot + Cl_tot* (x_CG_tot - x_AC_tot)/MAC_tot - Cm_T
    #translation equilibrium 
    eq1 = Eq(L_tot + L_T - weight + Fp,0)
    
    #rotation equilibrium 
    #eq2 = Eq(Cm0_tot + L_tot*(x_CG_tot-x_AC_tot)/(1/2*rho*speed**2 * MAC_tot*surf_tot)-L_T*(x_AC_tail - x_CG_tot)/(1/2*rho*speed**2 * MAC_tail*hor_tail_surf)*V_T+T*(z_CG_motors - z_CG_tot),0)
    eq2 = Eq(M0 + L_tot*(x_CG_tot-x_AC_tot)-L_T*(x_AC_tail - x_CG_tot)-T*(z_CG_motors - z_CG_tot) + Fp*(x_CG_tot-x_CG_motors),0) #+ drag_tail*(z_AC_tail-z_CG_tot)

    solution = solve((eq1,eq2),(L_tot,L_T))
    L_tot = solution[L_tot]
    L_T = solution[L_T]

    return L_tot,L_T

print("--------------------------CRUISE--------------------------------------------")
print("The new lift force generated by the body is",CL(config,fuel,Cm0_fus,Cm0_wing)[0],"[N]")
print("The new lift force generated by the tail is",CL(config,fuel,Cm0_fus,Cm0_wing)[1],"[N]")
print("----------------------------------------------------------------------")

def downwash(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing):
    #Cl_tot = CL(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing)[2]
    deps = 2*a/(np.pi*AR_tot)
    #eps = 2*Cl_tot/(np.pi*AR_tot)
    eps = 0
    deps = 0
    return eps, deps

print("--------------------------DOWNWASH--------------------------------------------")
print("The downwash effect eps equals", downwash(config,fuel,Cm0_fus,Cm0_wing)[0]," and deps/dalpha is equal to", downwash(config,fuel,Cm0_fus,Cm0_wing)[1])
print("----------------------------------------------------------------------")

##################################################################
######LONGITUDINAL STATIC STABILITY
##################################################################

def long_stat_stab_cruise(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing): #in the pitching plane
    #check the stability
    #neutral point : position of the cg in order to have the derivative equals 0
    x_CG_tot = CG_position(i,d)[0]
    engines_pos = CG_position(i,d)[2]
    deps = downwash(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing)[1]
    V_T = tail_eff(i,d)
    dalpha_prop = 1 - deps
    Fp = prop_force()
    eta = 0.9
    #hn = x_AC_tot/MAC_tot + V_T*a1_over_a*(1- deps) + dalpha_prop * Fp * (engines_pos-x_CG_tot)/MAC_tot#- 0.5*fus_width**2 * fus_length/(S_wing*a*MAC_wing)  #position of the neutral point  

    #page 416 Raymer
    q = 1/2*rho*speed**2
    hn = (a * x_AC_tot/MAC_tot +eta*a1*hor_tail_surf/surf_tot*(1-deps)*x_AC_tail/MAC_tot+Fp/(q*surf_tot)*(1-deps)*engines_pos)/(a + eta*hor_tail_surf/surf_tot*a1*(1-deps)+Fp/(q*surf_tot))

    Kn = hn - x_CG_tot/MAC_tot 

    if Kn >= 0.05 and Kn < 0.3 : 
        print("The static margin has a correct value and is equal to : ", (Kn*100), "% and the neutral point is positioned at",hn*MAC_tot,"from the nose, which represents",hn*MAC_tot*100/l_fus,"% of the total length.")

    else : 
        print("The static margin has to be changed and is equal to : ", (Kn*100), "% and the neutral point is positioned at",hn*MAC_tot,"m from the nose, which represents",hn*MAC_tot*100/l_fus,"% of the total length.")
    
    return 
print("--------------------------STATIC MARGIN AND NEUTRAL POINT--------------------------------------------")
long_stat_stab_cruise(config,fuel,Cm0_fus,Cm0_wing)
print("----------------------------------------------------------------------")

##################################################################
######CG POSITIONS RANGE
##################################################################

def get_CG(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing,Kn): 

    engines_pos = CG_position(i,d)[2]
    deps = downwash(i,d,Cm0_airfoil_fus,Cm0_airfoil_wing)[1]
    V_T = tail_eff(i,d)
    dalpha_prop = 1 - deps
    Fp = prop_force()
    eta = 0.9
    q = 1/2*rho*speed**2
    hn = (a * x_AC_tot/MAC_tot + eta*a1*hor_tail_surf/surf_tot*(1-deps)*x_AC_tail/MAC_tot+Fp/(q*surf_tot)*(1-deps)*engines_pos)/(a + eta*hor_tail_surf/surf_tot*a1*(1-deps)+Fp/(q*surf_tot))
    #x_CG = x_AC_tot/(1 + hor_tail_surf/surf_tot*a1_over_a*(1 - deps)) + x_AC_tail*a1_over_a*(1-deps)*hor_tail_surf/((1 + hor_tail_surf/surf_tot*a1_over_a*(1 - deps))*surf_tot)- Kn/(1 + hor_tail_surf/surf_tot*a1_over_a*(1 - deps))*MAC_tot
    x_CG = (hn*MAC_tot) - (Kn*MAC_tot)
    if Kn == 0.05:
        print("The center of gravity is positioned at",x_CG,"m from the nose when the static margin is equal to",Kn*100,"%. It is the maximal value of the range.")
    if Kn == 0.15:
        print("The center of gravity is positioned at",x_CG,"m from the nose when the static margin is equal to",Kn*100,"%. It is the minimal value of the range.")
    return 

print("--------------------------ACCEPTABLE POSITIONS FOR THE CENTER OF GRAVITY--------------------------------------------")
get_CG(config,fuel,Cm0_fus,Cm0_wing,0.05)
get_CG(config,fuel,Cm0_fus,Cm0_wing,0.15)
print("----------------------------------------------------------------------")