import numpy as np
import math



# -----------------
# Weight Estimation
# -----------------

def get_weight():
    #MTOM = m_e + m_fuel + m_payload
    # Random and not correct values
    MTOM = 117647 #106000  # Maximum Take-Off Mass in lbs 

    s_cabin = 4.622*9 + 9*(10.1+2.113-4.622)/2 # Cabin floor area in m²
    s_cabin_sqft = s_cabin* 10.764 # Cabin floor area in sq.ft

    s_aft = 9 * (16.8 - 10.1 - 2.113)  # Aft center-body area in m²
    s_aft_sqft = s_aft * 10.764  # Aft center-body area in sq.ft
    n_eng = 2  # Number of engines
    lamda_aft = 1  #0.536 ou 1 jsp trop  # Taper ratio of aft centre-body

    c_1 = 0.028  # for long range transport aircraft 
    b = 29  # Wingspan in m
    s = 66.44  # Wing area in m²
    Angle_25 = 32.123  # Quarter-chord sweep angle in degrees
    lamda = 0.161  # Wing taper ratio
    n = 2.5  #  Design normal acceleration factor
    v_cr = 340*0.9  # Cruise speed
    v_D = 1.25*v_cr  # Design dive speed in m/s (EAS : Equivalent Airspeed, airspeed at sea level that would produce the same dynamic pressure as the true airspeed at the aircraft's current altitude)


    tau = 0.086/4.114   #mean thickness of airfoil/mean aerodynamic chord of wing  # Average thickness to chord ratio.
    T_TO = 150 #?????? # Takeoff thrust per engine in kN
    T_pTO = 150 * 220.48089  # Takeoff thrust per engine in pounds of force
    BPR = 5  #?????? # Bypass ratio of engines
    n_fdcrew = 2  # Number of pilots
    n_pax = 8  # Number of passengers
    n_ccrew = 1  # Number of cabin crew
    klav =  3.9 # for long range airplane
    kbuf = 5.68 # for long range airplane
    p_c = 11  #entre 10.9 et 11.8psi d'habitude # Design cabin pressure in psi

    l_cab = 10.1  # Cabin length in meters
    l_cab_ft = l_cab * 3.281 # Cabin length in ft
    m_pax = 215  # Average passenger weight in lbs
    m_bgge = 50  # Baggage allowance per passenger in lbs
    m_container = 100 * 2.205  # Mass of the ULDs holding the baggages in kg
    wine = 12 * 40 #cases of wine [lbs]
    fop = 16  # Operational factor. For long range airline fop = 16 [110].
    fuel_frac = 0.45  #?????  # Fraction of MTOM used for fuel


    # Bradley prediction model for airframe mass :

    # the mass of the centre-body eq. (3.46)  [lbs]
    m_cab = 5.698865 * 0.3166422 * MTOM**(0.166552) * s_cabin_sqft**(1.0161158)   #OK

    #m_cab = 5.698865/450 * n_pax * 0.3166422 * MTOM**(0.166552) * s_cabin_sqft**(1.0161158)

    # The aft centre- body eq. (3.47)  [lbs]
    m_aft = (1 + 0.05 * n_eng) * 0.53 * s_aft_sqft * MTOM**(0.2) * (lamda_aft + 0.5)    #OK


    # Wing (from Howe's but with corrections for use of composites)  [kg]
    m_wing = 0.8 * c_1 * ((b*s/np.cos(np.pi/180*Angle_25)) * ((1 + 2* lamda)/(3 + 3*lamda)) * (MTOM/2.204*n/s)**(0.3) * ((v_D/tau))**(0.5) )**(0.9)   #Le mec l'utilise directement en lbs mais je pense qu'il faut convertir pcq l'article sur lequel c'est basé parle en kg
    "à convertir imo"
    m_wing = m_wing * 2.205

    # Landing Gear [lbs] 
    m_LG = 0.0445 * MTOM



    # Propulsion system [kg]
    m_eng = 2*2000*2.205

    # Nacelle group [kg] 
    m_nacgrp = 0.055 * T_pTO * n_eng




    # APU (auxiliary power unit) [lbs]
    m_APU = 0.001 * MTOM   #OK

    # Mass of engine instruments in [lbs] 
    m_enginst = n_eng * (5 + ((0.006 * MTOM)/1000))
    m_fltinst = 2 * (15 + (0.032 * MTOM)/1000)  # Mass of flight instruments in [lbs]
    m_otherinst = (0.15 * MTOM)/1000 + 0.012 * MTOM # Mass of other instruments in [lbs]

    # Instruments [lbs] 
    m_instr = (m_enginst + m_fltinst + m_otherinst)   #OK

    # Hydraulic system [lbs]
    m_hydr = 3.2 * (MTOM/2.204)**(0.5)  #OK
    m_hydr = m_hydr * 2.204

    n_pax_lux =10*n_pax
    # Furniture [lbs]   # si ça vient du livre mentionné, il parle tt le temps en lbs
    m_furn = ( (55 * n_fdcrew) + (32 * n_pax_lux) + (15 * n_ccrew) + klav * (n_pax_lux**(1.33)) + kbuf * (n_pax_lux**1.12) + 109 * ((n_pax_lux * (1 + p_c) / 100))**0.505 + 0.771 * ( MTOM / 1000) )
    #A predre avec pincettes
    #print(0.412*(n_pax+n_ccrew)**1.145 * MTOM**0.489)



    # AC [lbs]
    m_AC = (6.75 * 3.2808 * l_cab_ft**(1.28))   #OK

    """
    # Payload [kg?] -> maybe don't use the formula and use the mission requirements to determine the payload
    """

    n_pass = 8
    m_pass = 215
    m_bagage = 50 #on prend le worst case des deux dcp
    m_payload = n_pass * ( m_pass + m_bagage) + wine + m_container  #[lbs]
    #m_payload = n_pax*(m_pax + m_bgge) + wine + m_container

    # Operational Items [kg] ici kg, j'ai vérif
    m_ops = 85 * n_fdcrew + (fop + n_pax)
    "à modifier"
    m_ops = m_ops * 2.205

    # Electrical [lbs]
    m_elec = 0.75 * (MTOM/2.204)**(0.67)
    m_elec = m_elec * 2.204

    # Flight control (surface control) [lbs]
    m_fltcon = 0.11 * (MTOM/2.204)**(0.8)
    m_fltcon = m_fltcon * 2.204

    # Fuel [lbs]
    m_fuel = fuel_frac * MTOM


    m_prediction = m_cab + m_aft + m_wing + m_LG + m_eng + m_nacgrp + m_APU + m_instr + m_hydr + m_furn + m_AC + m_payload + m_ops + m_elec + m_fltcon + m_fuel

    # print(f"Mass of the cabin: {m_cab:.2f} lbs")
    # print(f"Mass of the aft centre-body: {m_aft:.2f} lbs")
    # print(f"Mass of the wing: {m_wing:.2f} lbs")
    # print(f"Mass of the landing gear: {m_LG:.2f} lbs")
    # print(f"Mass of the propulsion system: {m_eng:.2f} lbs")
    # print(f"Mass of the nacelle group: {m_nacgrp:.2f} lbs")
    # print(f"Mass of the auxiliary power unit (APU): {m_APU:.2f} lbs")
    # print(f"Mass of the engine instruments: {m_enginst:.2f} lbs")
    # print(f"Mass of the flight instruments: {m_fltinst:.2f} lbs")
    # print(f"Mass of other instruments: {m_otherinst:.2f} lbs")
    # print(f"Total mass of instruments: {m_instr:.2f} lbs")
    # print(f"Mass of the hydraulic system: {m_hydr:.2f} lbs")
    # print(f"Mass of the furniture: {m_furn:.2f} lbs")
    # print(f"Mass of the air conditioning (AC): {m_AC:.2f} lbs")
    # print(f"Mass of the payload: {m_payload:.2f} lbs")
    # print(f"Mass of operational items: {m_ops:.2f} lbs")
    # print(f"Mass of the electrical system: {m_elec:.2f} lbs")
    # print(f"Mass of the flight control system: {m_fltcon:.2f} lbs")
    # print(f"Mass of the fuel: {m_fuel:.2f} lbs")
    # print(f"Total predicted mass: {m_prediction:.2f} lbs")



    components = {
        "Cabin": m_cab,
        "Aft Centre-Body": m_aft,
        "Wing": m_wing,
        "Landing Gear": m_LG,
        "Engines": m_eng,
        "Nacelle Group": m_nacgrp,
        "APU": m_APU,
        "Engine Instruments": m_enginst,
        "Flight Instruments": m_fltinst,
        "Other Instruments": m_otherinst,
        "Total Instruments": m_instr,
        "Hydraulic System": m_hydr,
        "Furniture": m_furn,
        "AC": m_AC,
        "Payload": m_payload,
        "Operational Items": m_ops,
        "Electrical System": m_elec,
        "Flight Control System": m_fltcon,
        "Fuel": m_fuel
    }

    # Printing the mass of each component and its percentage
    for component, mass in components.items():
        percentage = (mass / m_prediction) * 100
        print(f"{component}: {mass:.2f} lbs, {percentage:.2f}%")

    print(f"Total predicted mass: {m_prediction:.2f} lbs")

    return m_cab, m_aft, m_wing,m_LG,m_eng,m_nacgrp,m_APU,m_enginst,m_fltinst,m_fltinst,m_otherinst,m_instr,m_hydr,m_furn,m_AC,m_payload,m_ops,m_elec,m_fltcon,m_fuel

