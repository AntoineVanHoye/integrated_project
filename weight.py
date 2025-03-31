import numpy as np
import math
from wings import detSurfac
from wings import wingGeometry

# Function to calculate air density using the ISA model
# Up to 11 km (Troposphere)
R = 287.05  # [m^2/s^2K]
gamma = 1.4
T_0 = 288.15 #[K]
rho_0 = 1.225  #[kg/m^3]
g = 9.81 #[m/s^2]
L =  0.0065
p_0 = 101325


def air_density(altitude):
    # Up to 11 km (Troposphere)
    if altitude <= 11000:
        T = T_0 - L * altitude  # Temperature [K]
        p = p_0 * (T / 288.15) ** 5.2561  # Pressure [Pa]
    else:
        # Simplification for stratosphere, constant T [K] above 11 km
        T = 216.65  # Constant temperature [K]
        p = 22632 * np.exp(-9.81 * (altitude - 11000) / (287.05 * T))
    rho = p / (287.05 * T)  # Air density [kg/m^3]
    return rho, T, p

def true_airspeed_at_altitude(M, altitude):
    T = air_density(altitude)[1]
    a = np.sqrt(gamma * R * T)  # [m/s] Speed of sound
    v = M * a  # [m/s] Aircraft velocity
    return v

# -----------------
# Weight Estimation
# -----------------

def get_weight():

    # Initial guess for MTOW
    MTOW = 104799  # Maximum Take-Off Mass in lbs (first converge value : 117647)
    tolerance = 0.1  # Convergence criterion
    max_iterations = 100  # Safety limit
    iteration = 0

    while iteration < max_iterations:
        iteration += 1

        # --- Constants ---
        s_cabin = 4.622*9 + 9*(10.1+2.113-4.622)/2 # Cabin floor area in m²
        s_cabin_sqft = s_cabin* 10.764 # Cabin floor area in sq.ft
    
        s_aft = 9 * (16.8 - 10.1 - 2.113)  # Aft center-body area in m²
        s_aft_sqft = s_aft * 10.764  # Aft center-body area in sq.ft
        n_eng = 2  # Number of engines
        lamda_aft = 1  #0.536 ou 1 jsp trop  # Taper ratio of aft centre-body
    
        c_1 = 0.028  # for long range transport aircraft 
        b = 28.95  # Wingspan in m
        Cl = 0.45
        sweep_LE_fus = 50
        sweep_quarter_wing = 29
        force = 557174.8525469044
        surface_wing_ideal, surf_fus, s, surf_tot = detSurfac(Cl, sweep_LE_fus, sweep_quarter_wing, force)
        Angle_25 = 29 # Quarter-chord sweep angle in degrees #change
        _, AR_wing, sweep_beta, sweep_beta_tot, c_root, lamda, sweep_quarter, c_tip, y, leading_edge, trailing_edge, quarter_line, c, h = wingGeometry(Cl,sweep_LE_fus, sweep_quarter_wing, force)
        n = 2.5  #  Design normal acceleration factor

        altitude = 12500
        rho_alt, T_alt, p_alt = air_density(altitude)
        v_cr = true_airspeed_at_altitude(0.9, altitude)
        v_cr = v_cr * np.sqrt(rho_alt/rho_0)
        v_D = 1.25 * v_cr  # Design dive speed in m/s (EAS : Equivalent Airspeed, airspeed at sea level that would produce the same dynamic pressure as the true airspeed at the aircraft's current altitude)
    
        tau = 0.077/3.476   #mean thickness of airfoil/mean aerodynamic chord of wing  # Average thickness to chord ratio. #change
        T_TO = 150 #?????? # Takeoff thrust per engine in kN
        T_pTO = 95*2* 220.48089  # Takeoff thrust per engine in pounds of force
        BPR = 5  #?????? # Bypass ratio of engines
        n_fdcrew = 2  # Number of pilots
        n_pax = 8  # Number of passengers
        n_ccrew = 1  # Number of cabin crew
        klav =  3.9 # for long range airplane
        kbuf = 5.68 # for long range airplane
        p_c = 11  #entre 10.9 et 11.8psi d'habitude # Design cabin pressure in psi
    
        l_cab = 10.1  # Cabin length in meters
        l_cab_ft = l_cab * 3.281 # Cabin length in ft
        #m_pax = 215  # Average passenger weight in lbs
        #m_bgge = 50  # Baggage allowance per passenger in lbs
        m_container = 100 * 2.205  # Mass of the ULDs holding the baggages in kg
        wine = 12 * 40 #cases of wine [lbs]
        fop = 16  # Operational factor. For long range airline fop = 16 [110].
        fuel_frac = 0.45  #?????  # Fraction of MTOM used for fuel


        # --- Weight Calculations ---
        # Bradley prediction model for airframe mass :
    
        # the mass of the centre-body eq. (3.46)  [lbs]
        m_cab = 5.698865 * 0.3166422 * MTOW**(0.166552) * s_cabin_sqft**(1.0161158)   #OK
    
        #m_cab = 5.698865/450 * n_pax * 0.3166422 * MTOM**(0.166552) * s_cabin_sqft**(1.0161158)
    
        # The aft centre- body eq. (3.47)  [lbs]
        m_aft = (1 + 0.05 * n_eng) * 0.53 * s_aft_sqft * MTOW**(0.2) * (lamda_aft + 0.5)    #OK
    
        # Wing (from Howe's but with corrections for use of composites)  [kg]
        m_wing = 0.8 * c_1 * ((b*s/np.cos(np.pi/180*Angle_25)) * ((1 + 2* lamda)/(3 + 3*lamda)) * (MTOW/2.204*n/s)**(0.3) * ((v_D/tau))**(0.5) )**(0.9)   #Le mec l'utilise directement en lbs mais je pense qu'il faut convertir pcq l'article sur lequel c'est basé parle en kg
        m_wing = m_wing * 2.205 # converted in lbs
    
        # Landing Gear [lbs] 
        m_LG = 0.0445 * MTOW
    
        # Propulsion system [lbs]
        m_eng = 2*2000*2.205
    
        # Nacelle group [lbs] 
        m_nacgrp = 0.055 * T_pTO * n_eng
    
        # APU (auxiliary power unit) [lbs]
        m_APU = 0.001 * MTOW   #OK
    
        # Mass of engine instruments in [lbs] 
        m_enginst = n_eng * (5 + ((0.006 * MTOW)/1000))
        m_fltinst = 2 * (15 + (0.032 * MTOW)/1000)  # Mass of flight instruments in [lbs]
        m_otherinst = (0.15 * MTOW)/1000 + 0.012 * MTOW # Mass of other instruments in [lbs]
    
        # Instruments [lbs] 
        m_instr = (m_enginst + m_fltinst + m_otherinst)   #OK
    
        # Hydraulic system [lbs]
        m_hydr = 3.2 * (MTOW/2.204)**(0.5)  #OK
        m_hydr = m_hydr * 2.204
    
        n_pax_lux =10*n_pax
        # Furniture [lbs]   # si ça vient du livre mentionné, il parle tt le temps en lbs
        m_furn = ( (55 * n_fdcrew) + (32 * n_pax_lux) + (15 * n_ccrew) + klav * (n_pax_lux**(1.33)) + kbuf * (n_pax_lux**1.12) + 109 * ((n_pax_lux * (1 + p_c) / 100))**0.505 + 0.771 * ( MTOW / 1000) )
        #A predre avec pincettes
        #print(0.412*(n_pax+n_ccrew)**1.145 * MTOM**0.489)
    
        # AC [lbs]
        m_AC = (6.75 * 3.2808 * l_cab_ft**(1.28))   #OK
    
        """
        # Payload [kg?] -> formula with the mission requirements
        """
        n_pass = 8
        m_pass = 215 # per passenger (value from the requirements)
        m_bagage = 50 # on prend le worst case des deux dcp
        #m_payload_1 = n_pass * ( m_pass + m_bagage) + wine + m_container  #[lbs]
        m_payload = n_pass * m_bagage + wine + m_container  #[lbs] (without taking account of the passenger)
        m_passenger = n_pass * m_pass # (Normally, it should be taken into account in the payload, but we need a separated value for stability)
    
        # Operational Items [kg] ici kg, j'ai vérif
        m_ops_kg = 85 * n_fdcrew + (fop + n_pax)
        m_ops = m_ops_kg * 2.205 # in [lbs]
    
        # Electrical [lbs]
        m_elec = 0.75 * (MTOW/2.204)**(0.67)
        m_elec = m_elec * 2.204
    
        # Flight control (surface control) [lbs]
        m_fltcon = 0.11 * (MTOW/2.204)**(0.8)
        m_fltcon = m_fltcon * 2.204
    
        # Fuel [lbs]
        m_fuel = fuel_frac * MTOW
    
        # New values corrected 
        m_eng = 4354*2.20462 # Rolls Royce Pearl 700 (no data for the 10X)   #  8377.566 (old engine value) 
        m_fuel = 26102.23125*0.8*2.20462 # Replace the value computed w/ a % of the MOTW by the one calculated in the propulsion part.
        
        # --- Compute new MTOW estimate --- 
        m_prediction = m_cab + m_aft + m_wing + m_LG + m_eng + m_nacgrp + m_APU + m_instr + m_hydr + m_furn + m_AC + m_payload + m_passenger + m_ops + m_elec + m_fltcon + m_fuel
        m_prediction_kg =  m_prediction / 2.20462
        # Check convergence
        if abs(m_prediction - MTOW) < tolerance:
            break
        MTOW = m_prediction
        
    # --- Print component weights and percentages ---
    components = {
        "Cabin": m_cab, "Aft Centre-Body": m_aft, "Wing": m_wing, "Landing Gear": m_LG,
        "Engines": m_eng, "Nacelle Group": m_nacgrp, "APU": m_APU, "Instruments": m_instr, "Hydraulic System": m_hydr, "Furniture": m_furn,
        "AC": m_AC, "Payload": m_payload, "Operational Items": m_ops,
        "Electrical System": m_elec, "Flight Control System": m_fltcon, "Fuel": m_fuel
    }

    for component, mass in components.items():
        percentage = (mass / m_prediction) * 100
        #print(f"{component}: {mass:.2f} lbs, {percentage:.2f}%")
    
    print(f"Total predicted mass (MOTW): {m_prediction:.2f} lbs")
    #print(f"Total predicted mass (MOTM): {m_prediction_kg:.2f} kg")
    #print("") # Vertical line space
    #print(f"The MOTW converged in {iteration} iterations.")

    return (m_cab, m_aft, m_wing, m_LG, m_eng, m_nacgrp, m_APU, m_enginst, m_instr, m_hydr,
            m_furn, m_AC, m_payload, m_ops, m_elec, m_fltcon, m_fuel, m_prediction, m_passenger)

print(get_weight())
