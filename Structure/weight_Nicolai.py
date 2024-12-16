import math

# ----------------------------
# Section 20.2.1.1 : Structure
# ----------------------------

# Credit for the following weight-estimation methods for conventional 
# metal aircraft goes to many sources in the aerospace industry. Th e weight 
# equations give the component weight in pounds

def wing_weight(S_w, AR, M_0, W_TO, N, t_c, Lambda_half_chord, taper) :
    
    # This wing weight equation is valid for an M0 range of 0.4–0.8, a t/c range of 0.08–0.15, and an aspect ratio (AR) range of 4–12.
    
    # Convert Lambda_half_chord to radians for the cosine function
   Lambda_half_chord_rad = math.radians(Lambda_half_chord)

   # Calculate the wing weight using the provided formula
   Wt = (0.00428 * (S_w **0.48) * (AR **1.0) * (M_0 **0.43) *
         ((W_TO * N) **0.84) * (taper **0.14) * ((100 * t_c) **-0.76) *
         (math.cos(Lambda_half_chord_rad) **-1.54))

   return Wt

# Parameters:
S_w = 1240.0 # Wing area in square feet (ft^2).
AR = 10.0 # Wing aspect ratio.
M_0 = 0.925 # Maximum Mach number at sea level.
W_TO = 104000.0 # Takeoff weight in pounds (lb).
N = 3.75 # Ultimate load factor.
t_c = 0.12 # Maximum thickness ratio (dimensionless).
Lambda_half_chord = 35.0 # Sweep of half-chord in degrees.
taper = 0.3 # Taper ratio


def horizontal_tail_weight(W_TO, N, S_HT, b_HT, t_HT, c_wing, L_t):
    """
    Calculate the weight of the horizontal tail.

    Parameters:
        W_TO (float): Takeoff weight in pounds (lb).
        N (float): Ultimate load factor.
        S_HT (float): Horizontal tail total planform area, including fuselage carry-through, in square feet.
        b_HT (float): Span of horizontal tail in feet.
        t_HT (float): Thickness of the horizontal tail at the root, in feet.
        c_wing (float): Mean aerodynamic chord of the wing, in feet.
        L_t (float): Tail moment arm in feet.

    Returns:
        float: Horizontal tail weight.
    """
    # Calculate gamma for horizontal tail
    gamma = (W_TO * N) ** 0.813 * S_HT ** 0.584 * ((b_HT / t_HT) ** 0.033) * (c_wing / L_t) ** 0.28
    
    # Calculate horizontal tail weight
    Wt = 0.0034 * (gamma ** 0.915) 
    return Wt

S_HT = 300.0     # Horizontal tail planform area in square feet (include fuselage carry-through)
b_HT = 38      # Span of horizontal tail in feet
t_HT = 2.2       # Thickness of the horizontal tail at the root, in feet
c_wing = 12.0    # Mac of the wing, in feet
L_t = 50       # Tail moment arm, in feet; distance from one-fourth wing mac to one fourth tail mac.

def vertical_tail_weight(W_TO, N, S_VT, M_0, L_t, AR_VT, lambda_vt, S_r, Lambda_VT, h_rh_lvt_ratio):
    """
    Calculate the weight of the vertical tail.

    Parameters:
        W_TO (float): Takeoff weight in pounds (lb).
        N (float): Ultimate load factor.
        S_VT (float): Vertical tail area in square feet.
        M_0 (float): Maximum Mach number at sea level.
        L_t (float): Tail moment arm in feet.
        AR_VT (float): Aspect ratio of vertical tail.
        lambda_vt (float): Vertical tail taper ratio.
        S_r (float): Rudder area in square feet (use S_r/S_VT = 0.3 if unknown).
        Lambda_VT (float): Sweep of vertical tail quarter-chord in degrees.
        h_rh_lvt_ratio (float): Ratio of horizontal tail height to vertical tail height.

    Returns:
        float: Vertical tail weight.
    """
    # Convert Lambda_VT to radians for cosine
    Lambda_VT_rad = math.radians(Lambda_VT)

    # Calculate gamma for vertical tail
    gamma = ((1 + h_rh_lvt_ratio) ** 0.5) * ((W_TO * N) ** 0.363) * (S_VT ** 1.089) * (M_0 ** 0.601) * (L_t ** -0.726) * ((1 + S_r / S_VT) ** 0.217) * ((AR_VT ** 0.337) * ((1 + lambda_vt) ** 0.363)) * (math.cos(Lambda_VT_rad) ** -0.484)

    # Calculate vertical tail weight
    Wt = 0.19 * (gamma ** 1.014)
    return Wt

h_rh_lvt_ratio = 1.0  # Ratio of horizontal tail height to vertical tail height. For a “T” tail this ratio is 1.0; for a fuselage-mounted horizontal tail this ratio is 0
S_VT = 220.0          # Vertical tail area in square feet
AR_VT = 1.7           # Aspect ratio of vertical tail
lambda_vt = 0.4      # Taper ratio of vertical tail
S_r = 50.0            # Rudder area in square feet
Lambda_VT = 35     # Sweep of vertical tail in degrees


def fuselage_weight_usaf(W_TO, K_INL, q, L, H):
    """
    Calculate fuselage weight using USAF and Commercial relation.

    Parameters:
        W_TO (float): Takeoff weight in pounds (lb).
        K_INL (float): Inlet factor (1.25 for fuselage inlets, 1.0 elsewhere).
        q (float): Maximum dynamic pressure in lb/ft^2.
        L (float): Fuselage length in feet.
        H (float): Maximum fuselage height in feet.

    Returns:
        float: Fuselage weight.
    """
    Wt = 10.43 * (K_INL ** 1.42) * ((q * 1e-2) ** 0.283) * ((W_TO * 1e-3) ** 0.95) * ((L / H) ** 0.71)
    return Wt

K_INL = 1.25  # Inlet factor; = 1.25 for inlets on fuselage; = 1.0 for inlets in wing root or elsewhere 
q = 600.0     # Dynamic pressure in lb/ft^2
L = 111     # Fuselage length in feet
H = 13      # Fuselage height in feet

def landing_gear_weight_usaf(W_TO):
    """
    Calculate the landing gear weight using USAF and Commercial relation.

    Parameters:
        W_TO (float): Takeoff weight in pounds (lb).

    Returns:
        float: Landing gear weight.
    """
    Wt = 62.21 * ((W_TO * 1e-3) ** 0.84)
    return Wt

# -----------------------------
# Section 20.2.1.2 : Propulsion
# -----------------------------

# Engine : 
    
# Engine weights should be based upon the engine manufacturer’s data 
#  and scaling factors. Assume the exhaust, cooling, turbo-supercharger, and 
#  lubrication systems weights are included in the engine weight.

# Propulsion Subsystems : Propulsion subsystem items are the air induction system, fuel system, 
#  engine controls, and starting system

def prop_sub_syst (N_i, L_d, A_i, P_2, F_gw, F_gf, K_eco, L_f, N_e, W_eng) : 
    
    # Air induction system : (uncomment the one you want to use for W_t_1)
    "JET ENGINE"
    # For  External Turbojet Cowl and Duct Weight:
    # W_t_1 = 3.00 * N_i * (A_i**0.5 * L_d * P_2)**0.731
        # NB :  This equation accounts for the exterior cowl or cover panels, ducting, 
        # and substructure such as rings, frames, stiff eners, and longerons, from the 
        # inlet lip to the engine compressor face, and should only be used for external 
        # engine installations.
    "TURBO FAN"
    # For External Turbofan Cowl and Duct Weight:
    W_t_1 = 7.435 * N_i * (L_d * A_i**0.5 * P_2)**0.731
         
    # Fuel system : 
    
    # Non self-sealing bladder cells : 
    W_t_2 = 23.10 * ((F_gw + F_gf)*10**-2)**0.758
        
    # Fuel System Bladder Cell Backing and Supports :
    W_t_3 = 7.91 * ((F_gw + F_gf)*10**-2)**0.854
        
    # Dump-and-Drain System : 
    W_t_4 = 7.38 * ((F_gw + F_gf)*10**-2)**0.458
        
    #  C.G. Control System (Transfer Pumps and Monitor):
    W_t_5 = 28.38 * ((F_gw + F_gf)*10**-2)**0.442
        
    # Engine Controls : 
    
    # Body- or Wing-Root-Mounted Jet: (in our case they are body mounted) 
    W_t_6 = K_eco * (L_f * N_e)**0.792
    
    # Starting systems :
    "JET ENGINE"   
    # One or Two Jet Engines—Cartridge and Pneumatic :
    # W_t_7_jet = 9.33 * (N_e * W_eng * 10**-3)**1.078    
    # One or Two Jet Engines— Electrical :
    # W_t_7 = W_t_7_jet + 38.98 * (N_e * W_eng * 10**-3)**0.918   
     
    "TURBO PROP"
    # Turboprop Engines— Pneumatic:
    W_t_7 = 12.05*(N_e * W_eng * 10**-3)**1.458
        
    W_t = W_t_1 + W_t_2 + W_t_3 + W_t_4 + W_t_5 + W_t_6 + W_t_7
    W_fs = W_t_2 + W_t_3 + W_t_4 + W_t_5 # Fuel system needed for anther formula

    return W_t, W_fs

N_i = 2 # Number of inlets, vehicle configuration
L_d = 8 # Subsonic duct length, per inlet, in feet (ft) 
A_i = 8 # Capture area per inlet, in square feet (ft2)
P_2 = 35 # Maximum static pressure at engine compressor face, in pounds per square inch absolute (psia)
F_gw = 4000 # Total wing fuel in gallons
F_gf = 2000 # Total fuselage fuel in gallons.
K_eco = 0.686 #  KECO = engine control engine-type coefficient ; = 0.686 for nonafterburning engines; = 1.080 for afterburning (A/B) engines
L_f = L # fuselage length, in feet (ft)
N_e = 2   # Number of engines 
W_eng = 3000 # Engine weight, in pounds per engine.

# ------------------------------------------------------------------
# Section 20.2.1.3 : Surface Controls Plus Hydraulics and Pneumatics
# ------------------------------------------------------------------

def surf_contrl_hydr_pneum(W_TO, q_max) : 
    
    W_t = 56.01 * (W_TO * q_max * 10**-5)**0.576
    
    return W_t

q_max = 600 # Maximum dynamic pressure, in pounds per square foot (lb/ft2).

# ------------------------------
# Section 20.2.1.4 : Instruments
# ------------------------------

def instruments (N_pil, W_TO, N_e) : 
    
    # Flight Instrument  Indicators :
    W_t_1 = N_pil * (15.0 + 0.032 * (W_TO*10**-3))
        
    # Engine Instrument Indicators (Turbine Engines):
    W_t_2 = N_e * (4.80 + 0.06 *(W_TO * 10**-3))
        
    # Miscellaneous Indicators :
    W_t_3 = 0.15*(W_TO* 10**-3)
    
    W_t = W_t_1 + W_t_2 + W_t_3
    
    return W_t

N_pil = 2 # Number of pilots

# ------------------------------------
# Section 20.2.1.5 : Electrical system
# ------------------------------------

# The weight prediction relationships are expressed in terms of the total 
# weight of the fuel system plus the total weight of the electronics system, the 
# prime users of electrical power on most aircraft.

def electrical_syst(W_fs,W_tron) : 
    
    W_t = 1162.66 * (W_fs * W_tron * 10**-3)**0.506
    
    return W_t

W_fs = prop_sub_syst (N_i, L_d, A_i, P_2, F_gw, F_gf, K_eco, L_f, N_e, W_eng)[1]   # Weight of fuel system, in pounds (lb)
W_tron = 1000 # Weight of electronics system, in pounds (lb)

# ------------------------------
# Section 20.2.1.6 : Furnishings
# ------------------------------

def furninshings (N_fds, N_pass, K_lav, K_buf, P_c, W_TO) :
    
    # Flight Deck Seats Executive and Commercial :
    W_t_1 = 54 * N_fds
    
    # Passenger Seats Executive and Commercial :
    W_t_2 = 32.03 * N_pass 
    
    # Lavatories and Water Provisions Executive and Commerical:
    W_t_3 = K_lav * N_pass**1.33
    
    #  Food Provisions Executive and Commercial :
    W_t_4 = K_buf * N_pass**1.12
    
    # Oxygen system : 
    W_t_5 = 109.33 * (N_pass * (1 + P_c)* 10**-2)**0.505
    
    # Baggage and Cargo Handling Provinsions 
    # I don't take that into account because I don't know what it is
    # W_t_ = K_cbc * N_pass**1.456 
    # K_cbc = 0.0646 without preload provisions; = 0.316 with preload provisions
    
    # Miscellaneous Furnishings and Equipment : 
    W_t_6 = 0.771 * (W_TO*10**-3)
    
    W_t = W_t_1 + W_t_2 + W_t_3 + W_t_4 + W_t_5 + W_t_6
    
    return W_t

N_fds = 1    # Number of flight deck stations
N_pass = 19   # Number of passengers
K_lav = 3.90 # = 3.90 for executive; = 1.11 for long-range commercia1 passenger; = 0.31 for short-range commercial passenger
K_buf = 5.68 # = 5.68 for long-range (707, 990, 737, 747, 757, 767, 777, etc.); = 1.02 for short-range (340, 202, Citation, Learjet, King Air, Jetstream, etc.)
P_c = 8.9      # Ultimate cabin pressure, in pounds per square inch (lb/in^2)

# --------------------------------------------------
# Section 20.2.1.7 : Air Conditioning and Anti-Icing
# --------------------------------------------------

def ac_anti_incing(V_pr, N_cr, N_att, N_pass) : 
    
    W_t = 469.30 * (V_pr * (N_cr + N_att + N_pass) * 10**-4)**0.419
    
    return W_t

V_pr = 3500 # Pressurized or occupied volume, in cubic feet (ft3)
N_att = 2 # Number of attendants

# -----------------------------------------
# Section 20.2.1.8 : Electronics (Avionics)
# -----------------------------------------

def avionics () : 
    
    W_t = 1200 # average for global 8000
    
    return W_t
# -----------------------------------------------
# Section 20.2.1.9 :  Landing Retardation Devices
# -----------------------------------------------

# See chapter 10.5 of the book -> assume to 0 for now 


# Compute individual weights
W_wing = wing_weight(S_w, AR, M_0, W_TO, N, t_c, Lambda_half_chord, taper)
W_horizontal_tail = horizontal_tail_weight(W_TO, N, S_HT, b_HT, t_HT, c_wing, L_t)
W_vertical_tail = vertical_tail_weight(W_TO, N, S_VT, M_0, L_t, AR_VT, lambda_vt, S_r, Lambda_VT, h_rh_lvt_ratio)
W_fuselage = fuselage_weight_usaf(W_TO, K_INL, q, L, H)
W_landing_gear = landing_gear_weight_usaf(W_TO)
W_surf_control = surf_contrl_hydr_pneum(W_TO, q_max)
W_instruments = instruments(N_pil, W_TO, N_e)
W_electrical = electrical_syst(W_fs, W_tron)
W_furnishings = furninshings(N_fds, N_pass, K_lav, K_buf, P_c, W_TO)
W_ac_anti_icing = ac_anti_incing(V_pr, N_cr=2, N_att=N_att, N_pass=N_pass)
W_avionics = avionics()

# Sum up all the weights
total_weight_lbs = (
    W_wing +
    W_horizontal_tail +
    W_vertical_tail +
    W_fuselage +
    W_landing_gear +
    W_surf_control +
    W_instruments +
    W_electrical +
    W_furnishings +
    W_ac_anti_icing +
    W_avionics
)

# Print individual weights
print(f"Wing Weight: {W_wing:.2f} lbs")
print(f"Horizontal Tail Weight: {W_horizontal_tail:.2f} lbs")
print(f"Vertical Tail Weight: {W_vertical_tail:.2f} lbs")
print(f"Fuselage Weight: {W_fuselage:.2f} lbs")
print(f"Landing Gear Weight: {W_landing_gear:.2f} lbs")
print(f"Surface Control + Hydraulics + Pneumatics Weight: {W_surf_control:.2f} lbs")
print(f"Instruments Weight: {W_instruments:.2f} lbs")
print(f"Electrical System Weight: {W_electrical:.2f} lbs")
print(f"Furnishings Weight: {W_furnishings:.2f} lbs")
print(f"Air Conditioning + Anti-Icing Weight: {W_ac_anti_icing:.2f} lbs")
print(f"Avionics Weight: {W_avionics:.2f} lbs")

# Print the total weight
print(f"Total weight of the aircraft: {total_weight_lbs:.2f} lbs")

# Convert total weight to kilograms
total_weight_kg = total_weight_lbs * 0.453592
# Print the total weight
print(f"Total weight of the aircraft: {total_weight_kg:.2f} kg")