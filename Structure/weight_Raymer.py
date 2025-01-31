# Formulas from Daniel Raymer book, results are weird !


import math

# Formulas 

# Wing Weight
def raymer_wing_weight(S_w, AR, sweep, taper_ratio, t_c, n_z, W_dg, S_csw):
    factor = 0.0051 * (W_dg * n_z)**0.557 * S_w**0.649 * AR**0.5 * t_c**-0.4
    AR_factor = (1 + taper_ratio)**0.1 * (math.cos(math.radians(sweep)))**-1.0 * S_csw**0.1
    return factor * AR_factor

# Horizontal Tail Weight
def raymer_ht_weight(n_z, W_dg, S_ht, sweep_ht, AR_ht, t_c, K_uht, F_w, B_h, L_t, K_y, S_e):
    factor = 0.0379 * K_uht * (1 + F_w / B_h)**-0.25 * W_dg**0.639 * n_z**0.10
    S_factor = S_ht**0.75 * L_t**-1.0 * K_y**0.704 * (math.cos(math.radians(sweep_ht)))**-1
    AR_factor = AR_ht**0.166 * (1 + S_e / S_ht)**0.1
    return factor * S_factor * AR_factor

# Vertical Tail Weight
def raymer_vt_weight(n_z, W_dg, S_vt, sweep, AR_vt, t_c_root, H_t, H_v, K_z, L_t):
    factor = 0.0026 * (1 + H_t / H_v)**0.225 * n_z**0.536 * W_dg**0.556 * L_t**-0.5
    S_factor = S_vt**0.5 * K_z**0.875 * (math.cos(math.radians(sweep)))**-1 * AR_vt**0.35 * t_c_root**-0.5
    return factor * S_factor

# Fuselage Weight
def raymer_fuselage_weight(S_f, n_z, W_dg, L, D, K_door, K_Lg, K_ws):
    factor = 0.3280 * K_door * K_Lg * (n_z * W_dg)**0.5 * L**0.25 * S_f**0.302
    geometry_correction = (1 + K_ws)**0.04 * ((L / D)**0.10)
    return factor * geometry_correction

# Main Landing Gear Weight
def raymer_mlg_weight(K_mp, W_l, N_l, L_m, N_mw, N_mss, V_stall):
    factor = 0.0106 * K_mp * (W_l)**0.888 * N_l**0.25 * L_m**0.4
    wheel_factor = N_mw**0.321 * N_mss**-0.5 * V_stall**0.1
    return factor * wheel_factor

# Nose Landing Gear Weight
def raymer_nlg_weight(K_np, W_l, N_l, L_n, N_nw):
    factor = 0.032 * K_np * W_l**0.646 * N_l**0.2
    length_factor = L_n**0.5 * N_nw**0.45
    return factor * length_factor

# Nacelle Group Weight
def raymer_nacelle_group_weight(K_ng, N_Lt, N_w, n_z, W_ec, N_eng, S_n):
    return 0.6724 * K_ng * N_Lt**0.1 * N_w**0.294 * n_z**0.119 * W_ec**0.611 * N_eng**0.984 * S_n**0.224

# Engine Controls Weight
def raymer_engine_controls(N_eng, L_ec):
    return 5 * N_eng + 0.8 * L_ec

# Engine Starter Weight (Pneumatic)
def raymer_engine_starter_pneumatic(N_eng, W_en):
    return 49.19 * ((N_eng * W_en) / 1000)**0.541

# Fuel System Weight
def raymer_fuel_system_weight(V_t, V_i, V_p, N_t):
    factor = 2.405 * V_t**0.606
    volume_ratio = (1 + V_i / V_t)**-1 * (1 + V_p / V_t)
    tank_factor = N_t**0.5
    return factor * volume_ratio * tank_factor

# Flight Control System Weight
def raymer_flight_control_weight(N_f, N_m, S_cs, I_yaw):
    factor = 145.9 * N_f**0.554
    mechanical_factor = (1 + N_m / N_f)**-1
    surface_area_factor = S_cs**0.2 * (I_yaw * 10**-6)**0.07
    return factor * mechanical_factor * surface_area_factor

# APU Weight
def raymer_APU_weight(W_APU_uninstalled):
    return 2.2 * W_APU_uninstalled

# Instruments Weight
def raymer_instrument_weight(K_r, K_tp, N_c, N_eng, L_f, B_w):
    return 4.509 * K_r * K_tp * N_c**0.541 * N_eng * (L_f + B_w)**0.5

# Hydraulic System Weight
def raymer_hydraulic_system_weight(N_f, L_f, B_w):
    return 0.2673 * N_f * (L_f + B_w)**0.937

# Avionics Weight
def avionics_weight(W_uav):
    return 1.73 * W_uav**0.983

# Electrical System Weight
def electrical_system_weight(R_kva, L_a, N_gen):
    return 7.291 * R_kva**0.782 * L_a**0.346 * N_gen**0.1

# Air-conditioning and Anti-icing Weight
def air_conditioning_weight(N_p, V_pr, W_uav, W_dg):
    return 62.36 * N_p**0.25 * (V_pr / 1000)**0.604 * W_uav**0.1 + 0.002 * W_dg

# Furnishings Weight
def raymer_furnishings_weight(N_c, W_c, S_f):
    return 0.0577 * N_c**0.1 * W_c**0.393 * S_f**0.75


# Main Function to compute all components
def main():
    # Inputs for a long-range business jet
    inputs = {
        "W_dg": 106000,       # Design gross weight (lb)
        "n_z": 1.5,           # Load factor
        "S_w": 715.15421,     # Wing area (ft²)
        "AR": 1.5,            # Wing aspect ratio
        "sweep": 32.123,      # Sweep angle (degrees) at 25% MAC
        "taper_ratio": 0.161,
        "t_c": 0.107,         # thickness-to-chord ratio (if not constant, use average of portion of wing inboard of C-bar) 
        "S_csw": 75,          # Control surface area, ft2 (wing-mounted, includes flaps)
        
        "S_ht": 34.2292400001119*2,         # Horizontal tail area (ft²)
        "AR_ht": 2.5,                       # Horizontal tail aspect ratio
        "sweep_ht": 30,                     # Horizontal tail sweep (degrees)
        "K_uht": 1.0,                       # 1.143 for unit (all-moving) horizontal tail; = 1.0 otherwise
        "F_w": 29.5276,                     # Fuselage witdh at intersection with the tail
        "B_h": 12.080052,                   # HT span
        "L_t": 30,                          # Tail length; wing quarter-MAC to tail quarter-MAC, ft 
        "K_y": 0.3*30,                      # Aircraft pitching radius of gyration, ft(= 0.3 Lt) 
        "S_e": (34.2292400001119*2*10/100), # Elevator area (ft^2) (j'approx 10% de la surface)
        "S_vt" : 100,                       # Vetical tail area (ft^2) 
        "AR_vt" : 2.5,                      # Aspect ratio of the vertical tail 
        "t_c_root" :  0.107,                # Thickness-to-chord ratio (pas précisé si c'est celui de la tail ou de l'aile donc je prends celui de l'aile => t_c)
        "H_t" :  27,                        # Horizontal tail height above fuselage, ft, valeur approximée (dans le livre y a aussi le ratio H_t / H_v = 0.0 for conventional tail; 1.0 for ''T" tail)
        "H_v" :  30,                        # Vertical tail height above fuselage, ft 
        "K_z" :  30,                        # Aircraft yawing radius of gyration, ft (~= Lt) 
        
        "S_f": 4080.759898,         # Fuselage wetted area (ft²)
        "L" : 55.104,               # Fuselage structural length, ft (excludes radome cowling, tail cap)  
        "D": 22.9659,               # Fuselage diameter (ft), j'approxime à 7m
        "K_door": 1.0,              # No cargo door
        "K_Lg": 1.12,               # Because fuselage mounted main landing gear
        "K_ws": 0.7,                # Wing sweep factor (There is a formula to implement !!!)
        
        "K_mp" : 1.0,               # 1.126 for kneeling gear; = 1.0 otherwise (I don't really know what's the kneeling gear)
        "W_l" : 106000,             # Landing design gross weight, lb (same as W_dg)
        "N_l" : 4.05,               # Ultimate landing load factor; = N_gear x 1.5 (Chat GPT : N_gear = The design gear load factor = 2.7. Here, = (2.7×1.5))
        "L_m" : 15 ,                # Extended length of main landing gear, (in). 
        "N_mw" : 4,                 # Number of main wheels
        "N_mss" : 2,                # Number of main gear shock struts (Shock struts absorb the impact during landing)
        "V_stall" : 185.7,          # Stall velocity in ft/s (assumed and converted from ~110 knots)
        "K_np" : 1.0,               # 1.15 for kneeling gear (C-5); = 1.0 otherwise
        "L_n" :  50,                # extended nose gear length, in. (aprroximation)
                  
        "V_t": 6000,     # Total fuel volume (gal)
        "V_i": 4500,     # Integral tank volume (gal)
        "V_p": 500,      # Protected tank volume (gal)
        "N_t": 3,        # Number of fuel tanks
        
        "W_APU_uninstalled": 400, # APU weight (lb)
        
        "N_m" : 0,               # Number of number of surface controls driven by mechanical actuation instead of hydraulics (must be <= N_f and is typically 0-3) 
        "S_cs" : 100,            # Total area of control surfaces, ft2 ()
        "I_yaw" : 12000000 ,     # Yawing moment of inertia, lb-ft2 (see Chap. 16, p.655) 
            
        "K_r" : 1.0,            # 1.133 if reciprocating engine; = 1.0 otherwise
        "K_tp" : 1.0,           # 0. 793 if turboprop; = 1.0 otherwise (we have turbofan)
        "N_eng" : 2.0,          # Number of engines 
        
        "N_f" : 7,          # Number of separate functions performed by surface cont. (typically 4-7)
        "L_f" : 55.11811,   # Total fuselage length
        "B_w" : 95.1444,    # Wing span (ft)
        
        "W_uav" : 1000,     # Uninstalled avionics weight, lb (typically = 800-1400 lb) 
        "N_c": 10,       # Number of seats
        "W_c": 80000,    # Weight of crew and furnishings (lb)
        "R_kva": 50,     # Electrical rating (kVA)
        "L_a": 100,      # Length of avionics components (ft)
        "N_gen": 2,      # Number of generators
        "N_p" : 11,      # Number of personnel on board (crew and passengers)
        "V_pr" : 500,    # volume of pressurized section, ft^3 
    }

    # Calculate all weights
    weights = {
        "Wing Weight": raymer_wing_weight(inputs["S_w"], inputs["AR"], inputs["sweep"], inputs["taper_ratio"], inputs["t_c"], inputs["n_z"], inputs["W_dg"], inputs["S_csw"]),
        "Horizontal Tail Weight": raymer_ht_weight(inputs["n_z"], inputs["W_dg"], inputs["S_ht"], inputs["sweep_ht"], inputs["AR_ht"], inputs["t_c"], inputs["K_uht"], inputs["F_w"], inputs["B_h"], inputs["L_t"], inputs["K_y"], inputs["S_e"]),
        "Vertical Tail Weight": raymer_vt_weight(inputs["n_z"], inputs["W_dg"], inputs["S_vt"], inputs["sweep"], inputs["AR_vt"], inputs["t_c_root"], inputs["H_t"], inputs["H_v"], inputs["K_z"], inputs["L_t"]),
        "Fuselage Weight": raymer_fuselage_weight(inputs["S_f"], inputs["n_z"], inputs["W_dg"], inputs["L"], inputs["D"], inputs["K_door"], inputs["K_Lg"], inputs["K_ws"]),
        "Main Landing Gear Weight": raymer_mlg_weight(inputs["K_mp"], inputs["W_dg"], inputs['N_l'], inputs["L_m"], inputs["N_mw"], inputs["N_mss"], inputs["V_stall"]),
        "Nose Landing Gear Weight": raymer_nlg_weight(inputs["K_np"], inputs["W_l"], inputs['N_l'], inputs["L_n"], inputs["N_mw"]),
        "APU Weight": raymer_APU_weight(inputs["W_APU_uninstalled"]),
        "Fuel System Weight": raymer_fuel_system_weight(inputs["V_t"], inputs["V_i"], inputs["V_p"], inputs["N_t"]),
        "Flight Control Weight": raymer_flight_control_weight(inputs["N_f"], inputs["N_m"], inputs["S_cs"], inputs["I_yaw"]),
        "Instrument Weight": raymer_instrument_weight(inputs["K_r"], inputs["K_tp"], inputs["N_c"], inputs['N_eng'], inputs["L_f"], inputs["B_w"]),
        "Hydraulic System Weight": raymer_hydraulic_system_weight(inputs["N_f"], inputs["L_f"], inputs["B_w"]),
        "Avionics Weight": avionics_weight(inputs["W_uav"]),
        "Electrical System Weight": electrical_system_weight(inputs["R_kva"], inputs["L_a"], inputs["N_gen"]),
        "Air Conditioning Weight": air_conditioning_weight(inputs["N_p"], inputs["V_pr"], inputs["W_uav"], inputs["W_dg"]),
        "Furnishings Weight": raymer_furnishings_weight(inputs["N_c"], inputs["W_c"], inputs["S_f"]),
    }

    # Total Weight Calculation
    total = sum(weights.values())
    weights["Total"] = total

    # Output all weights
    for name, weight in weights.items():
        print(f"{name}: {weight:.2f} lb")


if __name__ == "__main__":
    main()
