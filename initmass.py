import numpy as np
import math

# Design Variables
W0 = 106000  # MTOW (lbs)
c = 0.45  # TSFC (lb/hr/lb)
LbyD = 11.62  # Lift-to-Drag Ratio (assumed for BWB)
R = 8000  # Range (nautical miles)
V = 595.35  # Cruise Speed (knots)
W_fixed = 4100  # Weight of passengers + crew + baggage + avionics + sensors + refreshments

# Conversion factors
LBS_TO_NEWTON = 4.44822  # 1 lb = 4.44822 N

# Breguet Range Equation to compute fuel weight fraction
def breguet_range_equation(c, LbyD, R, V):
    """
    Calculate fuel weight fraction based on the Breguet range equation.
    """
    return math.exp((R * c) / (V * LbyD))

# Initialize weights and fuel weight calculation
def init_and_fuel_weight(W0, c, LbyD, R, V, W_fixed):
    """
    Calculate weights for various phases and fuel usage using the Nicolai method.
    """
    # Weight during various phases
    W1 = 0.97 * W0  # Startup and taxi
    W2 = 0.978 * W1  # Climb
    W3 = breguet_range_equation(c, LbyD, R, V)  # Cruise

    # Mission Fuel Fraction (Mff)
    Mff = (W1 / W0) * (W2 / W1) * (1 / W3)

    # Total fuel weight
    W_f = 1.06 * (1 - Mff) * W0

    # Nicolai method for empty weight fraction
    We_frac = 0.911 / (W0**0.053)  # Empty weight fraction (Nicolai method)
    We_a = W0 * (1 - W_f / W0) - W_fixed  # Available empty weight
    We_r = We_frac * W0  # Required empty weight

    # Operational empty weight (OEW) and zero fuel weight (ZFW)
    W_oe = We_a  # Operational empty weight
    W_zf = W0 - W_f  # Zero fuel weight

    # Maximum Payload Calculation
    MZFW = W_zf  # Maximum Zero Fuel Weight
    max_payload = MZFW - W_oe

    # Landing Weight Calculation
    landing_weight = W_zf  # Assuming no additional fuel burn after cruise

    return {
        "Available Empty Weight (We)": We_a,
        "Fuel Weight (Wf)": W_f,
        "Required Empty Weight (We)": We_r,
        "Operational Empty Weight (W_oe)": W_oe,
        "Zero Fuel Weight (W_zf)": W_zf,
        "Maximum Zero Fuel Weight (MZFW)": MZFW,
        "Maximum Payload": max_payload,
        "Landing Weight": landing_weight,
    }

# Calculate weights and display results
weights = init_and_fuel_weight(W0, c, LbyD, R, V, W_fixed)

print("Weight Breakdown (Nicolai Method):")
for key, value in weights.items():
    value_n = value * LBS_TO_NEWTON  # Convert to Newtons
    print(f"{key}: {value:.2f} lbs ({value_n:.2f} N)")
