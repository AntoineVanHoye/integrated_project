import numpy as np
import math

# Design Variables 
W0 = 106000  # MTOW (lbs)
c = 0.5  # TSFC (lb/hr/lb)
LbyD = 25  # Lift-to-Drag Ratio (assumes for BWB)
AR = 3  # Aspect Ratio (Assumed for BWB)
R = 9206.236  # Range (miles)
V = 647.026  # Cruise Speed (mph)
W_crew = 690  # Crew Weight (lbs)
A = 0.2678  # Roskam's A
B = 0.995  # Roskam's B
payloads = [2120, 1260, 2520]  # Payload variations (lbs) for the diffrent missions 


def breguet_range_equation(W_initial, c, LbyD, R, V):  #Breguet Range Equation to compute fuel weight fraction

    return W_initial * math.exp(-(R * c) / (V * LbyD))

def initandfuel_weight(W0, c, LbyD, R, V, W_crew, A, B): #Calculate weights for various phases and fuel usage
    
    # Weight during various phases
    W1 = 0.97 * W0  # Startup and taxi
    W2 = 0.985 * W1  # Climb
    W3 = breguet_range_equation(W2, c, LbyD, R, V)  # Cruise
    W4 = 0.995 * W3  # Landing and taxi

    # Mission Fuel Fraction (Mff)
    Mff = (W1 / W0) * (W2 / W1) * (W3 / W2) * (W4 / W3)

    # Fuel Weights
    W_used = (1 - Mff) * W0  # Used fuel weight
    W_res = 0.1 * W_used  # Reserved fuel weight
    Wf = W_used + W_res  # Total fuel weight

    # Landing Weight
    W_landing = W3  # Weight after cruise, before landing phase

    # Operational Empty Weight and Zero Fuel Weight
    We = 10 ** ((math.log10(W0) - A) / B)  # Empty weight (Roskam's method)
    W_oe = W0 - Wf  # Operational empty weight
    W_zf = W_oe + W_crew  # Zero fuel weight (includes crew weight)

    return {
        "Empty Weight (We)": We,
        "Fuel Weight (Wf)": Wf,
        "Used Fuel Weight": W_used,
        "Reserve Fuel Weight": W_res,
        "Landing Weight": W_landing,
        "Operational Empty Weight (W_oe)": W_oe,
        "Zero Fuel Weight (W_zf)": W_zf
    }

# Payload-Range Analysis
def payload_range_analysis(W0, c, LbyD, R, V, W_crew, payloads):
    """Analyze weights and fuel for different payloads."""
    results = []
    for payload in  payloads:
        W0_adj = W0 - payload - W_crew
        weights = initandfuel_weight(W0_adj, c, LbyD, R, V, W_crew, A, B)
        weights["Payload"] = payload
        results.append(weights)
    return results

weights = initandfuel_weight(W0, c, LbyD, R, V, W_crew, A, B)
print("General Weight Breakdown:")
for key, value in weights.items():
    print(f"{key}: {value:.2f} lbs")

print("\nPayload-Range Analysis:")
payload_results = payload_range_analysis(W0, c, LbyD, R, V, W_crew,  payloads)
for idx, result in enumerate(payload_results):
    print(f"\nPayload Case {idx + 1}:")
    for key, value in result.items():
        print(f"{key}: {value:.2f} lbs")
