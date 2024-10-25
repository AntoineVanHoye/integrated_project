import numpy as np
import math

#Starting Variable
W0 = 106000 #Pound
c = 0.45 #Pound/hours/Pound
LbyD = 18
AR = 8
R = 9206.236 #range in miles
V = 647.026 #speed in miles per hours (M = 0.85 )
Pl_1 = 2120 # payload 1
Pl_2 = 1260 # payload 2
Pl_3 = 2520 # payload 3
W_crew = 690 # crew weight
A = 0.2678 # Roskans values of A in the table
B = 0.9979 # Roskans values of B in the table

def initandfuel_weight():
    W1 = 0.97 * W0  #Textbook data (starup and taxi)
    W2 = 0.985* W1 #Textbook data (climb)
    W3 = W2*math.exp(-(R*c)/(V*LbyD)) #Breget Range Equation (cruise)
    W4 = W3 * 0.995 #Textbook data (landing and taxi)

    Mff = (W1/W0)*(W2/W1)*(W3/W2)*(W4/W3)

    W_used = (1- Mff)* W0 # used fuel weight
    W_res = 0.1*W_used # reserved fuel weight
    Wf = W_used + W_res # fuel weight

    #W_oe1 = W0-Wf-Pl_1 #operational emply weight (for mission 1 to 3)
    #W_oe2 = W0-Wf-Pl_2
    #W_oe3 = W0-Wf-Pl_3

    #W_e1 = W_oe1-W_crew
    #W_e2 = W_oe2-W_crew
    #W_e3 = W_oe3-W_crew

    We = 10**((math.log10(W0)-A)/B) #empty weight

    return We , Wf
