import numpy as np
import matplotlib as plt
import math

x_cg = 10.4
thetamin_lat_clear = 5
theta_LOF = 15

C_r = 4.622
C_t = 1.448
#C_R=  (S)*(C_r-C_t)/9

M = 2/3*(C_r + C_t-(C_r*C_t)/(C_r+C_t)) #for the standard straight tapered wing
#M = 1/2 * (C_r+C_t)
b = 20
H = b/8*(C_r**2+2*C_r*C_t+3*C_t**2)/(C_r**2+C_r*C_t+C_t**2)


print(M,H)


