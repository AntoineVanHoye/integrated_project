# -*- coding: utf-8 -*-
"""
Longitudinal state-space matrices for the F-4 aircraft
Source: M.V.Cook, Flight Dynamic Principles.
Written by Greg Dimitriadis, University of Liege, Belgium
February 2024
"""

import numpy as np
from scipy import linalg

# Flight Mach number
Mach=0.6
# Flight altitude in m   
altitude=35000*0.3048       
# Air density in kg/m^3  
rho=0.3809
# Airspeed in m/s
V0=178.0
# Angle of attack in rad
alphae=9.4*np.pi/180
# Flight path angle in rad
gamae=0.0
# Mass in kg
m=17642.0      
# Pitch moment of inertia in kg.m^2
Iy=165669.0  
# Wing area in m^2
S=49.239  
# Mean aerodynamic chord in m
cbarbar=4.889  
# Acceleration due to gravity in m/s^2
g=9.81                
# Dimensional stability derivatives
Xu=0.0076*(1/2*rho*V0*S)
Xw=0.0483*(1/2*rho*V0*S)
Xwdot=0.0*(1/2*rho*S*cbarbar)
Xq=0.0*(1/2*rho*V0*S*cbarbar)
Xeta=0.0618*(1/2*rho*V0**2*S)
Xtau=0.0
Zu=-0.7273*(1/2*rho*V0*S)
Zw=-3.1245*(1/2*rho*V0*S)
Zwdot=-0.3997*(1/2*rho*S*cbarbar)
Zq=-1.2109*(1/2*rho*V0*S*cbarbar)
Zeta=-0.3741*(1/2*rho*V0**2*S)
Ztau=0.0
Mu=0.0340*(1/2*rho*V0*S*cbarbar)
Mw=-0.2169*(1/2*rho*V0*S*cbarbar)
Mwdot=-0.5910*(1/2*rho*S*cbarbar**2)
Mq=-1.2732*(1/2*rho*V0*S*cbarbar**2)
Meta=-0.5581*(1/2*rho*V0**2*S*cbarbar)
Mtau=0.0*cbarbar
# Gravity angle
thetae=alphae+gamae;
# Equilibrium airspeeds
Ue=V0*np.cos(thetae);
We=V0*np.sin(thetae);
# Mass matrix
M=np.array([(m, -Xwdot, 0.0, 0.0),
            (0.0, m-Zwdot, 0.0, 0.0),
            (0.0,  -Mwdot, Iy, 0.0),
            (0.0, 0.0, 0.0, 1.0,)])
# Stiffness-damping matrix
K=np.array([(-Xu, -Xw, -(Xq-m*We), m*g*np.cos(thetae)),
            (-Zu, -Zw, -(Zq+m*Ue), m*g*np.sin(thetae)),
            (-Mu,  -Mw, -Mq, 0.0),
            (0.0, 0.0, -1.0, 0.0)])
# Forcing matrix
F=np.array([(Xeta, Xtau),
            (Zeta, Ztau),
            (Meta, Mtau),
            (0.0, 0.0)])
# State-space matrices
A=-linalg.solve(M,K)
B=linalg.solve(M,F)

print("Xu=", Xu)
print("Xw=", Xw)
print("Xwdot=", Xwdot)
print("Xq=", Xq)
print("Xeta=", Xeta)
print("Xtau=", Xtau)
print("Zu=", Zu)
print("Zw=", Zw)
print("Zwdot=", Zwdot)
print("Zq=", Zq)
print("Zeta=", Zeta)
print("Ztau=", Ztau)
print("Mu=", Mu)
print("Mw=", Mw)
print("Mwdot=", Mwdot)
print("Mq=", Mq)
print("Meta=", Meta)
print("Mtau=", Mtau)