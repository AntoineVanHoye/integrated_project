import numpy as np
import matplotlib.pyplot as plt

#geometry wings
beta =  #setting angle
(dxa,dya,dza) =  #distance of the aerodynamic center to the root wing
(dxw,dyw,dzw) =  #y placement of the center of gravity of the wing to the root of the regular wing
(dxa_emp,dya_emp,dza_emp) =  #distance of the aerodynamic center to the root wing of the empenage

#Forces
Lwt =   #lift of the wings tappered ( talk to the professor wether we need to separate the differents kind of lifts)
Lw =    #lift of the wing overall
Dwt =   #drag of the wings tappered ( talk to the professor wether we need to separate the differents kind of drags)
Dw =    #drag of the wing overall
Lfus =  #lift of the fuselage
Dfus =  #drag of the fuselage

P =     #lift empenage (add drag ?)
Ffin =  #lift fin (formulas of the course ?)
T =    #thrust (data from Amos)

#Weight
Wwt =    #weight of the wing tappered
Ww =     #weight of the wing overall
Wfus =  #weight of the fuselage
#Moments
Mwt =     #pitch down moment wings tappered
Mw =      #pitch down moment wing overall
Mfus =   #pitch down moment fuselage
Memp =   #pitch down moment empenage

def wing_loads_regular(n,alpha): #load at the root of the regular wing 


    Tx = (n*Wwt/2-Lwt/2)*np.sin(alpha+beta) + Dwt/2*cos(alpha+beta) #alpha is the angle of attack, beta is the setting angle
    Ty = 0
    Tz = (-n*Wwt/2-Lwt/2)*np.cos(alpha+beta) + Dwt/2*sin(alpha+beta)
    Mx = 1/2*(-n*(Wwt*dyw)+Lwt*dya)*np.cos(alpha+beta) + Dwt/2*dya*np.sin(alpha+beta)#dyw is the y placement of the center of gravity to the root of the regular wing , dya is the distance of the aerodynamic center to the root wing
    My = 1/2*(n*Wwt*dxw -Lwt*dxa +Dwt*dza)*np.cos(alpha) + 1/2(n*Wwt*dxw -Lwt*dxa - Dwt*dza)*np.sin(alpha+beta) + Mw/2 
    Mz = 1/2*(-n*(Wwt*dyw)+Lwt*dya)*np.sin(alpha+beta) - Dwt/2*dya*np.cos(alpha+beta)
    return (Tx,Ty,Tz,Mx,My,Mz)

def wing_loads_curved(n,alpha):

    Tx = (n*Ww/2-Lw/2)*np.sin(alpha+beta) + Dw/2*cos(alpha+beta) #alpha is the angle of attack, beta is the setting angle
    Ty = 0
    Tz = (-n*Ww/2-Lw/2)*np.cos(alpha+beta) + Dw/2*sin(alpha+beta)
    Mx = 1/2*(-n*(Ww*dyw)+Lw*dya)*np.cos(alpha+beta) + Dw/2*dya*np.sin(alpha+beta)#dyw is the y placement of the center of gravity to the root of the regular wing , dya is the distance of the aerodynamic center to the root wing
    My = 1/2*(n*Ww*dxw -Lw*dxa +D*dza)*np.cos(alpha) + 1/2(n*Ww*dxw -Lw*dxa - Dw*dza)*np.sin(alpha+beta) + Mw/2 
    Mz = 1/2*(-n*(Ww*dyw)+Lw*dya)*np.sin(alpha+beta) - Dw/2*dya*np.cos(alpha+beta)
    return 

def inner_wing_load(n,alpha):
    Tx = (n*(Wfus)+T+)*np.cos(alpha)
    Ty =
    Tz =
    Mx =
    My = 
    Mz =
    return