import numpy as np
import matplotlib.pyplot as plt

import wings as w
import static_stability as s

# ---------------------------------------- #
# Constants
config = 3
fuel = 2
Cm0_wing = w.getAirfoilWing()[4]
Cm0_fus = w.getAirfoilFus()[4]
# ---------------------------------------- #

AR = np.linspace(3, 5, 21)
AR = np.flip(AR)
sweep_fus = np.linspace(15, 50, 36)
sweep_wing = np.linspace(25, 35, 11)

results = np.array([0, 0, 0, 0, 0, 0, 0])
for i in range(len(AR)):
    for j in range(len(sweep_fus)):
        for k in range(len(sweep_wing)):

            stability, _ = s.long_stat_stab_cruise(config,fuel,Cm0_fus,Cm0_wing, AR[i], sweep_fus[j], sweep_wing[k]) 
            if 5 < (stability*100) < 15:
                stability, _ = s.long_stat_stab_cruise(config,1,Cm0_fus,Cm0_wing, AR[i], sweep_fus[j], sweep_wing[k])
                if 5 < (stability*100) < 15:

                    force = s.CL(config,fuel,Cm0_fus,Cm0_wing, AR[i], sweep_fus[j], sweep_wing[k])[0]
                    cl = w.getCl(AR[i], sweep_fus[j], force)
                    setting_angle = w.getCalageAngle(cl, AR[i], sweep_fus[j], sweep_wing[k])[0]

                    if (setting_angle*180/np.pi) <= 5.5:
                        print("Bad  ", [AR[i], sweep_fus[j], sweep_wing[k], stability*100, setting_angle*180/np.pi])

                    if (setting_angle*180/np.pi) <= 4.5:
                        tmp = [AR[i], sweep_fus[j], sweep_wing[k], stability*100, s.CL(config,fuel,Cm0_fus,Cm0_wing, AR[i], sweep_fus[j], sweep_wing[k])[0], s.CL(config,fuel,Cm0_fus,Cm0_wing, AR[i], sweep_fus[j], sweep_wing[k])[1], setting_angle*180/np.pi]
                        results = np.vstack([results, tmp])
                        print(tmp)

results = np.delete(results, 0, 0)
print(np.shape(results))
with open("config_possible.txt", "w") as file:
        for resutlats in results:  # Boucle sur les valeurs des arrays
            file.write(str(resutlats) + "\n")
#print(results)

# ------ Optimisation manuel -------- #

#Airfoil ? 
AR = 4
sweep_fus = 50
sweep_wing = 30
Cabin_width = 7
Cabin_length = 18.8
#Position tail ? 

force = 694123.11242362


#w.printFunction(AR, sweep_fus, sweep_wing, force)
#s.printFunction(AR, sweep_fus, sweep_wing)