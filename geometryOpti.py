import numpy as np
import matplotlib.pyplot as plt

import wings as w
import stability as s

# ---------------------------------------- #
# Constants
config = 3
fuel = 1
Cm0_wing = w.getAirfoilWing()[4]
Cm0_fus = w.getAirfoilFus()[4]
# ---------------------------------------- #

AR = np.linspace(4, 6, 20)
AR = np.flip(AR)
sweep_fus = np.linspace(30, 60, 30)
sweep_wing = np.linspace(20, 40, 20)

results = np.array([0, 0, 0, 0, 0, 0])
for i in range(len(AR)):
    for j in range(len(sweep_fus)):
        for k in range(len(sweep_wing)):
            stability = s.long_stat_stab_cruise(config,fuel,Cm0_fus,Cm0_wing, AR[i], sweep_fus[j], sweep_wing[k]) 
            if 5 < stability < 10:
                force = s.CL(config,fuel,Cm0_fus,Cm0_wing, AR[i], sweep_fus[j], sweep_wing[k])[0]
                cl = w.getCl(AR[i], sweep_fus[j], force)
                setting_angle = w.getCalageAngle(cl, AR[i], sweep_fus[j], sweep_wing[k])[0]
                
                if (setting_angle*180/np.pi) <= 6:
                    tmp = [AR[i], sweep_fus[j], sweep_wing[k], stability, s.CL(config,fuel,Cm0_fus,Cm0_wing, AR[i], sweep_fus[j], sweep_wing[k])[0], s.CL(config,fuel,Cm0_fus,Cm0_wing, AR[i], sweep_fus[j], sweep_wing[k])[1]]
                    results = np.vstack([results, tmp])
                    print(tmp)
results = np.delete(results, 0, 0)
print(np.shape(results))
#print(resutls)