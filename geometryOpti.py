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

AR = np.linspace(2, 6, 30)
sweep_fus = np.linspace(30, 60, 30)
sweep_wing = np.linspace(20, 40, 20)

resutls = np.array([0, 0, 0, 0, 0, 0])
for i in range(len(AR)):
    for j in range(len(sweep_fus)):
        for k in range(len(sweep_wing)):
            stability = s.long_stat_stab_cruise(config,fuel,Cm0_fus,Cm0_wing, AR[i], sweep_fus[j], sweep_wing[k]) 
            if 5 < stability < 10:
                tmp = [AR[i], sweep_fus[j], sweep_wing[k], stability, s.CL(config,fuel,Cm0_fus,Cm0_wing, AR[i], sweep_fus[j], sweep_wing[k])[0], s.CL(config,fuel,Cm0_fus,Cm0_wing, AR[i], sweep_fus[j], sweep_wing[k])[1]]
                resutls = np.vstack([resutls, tmp])
                print(tmp)
resutls = np.delete(resutls, 0, 0)
print(np.shape(resutls))
#print(resutls)