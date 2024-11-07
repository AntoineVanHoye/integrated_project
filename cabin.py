#estimation of the internal geometry of the cabin area
import math
import numpy as np
#Reference dimention (taken from the area of the bombardier cabin)

height = 1.88
width = 2.43
lenght = 16.59
theta = 2*math.acos((height-width/2)/(width/2))
#alpha=


Total_vol_estim = lenght*((np.pi-theta/2)*(width/2)**2+1/2*(width/2)**2*math.sin(theta))
Total_vol_rough = 16.59*2.43*1.88
print(Total_vol_estim)
print(Total_vol_rough)
print(Total_vol_estim/1.88)
