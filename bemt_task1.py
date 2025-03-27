import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from mat4py import loadmat

#from the Clark-Y polar function given 
def clarkypolarsRe(aoa, re):

  POLAR_DATA_FILE = 'utils/clarkypolars.mat'

  # Load hard-coded polarData
  tmp = loadmat(POLAR_DATA_FILE)
  polarData = tmp['polarData']

  aoa = np.arctan2(np.sin(aoa), np.cos(aoa))

  # 2D interpolation
  cl_interpolator = RegularGridInterpolator((np.array(polarData['aoa'])[:, 0], np.array(polarData['reynolds'])), np.array(polarData['cl']), bounds_error=False, fill_value=None)
  cd_interpolator = RegularGridInterpolator((np.array(polarData['aoa'])[:, 0], np.array(polarData['reynolds'])), np.array(polarData['cd']), bounds_error=False, fill_value=None)

  cl = cl_interpolator((aoa, re))
  cd = cd_interpolator((aoa, re))

  return cl, cd

#Given 
rho = 1.225 # kg/m^3
V_inf =20 #m/s
n = 150 #revs/s (9000 RPM)
Omega = 2* np.pi *n # rad/s rotational velocity
B = 2 # number of blades 
nu = 1.568e-5 # kinematic viscosity (m^2/s)

#propeller data CONVERTED TO m
props = {
  '9x4.5' : {'D':9*0.254, 'pitch': 4.5*0.0254},
  '9x6' : {'D':9*0.254, 'pitch': 6*0.0254},
  '11x7' : {'D':11*0.254, 'pitch': 7*0.0254},
  '11x10' : {'D':11*0.254, 'pitch': 10*0.0254}
}

#automatically load the geometies from the base directory 
def load_geometry (prop_name):
  file_path = f'apce_{prop_name}_geom.txt' #retrives the file from the base directory
  try:
    data =  np.loadtxt(file_path, skiprows =1)
    r_R = data[:,0] # was normalised
    c_in = data[:,1]
    twist_deg = data[:,2]

    R = props[prop_name]['D']/2
    r = r_R*R
    c = c_in *0.0254 # cpnverted to m
    twist =  twist_deg

    print(f"{prop_name}:loaded {len (r)} points, r[0]={r[0]:.4f} m, c[0]={c[0]:.4f} m, twist[0]={twist[0]:.2f} deg")
    return r,c, twist
  except FileNotFoundError:
        print(f"Error:{file_path} not found in directory")
        raise
  except Exception as e:
     print (f"Error loading {file_path}: {str(e)}")
     raise

# Bemt cmputation 
def bemt_analysis(D, pitch, r, c, twist):
   dr = r[1] - r[0]
   T_total = 0 
   Q_total = 0
   dT_dr = np.zeros_like(r)
   dQ_dr = np.zeros_like (r)

   for i in range (len(r)):
      a = 0.1 # this is an initial guess value till i get the main axial induce value via convergence 
      a_prime = 0.01 # this is an initial guess value till i get the main tangential induce value via convergence
      tol = 1e-6
      max_iter = 100

      for _ in range (max_iter) : 
         V_axial = V_inf * (1+a)
         V_tang = Omega*r[i]*(1-a_prime)
         V_r = np.sqrt(V_axial**2 + V_tang)
         phi = np.arctan2(V_axial, V_tang)
         alpha = np.radians(twist[i])- phi 

         Re = V_r*c[i]/nu 
         Cl, Cd = clarkypolarsRe(alpha, Re)

         dF_T = 0.5 *rho * V_r**2 * c[i] * (Cl *np.cos(phi)-Cd*np.sin(phi)) # thrust force per unit span 
         dF_Q = 0.5 * rho *V_r**2 * c[i] * (Cl * np.sin(phi) + Cd * np.cos(phi)) # torque force 

         sigma = B * c[i]/ (2 *np.pi * r[1]) # local solidity 
         a_new = dF_T/ (4*np.pi * r[i]* rho *V_inf**2 * (1+a))
         a_prime_new = dF_Q / (4 * np.pi * rho * r[i] * V_inf**2 * (1 + a) * Omega * r[i])
         if a_new > 0.5:
            a_new = 0.5
            a = 0.7 *a + 0.3*a_new # numerical method in solving implicit equations iteratively (1-w)*x + w*y
            a_prime = 0.7*a_prime + 0.3*a_prime_new 

            if abs(a - a_new) < tol and abs(a_prime - a_prime_new) < tol:
               break 
            
            dT_dr[i] = B * dF_T/dr
            dQ_dr[i] = B * dF_Q * r[i] / dr
            T_total += B * dF_T * dr 
            Q_total += B * dF_Q * r[i] * dr 

            if i == 0:
               print (f"{prop}: r= {r[i]*39.37:.2f} in, alpha={np.degrees (alpha):.2f} deg, Cl ={Cl:.3f}, Cd={Cd:.3f, dT/dr={dT_dr[i]:.2f} N/m}")

               P_total = Omega * Q_total
               eta = T_total * V_inf / P_total if P_total > 0 else 0 
               return T_total, P_total, eta, dT_dr, dQ_dr
            
            #Results 
            results = {}
            for prop in props:
               r, c, twist = load_geometry(prop)
               T, P, eta, dT_dr, dQ_dr = bemt_analysis(props[prop]['D'],props[prop]['pitch'], r, c, twist)
               results[prop] = {'T': T, 'P': P, 'eta': eta, 'r': r, 'dT_dr': dT_dr, 'dQ_dr': dQ_dr}
               print(f"{prop}: T = {T:.2f} N, P = {P:.2f} W, eta = {eta:.3f}")

               #Visual representation 
               plt.figure(figsize=(12,5))
               for prop in props:
                  plt.subplot(1, 2, 1)
                  plt.plot(results [prop]['r'] *39.37, results[prop]['dT_dr'],label=prop)
                  plt.xlabel('Radius (in)')
                  plt.ylabel('dT/dr (N/m)')
                  plt.title ('Thrust Distribution')
                  plt.legend()
                  plt.grid(True)

                  plt.subplot(1, 2, 2)
                  plt.plot(results[prop]['r'] * 39.37, results [prop] ['dQ_dr']* Omega, label=prop)
                  plt.xlabel('Radius (in)')
                  plt.ylabel('dP/dr (W/m)')
                  plt.title('Power Distribution ')
                  plt.legend ()
                  plt.grid(True)

                  plt.tight_layout()
                  plt.show()

               
            

         






  
  
  