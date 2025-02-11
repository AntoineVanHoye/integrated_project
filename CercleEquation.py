import numpy as np
from scipy.optimize import fsolve

def equation_cercle(x, y, R):
    """Équation du cercle centrée en (2R, R)"""
    return (x - 9.006)**2 + (y - R)**2 - R**2

def distance_horizontale(R, lambda_deg, y):
    """Calcule la distance horizontale entre la droite inclinée et le cercle pour une hauteur y"""
    lambda_rad = np.radians(lambda_deg)  # Conversion en radians
    
    # Coordonnée x sur la droite inclinée pour y donné
    x_droite = y / np.tan(lambda_rad)
    
    # Trouver l'intersection avec le cercle en résolvant l'équation
    x_cercle_solution = fsolve(equation_cercle, 2*R, args=(y, R))[0]
    
    # Distance horizontale cherchée
    d = x_cercle_solution - x_droite
    return d

h = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0] #, 2.2, 2.4, 2.6, 2.8, 3.0]#[0.5, 1, 1.5, 2, 2.5, 3]
R = 3.0
c = np.zeros(len(h))
for i in range(len(h)):
    c[i] = distance_horizontale(R,25, h[i])

#print(c)
print(c.tolist())