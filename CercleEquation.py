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


print(distance_horizontale(3.0, 40, 0.5)) 
print(distance_horizontale(3.0, 40, 1)) 
print(distance_horizontale(3.0, 40, 1.5)) 
print(distance_horizontale(3.0, 40, 2.0)) 
print(distance_horizontale(3.0, 40, 2.5)) 
print(distance_horizontale(3.0, 40, 3.0))
print(distance_horizontale(3.0, 40, 3.5)) 