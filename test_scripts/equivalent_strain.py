import numpy as np
from numpy import linalg as LA


#Strain tensor from ARAMIS nodemap export
eps_x = -0.048588741570711 / 100
eps_y = 1.085232496261597 / 100
eps_xy = -0.003783278865740
eps_eqv = 1.301692605018616 / 100

eps = np.array([[eps_x, eps_xy], [eps_xy, eps_y]])

# GOM Aramis implementation of equivalent strain using large strain theory
w, v = LA.eig(eps)
eps_1 = w[0]
eps_2 = w[1]
phi_1 = np.log(1 + eps_1)
phi_2 = np.log(1 + eps_2)
phi_3 = phi_1 + phi_2

phi_M = np.sqrt(2 / 3 * (phi_1 ** 2 + phi_2 ** 2 + phi_3 ** 2))
eps_M_large_strains = np.exp(phi_M) - 1

# small strains
w, v = LA.eig(eps)
eps_1 = w[0]
eps_2 = w[1]
eps_3 = eps_1 + eps_2

eps_M_small_strains = np.sqrt(2 / 3 * (eps_1 ** 2 + eps_2 ** 2 + eps_3 ** 2))

# Wikipedia
# https://en.wikipedia.org/wiki/Infinitesimal_strain_theory#Equivalent_strain
nu = 0.5  # !!!
eps_z = -nu/(1-nu)*(eps_x + eps_y)
eps = np.array([[eps_x, eps_xy, 0], [eps_xy, eps_y, 0], [0, 0, eps_z]])
eps_dev = eps - 1 / 3 * np.trace(eps) * np.eye(3)
eps_vm = np.sqrt(2 / 3 * np.trace(eps_dev @ eps_dev))  # @ is matrix multiplication

print(f"Wikipedia with nu={nu}:       {eps_vm}")
print(f"ARAMIS nodemap export:      {eps_eqv}")
print(f"ARAMIS (large strains):   {eps_M_large_strains}")
print(f"ARAMIS (small strains):  {eps_M_small_strains}")
