"""
Main script to run the phase-field simulation.
Initializes fields, applies boundary conditions, and iterates time steps.
"""
import numpy as np
from config import params
from compute.compute_laplacians import compute_laplacians
from compute.compute_driving import compute_driving
from compute.compute_anisotropy import compute_anisotropy
from compute.update_phi import update_phi
from compute.update_temp import update_temp
from io_utils import write_csv
from plot_utils import plot_field

# Allocate 2D fields

def alloc():
    return np.zeros((params['NX'], params['NY']))

phi      = alloc()
phi_new  = alloc()
temp     = alloc()
temp_new = alloc()

# Initialize phase and temperature fields
for i in range(1, params['NX']-1):
    for j in range(1, params['NY']-1):
        phi[i,j]  = params['PHI_FILL'] if i < params['NX']//4 else 1.0 - params['PHI_FILL']
        temp[i,j] = params['T_FILL']

# Zero-flux boundary conditions
phi[[0,-1],:] = phi[[1,-2],:]
phi[:,[0,-1]] = phi[:,[1,-2]]
temp[[0,-1],:] = temp[[1,-2],:]
temp[:,[0,-1]] = temp[:,[1,-2]]

# Initial output and plots
write_csv('phi_0.csv', phi)
write_csv('temp_0.csv', temp)
plot_field(phi, 'Phase: t=0')
plot_field(temp,'Temp: t=0')

# Time-stepping loop with interrupt handling
try:
    for t in range(1, params['TOTAL_TIMESTEPS']+1):
        for i in range(1, params['NX']-1):
            for j in range(1, params['NY']-1):
                lapX_phi, lapY_phi = compute_laplacians(phi, i, j)
                lapX_tmp, lapY_tmp = compute_laplacians(temp, i, j)
                drv = compute_driving(phi[i,j], temp[i,j])
                eta, e1, e2, th, dxy = compute_anisotropy(phi, i, j)
                phi_new[i,j] = update_phi(
                    phi[i,j], lapX_phi, lapY_phi, drv, eta, e1, e2, th, dxy
                )
                temp_new[i,j] = update_temp(
                    temp[i,j], lapX_tmp, lapY_tmp, phi_new[i,j], phi[i,j]
                )
        # Reapply zero-flux BCs to updated fields
        phi[[0,-1],:] = phi_new[[1,-2],:]
        phi[:,[0,-1]] = phi_new[:,[1,-2]]
        temp[[0,-1],:] = temp_new[[1,-2],:]
        temp[:,[0,-1]] = temp_new[:,[1,-2]]
        # Swap buffers
        phi, phi_new = phi_new, phi
        temp, temp_new = temp_new, temp

        # Periodic output and plotting
        if t % params['TIME_BREAK'] == 0:
            write_csv(f'phi_{t}.csv', phi)
            write_csv(f'temp_{t}.csv', temp)
            print(f"Output at t={t}")
            plot_field(phi, f'Phase: t={t}')
            plot_field(temp, f'Temp: t={t}')
except KeyboardInterrupt:
    print(f"Simulation interrupted at timestep {t}.")
    write_csv(f'phi_int_{t}.csv', phi)
    write_csv(f'temp_int_{t}.csv', temp)
    print("Saved interrupt state.")
