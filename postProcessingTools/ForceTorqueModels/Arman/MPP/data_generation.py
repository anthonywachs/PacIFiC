import numpy as np

# Input parameters________________________________________________
D = 0.001
mu = 1e-3
rho = 1000

Re = 10
phi= 0.1

u_c = (Re * mu) / (rho * D)
F_conv = 0.5 * rho * (u_c**2) * ((np.pi * D**2) / 4)
F_st = 3 * np.pi * mu * D * u_c
F_c = F_conv
T_c = F_c * D

# Domain length necessary for periodic mirroring
L = 25
n_closest = 1000

# Averaging window
n_avg = 5000
kernel_width = 4

# loading data________________________________________________________________
x_pos = np.loadtxt('insert_position_x.dat')
y_pos = np.loadtxt('insert_position_y.dat')
z_pos = np.loadtxt('insert_position_z.dat')

x_force = np.loadtxt('x_force')
y_force = np.loadtxt('y_force')
z_force = np.loadtxt('z_force')

x_torque = np.loadtxt('x_torque')
y_torque = np.loadtxt('y_torque')
z_torque = np.loadtxt('z_torque')

u_avg, v_avg, w_avg = (np.loadtxt(f'volumeAverage_velocity_width{kernel_width}dp_kernel0_withPorosity_comp{i}.res') for i in [0, 1, 2])

print(f'The drag array shape: \t{np.shape(x_force)}')
# taking the last row; removing the time column__________________________________
x_pos, y_pos, z_pos = x_pos[-1,1:], y_pos[-1,1:], z_pos[-1,1:]

# removing the time column
x_force, y_force, z_force = np.delete(x_force, [0], axis=1),\
                            np.delete(y_force, [0], axis=1),\
                            np.delete(z_force, [0], axis=1)
                            
x_torque, y_torque, z_torque = np.delete(x_torque, [0], axis=1),\
                               np.delete(y_torque, [0], axis=1),\
                               np.delete(z_torque, [0], axis=1)

u_avg, v_avg, w_avg =   np.delete(u_avg, [0], axis=1),\
                        np.delete(v_avg, [0], axis=1),\
                        np.delete(w_avg, [0], axis=1)

# reshape to have an array with one column______________________________________
x_pos, y_pos, z_pos =   x_pos.reshape(-1,1),\
                        y_pos.reshape(-1,1),\
                        z_pos.reshape(-1,1)

x_force, y_force, z_force = np.mean(x_force[-n_avg:,:], axis=0),\
                            np.mean(y_force[-n_avg:,:], axis=0),\
                            np.mean(z_force[-n_avg:,:], axis=0)
x_force, y_force, z_force = x_force.reshape(-1,1),\
                            y_force.reshape(-1,1),\
                            z_force.reshape(-1,1)
                            
x_torque, y_torque, z_torque = np.mean(x_torque[-n_avg:,:], axis=0),\
                               np.mean(y_torque[-n_avg:,:], axis=0),\
                               np.mean(z_torque[-n_avg:,:], axis=0)
x_torque, y_torque, z_torque = x_torque.reshape(-1,1),\
                               y_torque.reshape(-1,1),\
                               z_torque.reshape(-1,1)

u_avg, v_avg, w_avg = np.mean(u_avg[-5:,:], axis=0),\
                      np.mean(v_avg[-5:,:], axis=0),\
                      np.mean(w_avg[-5:,:], axis=0)
u_avg, v_avg, w_avg = u_avg.reshape(-1,1),\
                      v_avg.reshape(-1,1),\
                      w_avg.reshape(-1,1)

# non-dimensionalization_____________________________________________________
x_pos, y_pos, z_pos =   x_pos / D,\
                        y_pos / D,\
                        z_pos / D

u_avg, v_avg, w_avg =   u_avg / u_c,\
                        v_avg / u_c,\
                        w_avg / u_c

x_force, y_force, z_force =     x_force / F_c,\
                                y_force / F_c,\
                                z_force / F_c

x_torque, y_torque, z_torque =  x_torque / T_c,\
                                y_torque / T_c,\
                                z_torque / T_c

# printing information _____________________________________________________
F_d_F_st = x_force.mean() * F_c / F_st
print(f'F_d / F_st \t= \t{F_d_F_st:.2f}')
print(f'sigma_d / F_d \t= \t{(x_force.std() / x_force.mean()) * 100:.2f} %')
F_L = np.hstack((y_force, z_force))
F_L_mag = np.linalg.norm(F_L, axis=0)
print(f'F_L / F_st \t= \t{F_L.mean():.2g}')
ratio = F_L.std() / x_force.mean()
print(f'sigma_L / F_d \t= \t{ratio * 100:.2f} %')
                   
# Generate the data matrix with positions_______________________________
xyz = np.hstack((x_pos, y_pos, z_pos))

# The extension of periodic mirror___________________________________
delta = 0.5
D = 1

x_back  = xyz[xyz[:,0] > (L*delta + D), 0].reshape(-1,1) - L
x_front = xyz[xyz[:,0] < (L*delta + D), 0].reshape(-1,1) + L

y_back  = xyz[xyz[:,1] > (L*delta + D), 1].reshape(-1,1) - L
y_front = xyz[xyz[:,1] < (L*delta + D), 1].reshape(-1,1) + L

z_back  = xyz[xyz[:,2] > (L*delta + D), 2].reshape(-1,1) - L
z_front = xyz[xyz[:,2] < (L*delta + D), 2].reshape(-1,1) + L

xyz_x_back    = np.hstack((x_back, y_pos[xyz[:,0] > (L*delta + D)], z_pos[xyz[:,0] > (L*delta + D)]))
xyz_x_front   = np.hstack((x_front, y_pos[xyz[:,0] < (L*delta + D)], z_pos[xyz[:,0] < (L*delta + D)]))

xyz_y_back    = np.hstack((x_pos[xyz[:,1] > (L*delta + D)], y_back, z_pos[xyz[:,1] > (L*delta + D)]))
xyz_y_front   = np.hstack((x_pos[xyz[:,1] < (L*delta + D)], y_front, z_pos[xyz[:,1] < (L*delta + D)]))

xyz_z_back    = np.hstack((x_pos[xyz[:,2] > (L*delta + D)], y_pos[xyz[:,2] > (L*delta + D)], z_back))
xyz_z_front   = np.hstack((x_pos[xyz[:,2] < (L*delta + D)], y_pos[xyz[:,2] < (L*delta + D)], z_front))

num_p = np.size(xyz, axis=0)
NN_data = np.zeros( (num_p, 3 * (n_closest - 1) + 9 ) )

# here, I add the periodic images to the original xyz to be considered in the distances

xyz = np.vstack((xyz, xyz_x_back, xyz_x_front, xyz_y_back, xyz_y_front, xyz_z_back, xyz_z_front))
# the distance loop should be done over periodic particles as well
num_p_periodic = np.size(xyz, axis=0)

for m in np.arange(num_p):

    # Just a progress indicator
    if m%100 == 0 and m!=0:
        print('Particle '+str(m)+' is done!')

    dist  = np.zeros( (num_p_periodic, 4) )

    for i in range(num_p_periodic):
        dist[i,0] = np.sqrt( ( xyz[i, 0] - xyz[m, 0] )**2 +
                             ( xyz[i, 1] - xyz[m, 1] )**2 +
                             ( xyz[i, 2] - xyz[m, 2] )**2 )
        dist[i,1] = xyz[i, 0] - xyz[m, 0]
        dist[i,2] = xyz[i, 1] - xyz[m, 1]
        dist[i,3] = xyz[i, 2] - xyz[m, 2]

    dist = dist[ np.argsort( dist[:, 0] ) ]
    dist = dist[:n_closest, :]
    dist = dist[1:, 1:]
    dist = dist.flatten()
    dist = dist.reshape(1, -1)

    NN_data[m,:dist.size] = dist

    NN_data[m,-9] = u_avg[m, 0]
    NN_data[m,-8] = v_avg[m, 0]
    NN_data[m,-7] = w_avg[m, 0]
    
    NN_data[m,-6] = x_torque[m, 0]
    NN_data[m,-5] = y_torque[m, 0]
    NN_data[m,-4] = z_torque[m, 0]
    
    NN_data[m,-3] = x_force[m, 0]
    NN_data[m,-2] = y_force[m, 0]
    NN_data[m,-1] = z_force[m, 0]

np.savetxt("new_3D_Re10_phi01_M=1000", NN_data)
