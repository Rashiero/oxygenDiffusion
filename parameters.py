import numpy as np
######################
##### Parameters #####
######################

###### Box #####
# Nx = 100 # Box size => dx = x scale/Nx
# Ny = 100 # Box size => dy = y scale/Ny
# Nz = 100 # Box size => dz = z scale/Nz

# Nx = 133 #P8_Segmented
# Nx = 190 #P9_Segmented
# Nx = 217 #P10_Segmented
# Nx = 250 #P11_Segmented
# Nx = 257 #P12_Segmented
Ny = 512
Nz = 512

margin_x = 40
margin_y = 20
margin_z = 20

'''
<PVStateValue key="micronsPerPixel">
<IndexedValue index="XAxis" value="0.768008424418941"/>
<IndexedValue index="YAxis" value="0.768008424418941"/>
<IndexedValue index="ZAxis" value="1"/>
'''

d_x = 1e-4
d_y = 0.768008424418941e-4
d_z = 0.768008424418941e-4

Lx = d_x*Nx
Ly = d_y*Ny
Lz = d_z*Nz
##### Integration Steps #####
# Lx = 50e-4 #cm
# Ly = 50e-4 #cm
# Lz = 50e-4 #cm

# d_x = Lx/Nx # Lx
# d_y = Ly/Ny # Ly
# d_z = Lz/Nz # Lz

n_print = 1000

##### Constants #####
K = np.longfloat(6e9)
K_goal = np.longfloat(3e8)
#K = np.longfloat(3e12)
# K = 3e9 # mmHg cm s/ cm**2 O2 2.58 (3e8-9.3e7)
speed_const = 1/K # cm**3 O2 / mmHg um s
oxygen_degradation = 1e-2/60 # cm**3 O2 / cm**3 s mmHg <- (0.5-2.5)
D = 6e-10 # 6e-14 cm**3 O2 / um s mmHg
Co = .5 # cm**3 O2 / cm**3 mmHg
n = 3
P50 = 38. # mmHg
alfa = 3.1e-5 # cm**3 O2 / cm**3 mmHg
alfafa = 0.25
P0 = 100. # mmHg

##### Parameters ##### 
Hd = .33
Q =  1e-8/3 # cm**3/s

# Apparent viscosity
## mus 3.2mPa s, 2.4 mPa s
# mu = 3.2*7.50061683e-6
mu = 2.4*7.50061683e-6
# 1Pa = 0.00750061683 mmHg

# Pressures
PB0 = 120
PB1 = 50
