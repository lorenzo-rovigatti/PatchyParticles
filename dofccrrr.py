#!/usr/bin/env python

import sys
import math as m
import numpy as np
import scipy.linalg as linalg
import scipy

try:
    target_rho = float(sys.argv[1])
    nx, ny, nz = [int (x) for x in sys.argv[2:5]]
    real_N = int(sys.argv[5])
except:
    print >> sys.stderr, "# BCC, 2 particles per cell"
    print >> sys.stderr, "# supply rho, nx, ny, nz, N"
    print >> sys.stderr, "# good guess for rho: 1.0"
    print >> sys.stderr, "# let me suggest a command line: >", sys.argv[0], "1.1 14 6 6 500"
    sys.exit (-2)

cell_volume = 4. / target_rho
a = m.pow(cell_volume, 1./3.)

c1 = np.array([a, 0, 0])
c2 = np.array([0, a, 0])
c3 = np.array([0, 0, a])

# unit cell
rpmol = []
rpmol.append (np.array ([0., 0., 0.]))
rpmol.append (np.array ([0.5, 0.5, 0.]))
rpmol.append (np.array ([0.5, 0., 0.5]))
rpmol.append (np.array ([0., 0.5, 0.5]))

rmol = []
for l in range (nx):
	for m in range (ny):
		for n in range (nz):
			disp = l * c1 + m * c2 + n * c3
			for r in rpmol:
				mine = r[0] * c1 + r[1] * c2 + r[2] * c3 
				rmol.append (disp + mine)

Lx = nx * a
Ly = ny * a
Lz = nz * a

toremove = len(rmol) - real_N

## patch initialisation
f = 1.0 / np.sqrt(3.)
base_patches = []
base_patches.append(np.array([-f, -f, +f]))
base_patches.append(np.array([+f, -f, -f]))
base_patches.append(np.array([+f, +f, +f]))
base_patches.append(np.array([-f, +f, -f]))

Rs = []
al = np.pi / 4.
ct = np.cos(al)
st = np.sin(al)

# rotation matrix of alpha around z axis
Rs.append ([
    np.array([ ct, -st, 0.]),
    np.array([ st,  ct, 0.]),
    np.array([ 0.,  0., 1.])
    ])
Rs.append ([
    np.array([ ct, -st, 0.]),
    np.array([ st,  ct, 0.]),
    np.array([ 0.,  0., 1.])
    ])

al = -np.pi / 4.
ct = np.cos(al)
st = np.sin(al)
Rs.append ([
    np.array([ ct, -st, 0.]),
    np.array([ st,  ct, 0.]),
    np.array([ 0.,  0., 1.])
    ])
Rs.append ([
    np.array([ ct, -st, 0.]),
    np.array([ st,  ct, 0.]),
    np.array([ 0.,  0., 1.])
    ])

print >> sys.stderr, "(stderr) Generating FCC with N =", len(rmol) - toremove, "V =", Lx*Ly*Lz, "Box=", Lx, Ly, Lz
print >> sys.stderr, "(stderr) density of the crystal part", target_rho
print >> sys.stderr, "(stderr) distance between bonded neighs: ", a / np.sqrt(2.) 
print >> sys.stderr, "(stderr) Removing", toremove, "particles from a crystal with", len(rmol)
print >> sys.stderr, "(stderr) Producing file fcc.rrr"

out = open ('fcc.rrr', 'w')
print >> out, "0", len(rmol) - toremove, Lx, Ly, Lz
for i, r in enumerate(rmol):
    if i < toremove:
        continue
    print >> out, Rs[i%4][0][0], Rs[i%4][0][1], Rs[i%4][0][2]
    print >> out, Rs[i%4][1][0], Rs[i%4][1][1], Rs[i%4][1][2]
    print >> out, r[0] + 0.01, r[1] + 0.01 , r[2] + 0.01
out.close()

