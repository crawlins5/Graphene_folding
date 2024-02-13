#Called by scrolled_graphene_general for theta = 0 or 30
#Updated to remove edge atoms with only a single carbon nearby. 
import math
import numpy as np
from folded_general import theta
from folded_general import lx
from folded_general import ly
from folded_general import a
from folded_general import r1
from folded_general import r2
from folded_general import yextra
from folded_general import filename
from folded_general import writeresults
from folded_general import periodic

import math
pi = math.pi

nx = round(lx/(3*a));
ny = round(ly/(math.sqrt(3)*a));
#Converts dimensions of the sheet to number of unit cells in each direction

if (theta == 0):
    zigzagedge = True;
else:
    zigzagedge = False;

if zigzagedge:
    # Size of the unit cell
	A = 3*a;
	B = math.sqrt(3)*a;

    # Coordinates of the 4 atoms in the unit cell
	base = np.array([[0.0, 0.0, 0.0], [a/2, B/2, 0.0], [A/2, B/2, 0.0], [2*a,0.0,0.0]])
else:
	A = math.sqrt(3)*a;
	B = 3*a;
	base = np.array([[0.0, 0.0, 0.0], [0, 2*a, 0.0], [A/2, a/2, 0.0], [A/2, B/2,0.0]])
	nx = round(nx*math.sqrt(3))
	ny = round(ny/math.sqrt(3))

# Total number of atoms
N = max(base.shape)*nx*ny;
#Number of atoms if ny=1
Nx = max(base.shape)*nx;
# Calculate the coordinates of the atoms in the layer
coords = np.zeros((Nx,3));
spiralx_coords = np.zeros((Nx,3));
phi = 0
L = r2 + pi*r1
id = -1;
for ix in range(0,nx):
	for iatom in range(0,max(base.shape)):
		id = id + 1;
		coords[id,:] = base[iatom,:]+[ix*A,0,0];
		if (coords[id,0] < r2): #atom is in the top layer
			spiralx_coords[id,0] = r2 + L - coords[id,0]
			spiralx_coords[id,2] = 2.0*r1
		elif (coords[id,0] > L): #atom is on the bottom layer
			spiralx_coords[id,0] = coords[id,0]
			spiralx_coords[id,2] = coords[id,2]
		else:
			dx = abs(coords[id,0]-coords[id-1,0])
			if (coords[id-1,0] < r2):
				dxx = r2 - coords[id-1,0]
				dpsi = math.atan(dxx/r1)
				rr = math.sqrt(dxx*dxx + r1*r1)
				ddphi = math.acos((r1*r1+rr*rr-dx*dx)/(2.0*r1*rr))
				dphi = ddphi -dpsi
				phi = dphi
			elif (coords[id-1,0] > L):
				dxx = coords[id-1,0] - L
				dpsi =  math.atan(dxx/r1)
				rr = math.sqrt(dxx*dxx + r1*r1)
				ddphi = math.acos((r1*r1+rr*rr-dx*dx)/(2.0*r1*rr))
				dphi = ddphi -dpsi
				phi = pi - dphi
			else:
				dphi = math.acos(1.0 - dx*dx/(2.0*r1*r1))
				phi = phi+dphi
			spiralx_coords[id,0] = r1*math.cos(phi+pi/2.0) + L
			spiralx_coords[id,2] = r1*math.sin(phi+pi/2.0) + r1

spiralx_coords[:,1] = coords[:,1];

spiral_coords = np.zeros((N,3));

for iy in range(0,ny):
	if ((iy > yextra-1) and (iy < (ny-yextra))):
		spiral_coords[Nx*iy:Nx*(iy+1),:] = spiralx_coords[:,:]+[0,iy*B,0];
	else:
		spiral_coords[Nx*iy:Nx*(iy+1),:] = coords[:,:]+[0,iy*B,0];

CNoffset = 2
edge1 = np.zeros((2,3))
#print(spiral_coords[0,:], spiral_coords[Nx*(ny-1)+1,:])
if (theta == 0):
	edge1[0] += spiral_coords[0,:];
	edge1[1] += spiral_coords[Nx-1,:];
	spiral_coords[0:Nx-2,:] = spiral_coords[1:Nx-1,:];
	spiral_coords[Nx-2:N-2,:] = spiral_coords[Nx:N,:]
	spiral_coords[N-2,:] = edge1[0];
	spiral_coords[N-1,:] = edge1[1]

	#Will have to add something for yextra >0
elif (theta == 30):
	edge1[0] += spiral_coords[0,:];
	edge1[1] += spiral_coords[Nx*(ny-1)+1,:]
	spiral_coords[0:Nx*(ny-1),:] = spiral_coords[1:Nx*(ny-1)+1,:]
	spiral_coords[Nx*(ny-1): N-2,:] = spiral_coords[Nx*(ny-1)+2:N,:]
	spiral_coords[N-2,:] = edge1[0]
	spiral_coords[N-1,:] = edge1[1]

if writeresults:
	fid = open(filename,'w')
	fid.write('graphene a=%g\n' % (a))
	fid.write('%g atoms\n\n' % (N-CNoffset))
	if (yextra == 0):
		fid.write('1 atom types\n\n' % ())
	else:
		fid.write('2 atom types\n\n' % ()) #Adding in fixed graphene atoms
	fid.write('%g %g xlo xhi\n' % (0.0, 1.1*A * nx))
	if (periodic == 0):
		fid.write('%g %g ylo yhi\n' % (-10.0, 1.1*B * ny))
	else:
		fid.write('%g %g ylo yhi\n' % (0.0, B*ny))
	fid.write('%g %g zlo zhi\n\n' % (-2*r1, 4*r1))
	fid.write('Masses\n\n' % ())
	if (yextra == 0):
		fid.write('1 12.0107\n\n' % ())#free carbon atoms
	else:
		fid.write('2 12.0107\n\n' % ())#fixed carbon atoms 
	fid.write('Atoms\n\n' % ())
	for i in range(0,N-CNoffset):
		if ((spiral_coords[i,0] < L) and ((spiral_coords[i,1] < float(yextra)*B) or (spiral_coords[i,1] > float(ny-yextra)*B))):#Fixed atoms
			fid.write('%g 2 0 %g %g %g\n' % (i+1, spiral_coords[i,0], spiral_coords[i,1], spiral_coords[i,2]))
		else:#Free atoms
			fid.write('%g 1 0 %g %g %g\n' % (i+1, spiral_coords[i,0], spiral_coords[i,1], spiral_coords[i,2]))

	fid.close()


