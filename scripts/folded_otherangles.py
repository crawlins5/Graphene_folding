#Called by folded_graphene_general for theta not = 0 or 30
#Updated to remove atoms on the edge with a CN of 1
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

theta = theta*math.pi/180;
pi = math.pi

nx = round(lx/(3*a));
ny = round(ly/(math.sqrt(3)*a));

# Size of the unit cell
A = 3*a;
B = math.sqrt(3)*a;

# Coordinates of the 4 atoms in the unit cell
base = np.array([[0.0, 0.0, 0.0], [a/2, B/2, 0.0], [A/2, B/2, 0.0], [2*a,0.0,0.0]])

if (periodic == 1):
	width = B*ny*math.cos(theta); #Diagonal width
	i = width*math.sin(theta)*2.0/A
	
	i_low = math.floor(i);
	if ((i_low % 2) == 0):#Check if integer is odd or even
		i_low_even = 1
	else:
		i_low_even = 0
	
	width_low = i_low/(math.sin(theta))*A/2.0
	j_low_raw = width_low*math.cos(theta)*2.0/B
	j_low_up = math.ceil(j_low_raw)
	j_low_down = math.floor(j_low_raw)
	if (i_low_even == 1):
		if ((j_low_up % 2) ==0):
			j_low = j_low_up
		else:
			j_low = j_low_down
	else:
		if ((j_low_up % 2) == 0):
			j_low = j_low_down
		else:
			j_low = j_low_up
	theta_low = math.atan(i_low*math.sqrt(3)/j_low)
	width_low2 = B*j_low/(2.0*math.cos(theta_low))

	i_high = math.ceil(i)
	if (i_low_even == 1):
		i_high_even = 0
	else:
		i_high_even = 1
	
	width_high = i_high/(math.sin(theta))*A/2.0
	j_high_raw = width_high*math.cos(theta)*2.0/B
	j_high_up = math.ceil(j_high_raw)
	j_high_down = math.floor(j_high_raw)
	if (i_high_even == 1):
		if ((j_high_up % 2)==0):
			j_high = j_high_up
		else:
			j_high = j_high_down
	else:
		if ((j_high_up % 2) == 0):
			j_high = j_high_down
		else:
			j_high = j_high_up

	theta_high = math.atan(i_high*math.sqrt(3)/j_high)
	width_high2 = B*j_high/(2.0*math.cos(theta_high))
	dw_low = width - width_low2
	dw_high = width_high2 - width
	print("theta=",theta_low*180/math.pi, "width=", width_low)
	print("theta=", theta_high*180/math.pi, "width=", width_high)
	if (dw_high > dw_low):
		theta = theta_low
		ny = j_low/2.0
	else:
		theta = theta_high
		ny = j_high/2.0
	print("theta changed to", theta*180/math.pi, "and ny changed to", ny)

#End of periodic, may need to consider resizing plan below before continuing
		
#increase number of atoms to account for rotation
#dy = nx*A*math.sin(theta)
#dx = math.sqrt((nx*A)**2+(ny*B)**2)*math.sin(theta+math.atan(nx*A/(ny*B)))

lxf = lx/math.cos(theta) + math.sin(theta)*ly
lyf = ly + lx*math.tan(theta)

nx2 = math.ceil(lxf/A)
ny2 = math.ceil(lyf/B)
print(nx2,ny2)
# Total number of atoms
N = max(base.shape)*nx2*ny2;
# Calculate the coordinates of the atoms in the layer
coords = np.zeros((N,3));
id = 0;
for iy in range(0,ny2):
	for ix in range(0,nx2):
		for iatom in range(0,max(base.shape)):
			#id = id + 1;
			coords[id,:] = base[iatom,:]+[ix*A,iy*B,0];
			id += 1
#rotate graphene
rot_coords = np.zeros((N,3))
rot_coords[:,0] = coords[:,0]*math.cos(theta)+(coords[:,1]-B*ny)*math.sin(theta)
rot_coords[:,1] = -coords[:,0]*math.sin(theta)+(coords[:,1]-B*ny)*math.cos(theta)+B*ny

#Remove atoms which don't fit within the original dimensions
#Reset array size
new_coords = np.zeros((N,3))
Ncount = 0 #counter for atoms
for i in range(N):
	#print(rot_coords[i,0], rot_coords[i,1])
	if ((rot_coords[i,0] >= 0.0) and (rot_coords[i,0] < lx) and (rot_coords[i,1] >= 0.0) and (rot_coords[i,1] < ly)):
		#print("yes")
		new_coords[Ncount,:] = rot_coords[i,:]
		Ncount += 1
#print(N,Ncount)
rot_coords = np.zeros((Ncount,3))
rot_coords[:,:] = new_coords[0:Ncount,:]
N=Ncount

#Remove certain edge atoms from the flat sheet.
#First, identify edge atoms
max_x = np.amax(rot_coords[:,0])
min_x = np.amin(rot_coords[:,0])
max_y = np.amax(rot_coords[:,1])
min_y = np.amin(rot_coords[:,1])
#print(max_x,min_x,max_y, min_y)
x_rot = rot_coords[:,0]
y_rot = rot_coords[:,1]
#print(x_rot,y_rot)
if (periodic == 0):
	edge_atom = np.where((abs(x_rot-max_x) < a) | (abs(x_rot-min_x) < a) | (abs(y_rot-max_y) < a) | (abs(y_rot-min_y) < a))[0]
else:
	edge_atom = np.where((abs(x_rot-max_x) < a) | (abs(x_rot-min_x) < a) )[0]

j=0
new_rot = rot_coords
print(len(rot_coords),len(edge_atom))
for i, ed in enumerate(edge_atom):
	dx = rot_coords[:,0] - rot_coords[ed,0]
	dy = rot_coords[:,1] - rot_coords[ed,1]

	radius2 = dx*dx+dy*dy
	CN_check = np.where(radius2 < 1.1*a*a)[0]
	CN = len(CN_check)
	#print(i,ed,CN)
	#Have to check CN=1 and 2 since it the atom is checked against itself
	if ((CN == 1) or (CN == 2)):
		print("removed atom at coords",CN,ed,  rot_coords[ed,0], rot_coords[ed,1])
		new_rot = np.delete(new_rot, ed-j,axis=0)
		j+=1 
if (j > 0):
	rot_coords = new_rot
print(j,'atoms removed')
#Second edge atoms check 
max_x = np.amax(rot_coords[:,0])
min_x = np.amin(rot_coords[:,0])
max_y = np.amax(rot_coords[:,1])
min_y = np.amin(rot_coords[:,1])
#print(max_x,min_x,max_y, min_y)
x_rot = rot_coords[:,0]
y_rot = rot_coords[:,1]
#print(x_rot,y_rot)
if (periodic == 0):
	edge_atom = np.where((abs(x_rot-max_x) < a) | (abs(x_rot-min_x) < a) | (abs(y_rot-max_y) < a) | (abs(y_rot-min_y) < a))[0]
else:
	edge_atom = np.where((abs(x_rot-max_x) < a) | (abs(x_rot-min_x) < a) )[0]

jj=0
new_rot = rot_coords
print(len(rot_coords),len(edge_atom))
for i, ed in enumerate(edge_atom):
	dx = rot_coords[:,0] - rot_coords[ed,0]
	dy = rot_coords[:,1] - rot_coords[ed,1]

	radius2 = dx*dx+dy*dy
	CN_check = np.where(radius2 < 1.1*a*a)[0]
	CN = len(CN_check)
	#print(i,ed,CN)
	#Have to check CN=1 and 2 since it the atom is checked against itself
	if ((CN == 1) or (CN == 2)):
		print("removed atom at coords",CN,ed,  rot_coords[ed,0], rot_coords[ed,1])
		new_rot = np.delete(new_rot, ed-jj,axis=0)
		jj+=1 
if (jj > 0):
	rot_coords = new_rot
print(jj,'atoms removed')


#print(len(rot_coords))
N =N-j-jj
L = r2 +pi*r1
spiral_coords = np.zeros((N,3));

for i in range(0,N):
	if (rot_coords[i,0] < r2): #atom is in the top layer
		spiral_coords[i,0] = r2 + L - rot_coords[i,0]
		spiral_coords[i,2] = 2.0*r1
	elif (rot_coords[i,0] > L): #atom is on the bottom layer
		spiral_coords[i,0] = rot_coords[i,0]
		spiral_coords[i,2] = rot_coords[i,2]
	else:
		dx = abs(rot_coords[i,0]-rot_coords[i-1,0])
		if (rot_coords[i-1,0] < r2):
			dxx = r2 - rot_coords[i-1,0]
			dpsi = math.atan(dxx/r1)
			rr = math.sqrt(dxx*dxx + r1*r1)
			ddphi = math.acos((r1*r1+rr*rr-dx*dx)/(2.0*r1*rr))
			dphi = ddphi -dpsi
			phi = dphi
		elif (rot_coords[i-1,0] > L):
			if (dx > nx2*A/2.0):
				phi = (rot_coords[i,0]-r2)/r1
			else:
				dxx = rot_coords[i-1,0] - L
				dpsi =  math.atan(dxx/r1)
				rr = math.sqrt(dxx*dxx + r1*r1)
				print(r1,rr,dx)
				ddphi = math.acos((r1*r1+rr*rr-dx*dx)/(2.0*r1*rr))
				dphi = ddphi -dpsi
				phi = pi - dphi
		else:
			dphi = math.acos(1.0 - dx*dx/(2.0*r1*r1))
			phi = phi+dphi
		spiral_coords[i,0] = r1*math.cos(phi+pi/2.0) + L
		spiral_coords[i,2] = r1*math.sin(phi+pi/2.0) + r1

spiral_coords[:,1] = rot_coords[:,1]


print(N)
for j in range(N):
	if (((rot_coords[j,1] < (B*yextra)) or (rot_coords[j,1] > (B*(ny2-yextra))))):
		spiral_coords[j,0] = rot_coords[j,0]
		spiral_coords[j,2] = rot_coords[j,2]

if writeresults:
	fid = open(filename,'w')
	fid.write('graphene a=%g, theta=%g ,ny=%g\n' % (a,theta*180/math.pi,ny))
	fid.write('%g atoms\n\n' % (N))
	if (yextra == 0):
		fid.write('1 atom types\n\n' % ())
	else:
		fid.write('2 atom types\n\n' % ()) #Adding in fixed graphene atoms
	fid.write('%g %g xlo xhi\n' % (0.0, 1.1*lx))
	if (periodic == 0):
		fid.write('%g %g ylo yhi\n' % (-10.0, 1.2*ly))
	else:
		fid.write('%g %g ylo yhi\n' % (0.0, ly))
	fid.write('%g %g zlo zhi\n\n' % (-2*r1, 4*r1))
	fid.write('Masses\n\n' % ())
	if (yextra == 0):
		fid.write('1 12.0107\n\n' % ())#free carbon atoms
	else:
		fid.write('2 12.0107\n\n' % ())#fixed carbon atoms 
	fid.write('Atoms\n\n' % ())
	for i in range(0,N): 
		if (rot_coords[i,0] < L and ((rot_coords[i,1] < (B*yextra)) or (rot_coords[i,1] > (B*(ny2-yextra))))):
			fid.write('%g 2 0 %g %g %g\n' % (i+1, spiral_coords[i,0], spiral_coords[i,1], spiral_coords[i,2]))
		else:
			fid.write('%g 1 0 %g %g %g\n' % (i+1, spiral_coords[i,0], spiral_coords[i,1], spiral_coords[i,2]))

	fid.close()

