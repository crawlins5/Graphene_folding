#Called by folded_graphene_general for theta not = 0 or 30
#Updated to remove atoms on the edge with a CN of 1
import math
import numpy as np
from folded_general import theta
from folded_general import lx
from folded_general import ly
from folded_general import ftheta
from folded_general import r1
from folded_general import r2
from folded_general import filename
from folded_general import writeresults
from folded_general import periodic
#Feb 27 Added substrate dimensions
from folded_general import s_lx
from folded_general import s_ly
from folded_general import s_lz
from folded_general import y_extra
from folded_general import h2o
from folded_general import Hatoms

import math
a=1.4
CHbond = 1.1
theta = theta*math.pi/180;
pi = math.pi
ftheta = ftheta*pi/180
ly_full = ly + 2.0*y_extra

nx = round(lx/(3*a));
ny = round(ly_full/(math.sqrt(3)*a));
ny_fold = round(ly/(math.sqrt(3)*a))
ny_extra = math.floor((ny - ny_fold)/2)

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
    ly_full = ny*B

#End of periodic, may need to consider resizing plan below before continuing
		
#increase number of atoms to account for rotation
#dy = nx*A*math.sin(theta)
#dx = math.sqrt((nx*A)**2+(ny*B)**2)*math.sin(theta+math.atan(nx*A/(ny*B)))

lxf = lx/math.cos(theta) + math.sin(theta)*ly_full
lyf = ly_full + lx*math.tan(theta)

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
	if ((rot_coords[i,0] >= 0.0) and (rot_coords[i,0] < lx) and (rot_coords[i,1] >= 0.0) and (rot_coords[i,1] < ly_full)):
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
    edge_atom = np.where((abs(x_rot-max_x) < 3*a) | (abs(x_rot-min_x) < 3*a) | (abs(y_rot-max_y) < 3*a) | (abs(y_rot-min_y) < 3*a))[0]
else:
    edge_atom = np.where(((abs(x_rot-max_x) < 3*a) | (abs(x_rot-min_x) < 3*a)) & ((abs(y_rot-max_y) > a) &  (abs(y_rot-min_y) > a) ))[0]


j=0
new_rot = rot_coords
print(len(rot_coords),len(edge_atom))
CN2count = 0
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
    elif (CN == 3):
        CN2count += 1

if (j > 0):
    rot_coords = new_rot
print(j,'atoms removed')

HatomID = np.zeros(CN2count,dtype=int)
NHatoms = 0
Hid = 0
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
    edge_atom = np.where((abs(x_rot-max_x) < 3*a) | (abs(x_rot-min_x) < 3*a) | (abs(y_rot-max_y) < 3*a) | (abs(y_rot-min_y) < 3*a))[0]
else:
    edge_atom = np.where(((abs(x_rot-max_x) < 3*a) | (abs(x_rot-min_x) < 3*a)) & ((abs(y_rot-max_y) > a) &  (abs(y_rot-min_y) > a) ))[0]


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
    elif ((CN == 3) and (Hatoms == 1)):
        #Label CN2 carbon atom to have an additional H atoms
        if (NHatoms >= CN2count):
            HatomID = np.append(HatomID, ed-jj)
        else:
            HatomID[Hid] = ed-jj
        Hid += 1
        NHatoms +=1

if (jj > 0):
    rot_coords = new_rot
print(jj,'atoms removed')


#print(len(rot_coords))
N =N-j-jj
#L = r2 +pi*r1
L = r2 +ftheta*r1
spiral_coords = np.zeros((N,3));
x0 = L-r2*math.cos(ftheta)
z0 = r1*(1-math.cos(ftheta))+r2*math.sin(ftheta)
for i in range(0,N):
	if ((rot_coords[i,0] >= L) or (rot_coords[i,1] < y_extra) or (rot_coords[i,1] > (y_extra+ly))): #atom is on the bottom layer or part of the y_extra unfolded segments for tearing
		spiral_coords[i,0] = rot_coords[i,0]
		spiral_coords[i,2] = rot_coords[i,2]
	elif (rot_coords[i,0] < r2): #atom is in the top layer
		#spiral_coords[i,0] = r2 + L - rot_coords[i,0]
		#spiral_coords[i,2] = 2.0*r1
		spiral_coords[i,0] = x0 + rot_coords[i,0]*math.cos(ftheta)
		spiral_coords[i,2] = z0 - rot_coords[i,2]*math.sin(ftheta)
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
CNoffset = 0
if (y_extra != 0): # remove edge atoms at tears
	tmp_coords = spiral_coords
	for atom in range(N):
		if ( (abs(spiral_coords[atom,1] - y_extra) < a) or (abs(spiral_coords[atom,1] - (ly+y_extra)) < a)):
			radius2 = (spiral_coords[atom,0] - spiral_coords[:,0])**2 + (spiral_coords[atom,1] - spiral_coords[:,1])**2+ (spiral_coords[atom,2] - spiral_coords[:,2])**2
			CN_atoms = np.where(radius2 < (a*1.1)**2)[0]
			CN = len(CN_atoms)
			if (CN < 3): #Has a co-ordination number less than 2
				print("removed atom at coords", CN, spiral_coords[atom,0],spiral_coords[atom,1], spiral_coords[atom,2])
				tmp_coords = np.delete(tmp_coords, atom - CNoffset,axis=0)
				CNoffset += 1
				#remove atom

	if (CNoffset > 0):
		spiral_coords = tmp_coords
	print(CNoffset,'atoms removed')

N=N-CNoffset


if (Hatoms == 1):
    Hatomxyz = np.zeros((NHatoms,3))
    Hatomxyz[:,0] = np.where(abs(spiral_coords[HatomID,0]-max_x) < 3*a, spiral_coords[HatomID,0]+CHbond, spiral_coords[HatomID,0]) #3b
    Hatomxyz[:,0] = np.where((abs(rot_coords[HatomID,0]-min_x) < 3*a) & (spiral_coords[HatomID,2] == 2*r1), spiral_coords[HatomID,0]+CHbond, Hatomxyz[:,0]) #3a
    Hatomxyz[:,0] = np.where((abs(spiral_coords[HatomID,0]-min_x) < 3*a), spiral_coords[HatomID,0]-CHbond, Hatomxyz[:,0]) #4
    if (periodic == 0):
        Hatomxyz[:,1] = np.where(abs(spiral_coords[HatomID,1]-max_y) < 3*a, spiral_coords[HatomID,1]+CHbond, spiral_coords[HatomID,1]) #2
        Hatomxyz[:,1] = np.where(abs(spiral_coords[HatomID,1]-min_y) < 3*a, spiral_coords[HatomID,1]-CHbond, Hatomxyz[:,1]) #1
    else:
        Hatomxyz[:,1] = spiral_coords[HatomID,1]

    Hatomxyz[:,2] = spiral_coords[HatomID,2]

#Split into 4 categories:
#1) Atoms along long edge at the y=0 line (-CHbond to y co-ord)
#2) Atoms along long edge at the y=ly line (+CHbond to y co-ord)
#3) Atoms along the short edges (both top 3a and bottom 3b fold) (+CHbond to x co-ord)
#4) If yextra !=0 , the 2 short edges from the flat part of the sheet (-CHbond to x co-ord)
#For periodic set-ups, only #3 and #4 are in the lists
#We currently don't add H atoms along the lines of tearing.


#Substrate info----------------------------------------------------------------------------------------------------------------------
#Lattice parameter of Si cell
l = 5.43 
#Interatomic distance
#If using a periodic system with Si, then the lattice parameter needs to be rescaled to make the alignments work. For a sufficiently wide box the changes should be minimal
if (periodic == 1):
    nSi = round(ly_full/l)
    print("Old Si lattice parameter", l)
    l = ly_full/nSi
    print("New Si lattice parameter", l, "error percentage", (5.43-l)/5.43*100,"%")
    s_ly = ly_full


a_Si = l*math.sqrt(3)/4
#Spacing between Si and C
CSi = 2.0 #Approximately the graphene interlayer spacing
x_offset = -10

mx = round(s_lx/l);
my = round(s_ly/l);
mz = round(s_lz/l);
#Converts dimensions of the sheet to number of unit cells in each direction
print("Old Si dimensions x=",s_lx," y=",s_ly,"z=",s_lz)
s_lx = l*mx
s_ly = l*my
s_lz = l*mz
#Calculates real box size
print("New Si dimensions x=",s_lx," y=",s_ly,"z=",s_lz)


# Coordinates of the 8 atoms in the unit cell
s_base = np.array([[0.0, 0.0, 0.0], [l/2, l/2, 0.0], [0, l/2, l/2], [l/2,0,l/2], [a_Si/math.sqrt(3), a_Si/math.sqrt(3), a_Si/math.sqrt(3)], [l-a_Si/math.sqrt(3), l-a_Si/math.sqrt(3), a_Si/math.sqrt(3)], [a_Si/math.sqrt(3), l-a_Si/math.sqrt(3), l-a_Si/math.sqrt(3)],[l-a_Si/math.sqrt(3), a_Si/math.sqrt(3), l-a_Si/math.sqrt(3)]])

# Total number of Si atoms
Ns = max(s_base.shape)*mx*my*mz;
if (mz==0):
	Ns1 = 0
else:
	Ns1 = max(s_base.shape)*mx*my

#Number of Si atoms in a single layer
# Calculate the coordinates of the atoms in the layer
if (Ns !=0):
	s_coords = np.zeros((Ns,3));
	id = -1;
	for iz in range(0,mz):
		for iy in range(0,my):
			for ix in range(0,mx):
				for iatom in range(0,max(s_base.shape)):
					id = id + 1;
					s_coords[id,:] = s_base[iatom,:]+[ix*l,iy*l,iz*l];
	#Shift substrate to be in-line with graphene
	s_coords[:,0] = s_coords[:,0] +x_offset
	s_coords[:,1] = s_coords[:,1] +(ly_full-s_ly)/2
	s_coords[:,2] = s_coords[:,2] -(s_lz+CSi)

#End Substrate info ------------------------------------------------------------------------------------------------------------------------

#Box dimensions
scale = 1.1
if (periodic == 1):
	scale=1.0

if ((s_lx > lx) and (Ns !=0)):
	length = s_lx*1.1
else:
	length = 1.1*lx
if ((s_ly > ly) and (Ns !=0)):
	width0 = (ly_full-s_ly*scale)/2
	width = (ly_full+s_ly*scale)/2
else:
	width0 = -100.0*(scale-1.0)
	width = ly_full*scale


if writeresults:
    fid = open(filename,'w')
    fid.write('graphene a=%g and Si a=%g \n' % (a,a_Si))
    fid.write('%d atoms\n\n' % (N+Ns+NHatoms))
    natomtypes = 1
    if (Ns != 0):
        natomtypes +=2
    if (y_extra !=0):
        natomtypes +=1
    if ((NHatoms) !=0):
        natomtypes +=1
    fid.write('%d atom types\n\n' %(natomtypes))

    fid.write('%g %g xlo xhi\n' % (x_offset, length))
    if (periodic == 0):
        fid.write('%g %g ylo yhi\n' % (width0, width))
    else:
        fid.write('%g %g ylo yhi\n' % (0.0, B*ny))
    fid.write('%g %g zlo zhi\n\n' % ((-s_lz-CSi)*1.1, 4*r1+1.0))
    fid.write('Masses\n\n' % ())
    fid.write('1 12.0107\n' % ())#free carbon atoms
    if (y_extra > 0):
        fid.write('2 12.0107\n' % ()) #fixed carbon atoms
        if (Ns >0):
            Si_free = 3
            Si_fixed = 4
            Hlabel = 5
        else:
            Hlabel = 3
    else:
        if (Ns >0):
            Si_free = 2
            Si_fixed = 3
            Hlabel = 4
        else:
            Hlabel = 2
    if (Ns >0):
        fid.write('%g 28.0855\n' % (Si_free))# Silicon atoms
        fid.write('%g 28.0855\n' % (Si_fixed)) #Fixed Si atoms
    if ((NHatoms) >0):
        fid.write('%g 1.008\n' % (Hlabel))
    fid.write('\nAtoms\n\n' % ())
    for i in range(0,N):
        if ((spiral_coords[i,1] < (y_extra/2)) or (spiral_coords[i,1] > (ly+(3*y_extra/2)))):# Fixed atoms on y_extra strips
            if (h2o == 1):
                fid.write('%d 2 1 0 %g %g %g\n' % (i+1, spiral_coords[i,0], spiral_coords[i,1], spiral_coords[i,2]))
            else:
                fid.write('%d 2 0 %g %g %g\n' % (i+1, spiral_coords[i,0], spiral_coords[i,1], spiral_coords[i,2]))
        else:
            if (h2o == 1):
                fid.write('%d 1 1 0 %g %g %g\n' % (i+1, spiral_coords[i,0], spiral_coords[i,1], spiral_coords[i,2]))
            else:
                fid.write('%d 1 0 %g %g %g\n' % (i+1, spiral_coords[i,0], spiral_coords[i,1], spiral_coords[i,2]))
    for j in range(0,Ns1):
        fid.write('%d %g 0 %g %g %g\n' % (N+j+1, Si_fixed, s_coords[j,0], s_coords[j,1], s_coords[j,2]))
    for k in range(Ns1,Ns):
        fid.write('%d %g 0 %g %g %g\n' % (N+k+1, Si_free, s_coords[k,0], s_coords[k,1], s_coords[k,2]))
    for hh in range(0,NHatoms):
        fid.write('%d %g 0 %g %g %g\n' % (N+Ns+hh+1, Hlabel, Hatomxyz[hh,0], Hatomxyz[hh,1], Hatomxyz[hh,2]))
    fid.close()


