#Called by scrolled_graphene_general for theta = 0 or 30
#Updated to remove edge atoms with only a single carbon nearby. 
#Tearing added, fixed atoms with x < 2*a to aid with tearing
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
#Feb 24 Added substrate dimensions
from folded_general import s_lx
from folded_general import s_ly
from folded_general import s_lz
from folded_general import y_extra
from folded_general import h2o
from folded_general import Hatoms

import math
a=1.4
pi = math.pi
CHbond=1.1

ly_full = ly + 2.0*y_extra

nx = round(lx/(3*a));
ny = round(ly_full/(math.sqrt(3)*a));
ny_fold = round(ly/(math.sqrt(3)*a))
ny_extra = math.floor((ny - ny_fold)/2)
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
    ly_full = ny*math.sqrt(3)*a
else:
    A = math.sqrt(3)*a;
    B = 3*a;
    base = np.array([[0.0, 0.0, 0.0], [0, 2*a, 0.0], [A/2, a/2, 0.0], [A/2, B/2,0.0]])
    nx = round(nx*math.sqrt(3))
    ny = round(ny/math.sqrt(3))
    ny_fold = round(ny_fold/math.sqrt(3))
    ny_extra = round(ny_extra/math.sqrt(3))
    ly_full = ny*3.0*a

# Total number of atoms
N = max(base.shape)*nx*ny;
#Number of atoms if ny=1
Nx = max(base.shape)*nx;
# Calculate the coordinates of the atoms in the layer
coords = np.zeros((Nx,3));
spiralx_coords = np.zeros((Nx,3));
phi = 0
L = r2 + (ftheta*pi/180)*r1
x0 = L+r2*(1-math.cos(ftheta*pi/180))
z0 = r1*(1-math.cos(ftheta*pi/180))+r2*math.sin(ftheta*pi/180)
id = -1;
for ix in range(0,nx):
	for iatom in range(0,max(base.shape)):
		id = id + 1;
		coords[id,:] = base[iatom,:]+[ix*A,0,0];
		if (coords[id,0] < r2): #atom is in the top layer
			#spiralx_coords[id,0] = r2 + L - coords[id,0]
			#spiralx_coords[id,2] = 2.0*r1
			spiralx_coords[id,0] = x0 - coords[id,0]*math.cos(ftheta*pi/180)
			spiralx_coords[id,2] = z0 - coords[id,2]*math.sin(ftheta*pi/180)
		elif (coords[id,0] >= L): #atom is on the bottom layer
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

for iy in range(0,ny_extra):
	spiral_coords[Nx*iy:Nx*(iy+1),:] = coords[:,:] +[0,iy*B,0];

for iy in range(ny_extra,ny_extra+ny_fold):
	spiral_coords[Nx*iy:Nx*(iy+1),:] = spiralx_coords[:,:]+[0,iy*B,0];

for iy in range(ny_extra+ny_fold, ny):
	spiral_coords[Nx*iy:Nx*(iy+1),:] = coords[:,:] +[0,iy*B,0];



#Check for edge atoms with co-ordination number less than 2 and remove them

#April 2024, generalised to account for tearing
#Edge atoms will be along lines y=0, ny_extra, ny_extra+ny_fold and ny. For ny_extra=0 it will just be y=0 and y=ny.
#For these symmetric theta values, we can ignore the "short edges" as these should by definition have the correct CN number except at the corners, which the above should cover.
#All other atoms will be "safe"

tmp_coords = spiral_coords
CNoffset=0
CN2count=0
for atom in range(N):
    if ((periodic == 0) and ((spiral_coords[atom,1] < a) or (spiral_coords[atom,1] > (ly-a)) )):
        radius2 = (spiral_coords[atom,0] - spiral_coords[:,0])**2 + (spiral_coords[atom,1] - spiral_coords[:,1])**2+ (spiral_coords[atom,2] - spiral_coords[:,2])**2
        CN_atoms = np.where(radius2 < (a*1.1)**2)[0]
        CN = len(CN_atoms)
        if (CN < 3):
            #Has a co-ordination number less than 2
            print("removed atom at coords", CN, spiral_coords[atom,0],spiral_coords[atom,1], spiral_coords[atom,2])
            tmp_coords = np.delete(tmp_coords, atom - CNoffset,axis=0)
            CNoffset += 1
            #remove atom
        elif (CN == 3):
            CN2count += 1
    elif (ny_extra > 0):
        if ((abs(spiral_coords[atom,1] - y_extra) < a) or (abs(spiral_coords[atom,1] - (ly+y_extra)) < a)):
            radius2 = (spiral_coords[atom,0] - spiral_coords[:,0])**2+(spiral_coords[atom,1] - spiral_coords[:,1])**2+(spiral_coords[atom,2] - spiral_coords[:,2])**2
            CN_atoms = np.where(radius2 < (a*1.1)**2)[0]
            CN = len(CN_atoms)
            if (CN < 3): #Has a co-ordination number less than 2
                print("removed atom at coords", CN, spiral_coords[atom,0],spiral_coords[atom,1], spiral_coords[atom,2])
                tmp_coords = np.delete(tmp_coords, atom - CNoffset,axis=0)
                CNoffset += 1
	        #remove atom
            elif (CN == 3):
                CN2count+=1
if (CNoffset > 0):
	spiral_coords = tmp_coords
print(CNoffset,'atoms removed')

#Repeat the above just in case

tmp_coords = spiral_coords
CNoffset2=0

HatomCoords = np.zeros((2*CN2count,3))
Hid = 0
for atom in range(N-CNoffset):
    if ((periodic == 0) and ((spiral_coords[atom,1] < a) or (spiral_coords[atom,1] > (ly-a)) )):
        radius2 = (spiral_coords[atom,0] - spiral_coords[:,0])**2 + (spiral_coords[atom,1] - spiral_coords[:,1])**2+ (spiral_coords[atom,2] - spiral_coords[:,2])**2
        CN_atoms = np.where(radius2 < (a*1.1)**2)[0]
        CN = len(CN_atoms)
        if (CN < 3):
            #Has a co-ordination number less than 2
            print("removed atom at coords", CN, spiral_coords[atom,0],spiral_coords[atom,1], spiral_coords[atom,2])
            tmp_coords = np.delete(tmp_coords, atom - CNoffset2,axis=0)
            CNoffset2 += 1
	    #remove atom
        elif (CN == 3):
            if (Hatoms == 1): #Add H atoms to CN2 edge carbon atoms
                HatomCoords[Hid,0] = spiral_coords[atom,0]
                if (spiral_coords[atom,1] <a): #Near bottom
                    HatomCoords[Hid,1] = spiral_coords[atom,1] -CHbond
                else:
                    HatomCoords[Hid,1] = spiral_coords[atom,1] +CHbond
                HatomCoords[Hid,2] = spiral_coords[atom,2]
                Hid +=1
    elif (ny_extra > 0):
        if ((abs(spiral_coords[atom,1] - y_extra) < a) or (abs(spiral_coords[atom,1] - (ly+y_extra)) < a)):
            radius2 = (spiral_coords[atom,0] - spiral_coords[:,0])**2 + (spiral_coords[atom,1] - spiral_coords[:,1])**2+ (spiral_coords[atom,2] - spiral_coords[:,2])**2
            CN_atoms = np.where(radius2 < (a*1.1)**2)[0]
            CN = len(CN_atoms)
            if (CN < 3): #Has a co-ordination number less than 2
                print("removed atom at coords", CN, spiral_coords[atom,0],spiral_coords[atom,1], spiral_coords[atom,2])
                tmp_coords = np.delete(tmp_coords, atom - CNoffset2,axis=0)
                CNoffset2 += 1
				#remove atom
if (CNoffset2 >0):
	spiral_coords = tmp_coords
print(CNoffset2, 'atoms removed')
N = N- CNoffset-CNoffset2

#Add H atoms to short edge on top/bottom layer of fold
NHatoms = 0
if (Hatoms == 1):
    #Bottom layer edge atoms
    xmax = np.amax(spiral_coords[:,0])
    xmaxID = np.where(spiral_coords[:,0] == xmax)[0]
    #Top layer edge atoms
    zmax = 2*r1
    zmaxID = np.where((spiral_coords[:,2] == zmax))[0]
    xmaxzmax = np.amax(spiral_coords[zmaxID,0])
    zmaxedge = np.where((spiral_coords[:,2] == zmax) & (spiral_coords[:,0] == xmaxzmax))[0]
    NHatoms = len(zmaxedge) + len(xmaxID)
    HedgeID = np.concatenate((zmaxedge, xmaxID),axis=0)
    Hatomxyz = np.zeros((NHatoms,3))
    Hatomxyz[:,0] = spiral_coords[HedgeID,0] +CHbond
    Hatomxyz[:,1] = spiral_coords[HedgeID,1]
    Hatomxyz[:,2] = spiral_coords[HedgeID,2]


if (Hid > 0):
    NHatoms += Hid-1
    Hatomxyz = np.concatenate((Hatomxyz,HatomCoords),axis=0)


#Substrate info----------------------------------------------------------------------------------------------------------------------
#Lattice parameter of Si cell
l = 5.43 
#If using a periodic system with Si, then the lattice parameter needs to be rescaled to make the alignments work. For a sufficiently wide box the changes should be minimal
if (periodic == 1):
    nSi = round(ly_full/l)
    print("Old Si lattice parameter", l)
    l = ly_full/nSi
    print("New Si lattice parameter", l, "error percentage", (5.43-l)/5.43*100,"%")
    s_ly = ly_full

#Interatomic distance
a_Si = l*math.sqrt(3)/4
#Spacing between Si and C
#CSi = 3.5 #Approximately the graphene interlayer spacing
CSi = 2.0
x_offset = -10 #Shift substrate off the x=0 line
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
	s_coords[:,0] = s_coords[:,0] + x_offset
	s_coords[:,1] = s_coords[:,1] +(ly_full-s_ly)/2
	s_coords[:,2] = s_coords[:,2] -(s_lz+CSi)

#End Substrate info ------------------------------------------------------------------------------------------------------------------------

#Box dimensions
scale = 1.1
if (periodic == 1):
	scale = 1.0

if ((s_lx > lx) and (Ns !=0)):
	length = s_lx*1.1
else:
	length = 1.1*lx
if ((s_ly > ly_full) and (Ns !=0)):
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


