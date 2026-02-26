#Called by folded_graphene_general for all theta
import math
import numpy as np
from folded_general import theta
from folded_general import lx
from folded_general import ly
from folded_general import ftheta
from folded_general import r1
from folded_general import r2
from folded_general import filename
from folded_general import filename2

from folded_general import writeresults
from folded_general import periodic_y
#Sept 25 For periodic layers
from folded_general import periodic_x
#Feb 24 Added substrate dimensions
from folded_general import s_lx
from folded_general import s_ly
from folded_general import s_lz
from folded_general import y_extra
from folded_general import h2o
from folded_general import Hatoms
from folded_general import layer_num

import math
a=1.4
CHbond = 1.1
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
if ((theta==0) or (theta==30)):
	ly_full = ny*B
	lx = nx*A



theta = theta*math.pi/180;
pi = math.pi
ftheta = ftheta*pi/180

if ((periodic_y == 1) and (theta!=0.0)):
	width = B*ny#*math.cos(theta); #Diagonal width
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
	print("theta=",theta_low*180/math.pi, "width=", width_low, width_low2)
	print("theta=", theta_high*180/math.pi, "width=", width_high, width_high2)
	if (dw_high > dw_low):
		theta = theta_low
		ly_full = width_low2
		ny = j_low/2.0
	else:
		theta = theta_high
		ly_full = width_high2
		ny = j_high/2.0
	print("theta changed to", theta*180/math.pi, "and width changed to", ly_full)
	if (periodic_x == 1):
		lx_scale = lx/math.sqrt(3.0)
		scaling = lx_scale/ly_full*ny
        #print(scaling)
		lx_new = round(scaling)/ny*ly_full*math.sqrt(3.0)
		print("Length changed to", lx_new)
		lx = lx_new
ly = ly_full - 2.0*y_extra
    #ly_full = ny*B

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
z0=0.0
if (N != 0):
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
	rot_coords[:,0] = coords[:,0]*math.cos(theta)+(coords[:,1]-ly_full)*math.sin(theta)
	rot_coords[:,1] = -coords[:,0]*math.sin(theta)+(coords[:,1]-ly_full)*math.cos(theta)+ly_full

    #Remove atoms which don't fit within the original dimensions
    #Reset array size
	new_coords = np.zeros((N,3))
	Ncount = 0 #counter for atoms
	for i in range(N):
        #print(rot_coords[i,1],ly_full, ly_full-rot_coords[i,1])
            #print(rot_coords[i,0], rot_coords[i,1])
		if ((rot_coords[i,0] >= 0.0) and (rot_coords[i,0] < lx) and (rot_coords[i,1] >= 0.0) and (rot_coords[i,1] < ly_full) and (abs(rot_coords[i,0] - lx) > 0.001) and (abs(rot_coords[i,1] - ly_full) > 0.001)):
                    #print("yes")
			new_coords[Ncount,:] = rot_coords[i,:]
			Ncount += 1
    #print(N,Ncount)
    
    #Remove certain edge atoms from the flat sheet.
	nrm = 0
	j=-1
	CNH = np.zeros(3,dtype=int)
	    #Cycles the check for removal of atoms until the number of removed atoms is zero
	while (j !=0):
		#The below two settings only matter for last while loop
		NHatoms=0
		HatomID = np.zeros((0,3),dtype=int)

	    #First, identify edge atoms
		max_x = np.amax(new_coords[:Ncount,0])
		min_x = np.amin(new_coords[:Ncount,0])
		max_y = np.amax(new_coords[:Ncount,1])
		min_y = np.amin(new_coords[:Ncount,1])
		x_rot = new_coords[:Ncount,0]
		y_rot = new_coords[:Ncount,1]
	    #print(x_rot,y_rot)
		if (periodic_y == 0):
			if (periodic_x == 0):
				edge_atom = np.where((abs(x_rot-max_x) < 3*a) | (abs(x_rot-min_x) < 3*a) | (abs(y_rot-max_y) < 3*a) | (abs(y_rot-min_y) < 3*a))[0]
			else:
				edge_atom = np.where(((x_rot > (min_x+A/2.1)) & (x_rot < (max_x-A/2.1))) & ((abs(y_rot-max_y) < 3*a) |  (abs(y_rot-min_y) < 3*a) ))[0]
		else:
			if (periodic_x == 0):
				edge_atom = np.where(((abs(x_rot-max_x) < 3*a) | (abs(x_rot-min_x) < 3*a)) & ((y_rot >(min_y+B/2.1)) & (y_rot < (max_y-B/2.1) )))[0]
			else:
				edge_atom = np.empty(0)

		j=0
		new_rot = new_coords[:Ncount,:]
	    #print(len(new_coords[:Ncount,:]),len(edge_atom))
		for i, ed in enumerate(edge_atom):
			dx = new_coords[:Ncount,0] - new_coords[ed,0]
			dy = new_coords[:Ncount,1] - new_coords[ed,1]

			radius2 = dx*dx+dy*dy
			CN_check = np.where((radius2 < 1.1*a*a) & (radius2 > 0.1))[0]
			CN = len(CN_check)
		#print(i,ed,CN)
			if ((CN == 0) or (CN == 1)):
				print("removed atom at coords",CN,ed,  new_coords[ed,0], new_coords[ed,1])
				new_rot = np.delete(new_rot, ed-j,axis=0)
				j+=1
			elif (CN == 2):
				CNH[0] = ed
				CNH[1:] = CN_check
				HatomID = np.append(HatomID, [CNH],axis=0)
				NHatoms +=1
	    
		if (j > 0):
			new_coords = new_rot
		nrm +=j
		print(nrm,j,'atoms removed')
		
    
	Ncount = Ncount-nrm
	rot_coords = np.zeros((Ncount*layer_num,3))
	rot_coords[:Ncount,:] = new_coords[0:Ncount,:]#Top layer
	rot_coords[:Ncount,2] = rot_coords[:Ncount,2]+(2*r1)*(layer_num-1) #Shift top layer up
	rot_coords[Ncount:Ncount*layer_num,:] = new_coords[0:Ncount*(layer_num-1),:] #Bottom layer if layer_num=2
	rot_coords[Ncount:Ncount*layer_num,1] = ly_full-rot_coords[Ncount:Ncount*layer_num,1] #Flip bottom layer
	N=Ncount*layer_num

    #print(len(rot_coords))
    #L = r2 +pi*r1
	L = r2 +ftheta*r1
	spiral_coords = np.zeros((N,3));
	x0 = L-r2*math.cos(ftheta)-r1*math.sin(ftheta)
	z0 = r1*(1-math.cos(ftheta))+r2*math.sin(ftheta)
	for i in range(0,N):
		if ((rot_coords[i,0] >= L) or (rot_coords[i,1] < y_extra) or (rot_coords[i,1] > (y_extra+ly))): #atom is on the bottom layer or part of the y_extra unfolded segments for tearing
			spiral_coords[i,0] = rot_coords[i,0]
			spiral_coords[i,2] = rot_coords[i,2]
			if (rot_coords[i,0]< L):
				if (rot_coords[i,0] < r2):
					phi = 0
				else:
					phi = (rot_coords[i,0]-r2)/r1

		elif (rot_coords[i,0] < r2): #atom is in the top layer
			#spiral_coords[i,0] = r2 + L - rot_coords[i,0]
 			#spiral_coords[i,2] = 2.0*r1
			spiral_coords[i,0] = x0 + rot_coords[i,0]*math.cos(ftheta)
			spiral_coords[i,2] = z0 - rot_coords[i,0]*math.sin(ftheta)+(layer_num-1)*2.0*r1
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
					ddphi = math.acos((r1*r1+rr*rr-dx*dx)/(2.0*r1*rr))
					dphi = ddphi -dpsi
					phi = pi - dphi
			else:
                            #print(rot_coords[i-1,0],rot_coords[i-1,1], rot_coords[i-1,2])
                            #print(rot_coords[i,0],rot_coords[i,1], rot_coords[i,2])
				dphi = math.acos(1.0 - dx*dx/(2.0*r1*r1))
				phi = phi+dphi
			spiral_coords[i,0] = L - r1*math.sin(phi+pi-ftheta) 
			spiral_coords[i,2] = r1 + r1*math.cos(phi+pi-ftheta)+(layer_num-1)*2.0*r1

	spiral_coords[:,1] = rot_coords[:,1]
	if (y_extra != 0): # remove edge atoms at tears
		CNoffset = -1
		nrm2 = 0
		while (CNoffset !=0):
			NHatoms2 = 0
			HatomID2 = np.zeros((0,3),dtype=int)
			tmp_coords = spiral_coords
			CNoffset = 0
			for atom in range(N):
				if ( (abs(spiral_coords[atom,1] - y_extra) < 3*a) or (abs(spiral_coords[atom,1] - (ly+y_extra)) < 3*a)):
					radius2 = (spiral_coords[atom,0] - spiral_coords[:,0])**2 + (spiral_coords[atom,1] - spiral_coords[:,1])**2+ (spiral_coords[atom,2] - spiral_coords[:,2])**2
					CN_atoms = np.where((radius2 < (a*a*1.1))& (radius2>0.1))[0]
					CN = len(CN_atoms)
					if ((CN == 0) or (CN == 1)):
						print("removed atom at coords", atom, spiral_coords[atom,0],spiral_coords[atom,1], spiral_coords[atom,2])
						tmp_coords = np.delete(tmp_coords, atom - CNoffset,axis=0)
						if (np.any(HatomID == (atom-nrm2-CNoffset))):
							rmH = np.where(HatomID == (atom-CNoffset-nrm2))[0]
							print(atom, rmH,len(rmH))
							for k in range(len(rmH)):
								print(rmH[k]-k, HatomID[rmH[k]-k])
								HatomID = np.delete(HatomID, rmH[k]-k,axis=0)
							NHatoms -= len(rmH)
						HatomID = np.where(HatomID > atom, HatomID - 1,HatomID)
						#Adjust index of Hatom IDs to compensate for removed atoms
	
					    #remove atom
						CNoffset +=1
					elif ((CN == 2) and (spiral_coords[atom,2]>0.0)):
				    #Label CN2 carbon atom to have an additional H atoms
				    #Only add to carbon atoms on the ribbon, not on the flat portion
						CNH[0] = atom
						CNH[1:] = CN_atoms
						HatomID2 = np.append(HatomID2, [CNH],axis=0)
						NHatoms2 +=1

			nrm2 += CNoffset			    
			if (CNoffset > 0):
				spiral_coords = tmp_coords
			print(nrm2,CNoffset,'atoms removed')
			N = N - CNoffset

	#NHatoms = NHatoms+NHatoms2
	#HatomID = np.append(HatomID,HatomID2, axis=0)
	#print(N,NHatoms, HatomID)
	if (Hatoms == 1):
		Hatomxyz = np.zeros((NHatoms*layer_num,3))
		for i, ed in enumerate(HatomID[:,0]):
			atm1 = HatomID[i,1]
			atm2 = HatomID[i,2]
			atm3x = spiral_coords[ed,0]
			atm3y = spiral_coords[ed,1]
			atm3z = spiral_coords[ed,2]

			xave = (spiral_coords[atm1,0] + spiral_coords[atm2,0])/2.0
			yave = (spiral_coords[atm1,1] + spiral_coords[atm2,1])/2.0
			zave = (spiral_coords[atm1,2] + spiral_coords[atm2,2])/2.0

			dx = atm3x - xave
			dy = atm3y - yave
			dz = atm3z - zave
			Hx = dx*CHbond*2.0/a
			Hy = dy*CHbond*2.0/a
			Hz = dz*CHbond*2.0/a
			Hatomxyz[i,0] = atm3x+Hx
			Hatomxyz[i,1] = atm3y+Hy
			Hatomxyz[i,2] = atm3z+Hz

			if (layer_num == 2):
				Hatomxyz[i+NHatoms,0] = Hatomxyz[i,0]
				Hatomxyz[i+NHatoms,1] = ly_full - Hatomxyz[i,1]
				Hatomxyz[i+NHatoms,2] = Hatomxyz[i,2]-2*r1

            #print(i,ed, HatomID[i,:])
            #print(xave,yave,zave, dx,dy,dz)

#Substrate info----------------------------------------------------------------------------------------------------------------------
#Lattice parameter of Si cell
l = 5.43 
#Interatomic distance
#If using a periodic system with Si, then the lattice parameter needs to be rescaled to make the alignments work. For a sufficiently wide box the changes should be minimal
if (periodic_y == 1):
	if (N == 0):
		nSi = round(s_ly/l)
		print("Old Si lattice parameter", l)
		l = s_ly/nSi
		print("New Si lattice parameter", l, "error percentage", (5.43-l)/5.43*100,"%")
	else:
		nSi = round(ly_full/l)
		print("Old Si lattice parameter", l)
		l = ly_full/nSi
		print("New Si lattice parameter", l, "error percentage", (5.43-l)/5.43*100,"%")
		s_ly = ly_full


a_Si = l*math.sqrt(3)/4
#Spacing between Si and C
CSi = 2.0 #Approximately the graphene interlayer spacing
x_offset = -10
if (periodic_x == 1):
	x_offset = 0.0

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
scale_y = 1.1
scale_x = 1.1
if (periodic_y == 1):
	scale_y=1.0

if (periodic_x == 1):
	scale_x=1.0

if ((s_lx > lx) and (Ns !=0)):
	length = s_lx*scale_x
else:
	length = scale_x*lx

if ((s_ly > ly_full) and (Ns !=0)):
	width0 = (ly_full-s_ly*scale_y)/2
	width = (ly_full+s_ly*scale_y)/2
else:
	width0 = -100.0*(scale_y-1.0)
	width = ly_full*scale_y

if writeresults:
	fid = open(filename,'w')
	fid2 = open(filename2, 'w')
	fid.write('graphene a=%g and Si a=%g \n' % (a,a_Si))
	fid.write('%d atoms\n\n' % (N+Ns+NHatoms*layer_num))
	fid2.write('%d\n' % (N+Ns+NHatoms*layer_num))
	fid2.write('Layers=%d \t length = %g \t width = %g \n' % (layer_num, length, (width-width0)))
	natomtypes = 0
	if (N != 0):
		natomtypes += 1
	if (Ns != 0):
		natomtypes +=2
	if (y_extra !=0):
		natomtypes +=1
	if ((NHatoms) !=0):
		natomtypes +=1
	fid.write('%d atom types\n\n' %(natomtypes))

	fid.write('%g %g xlo xhi\n' % (x_offset, length))
	fid.write('%g %g ylo yhi\n' % (width0, width))
	fid.write('%g %g zlo zhi\n\n' % ((-s_lz-CSi)*1.1, z0+1.0+(layer_num-1)*2.0*r1))
	fid.write('Masses\n\n' % ())
	if (N > 0):
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
	else:#Just substrate
		Si_free = 1
		Si_fixed = 2
	if (Ns >0):
		fid.write('%g 28.0855\n' % (Si_free))# Silicon atoms
		fid.write('%g 28.0855\n' % (Si_fixed)) #Fixed Si atoms
	if ((NHatoms) >0):
		fid.write('%g 1.008\n' % (Hlabel))
	fid.write('\nAtoms\n\n' % ())
	for i in range(0,N):
		fid2.write('C \t %g \t %g \t %g \n' % (spiral_coords[i,0], spiral_coords[i,1], spiral_coords[i,2]))
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
		fid2.write('Si \t %g \t %g \t %g \n' % (s_coords[j,0], s_coords[j,1], s_coords[j,2]))

		fid.write('%d %g 0 %g %g %g\n' % (N+j+1, Si_fixed, s_coords[j,0], s_coords[j,1], s_coords[j,2]))
	for k in range(Ns1,Ns):
		fid2.write('Si \t %g \t %g \t %g \n' % (s_coords[k,0], s_coords[k,1], s_coords[k,2]))
		fid.write('%d %g 0 %g %g %g\n' % (N+k+1, Si_free, s_coords[k,0], s_coords[k,1], s_coords[k,2]))
	for hh in range(0,NHatoms*layer_num):
		fid.write('%d %g 0 %g %g %g\n' % (N+Ns+hh+1, Hlabel, Hatomxyz[hh,0], Hatomxyz[hh,1], Hatomxyz[hh,2]))
		fid2.write('H \t  %g \t %g \t  %g\n' % (Hatomxyz[hh,0], Hatomxyz[hh,1], Hatomxyz[hh,2]))

	fid.close()
	fid2.close()

