#!/usr/bin/env python
# coding: utf-8
#November 2023: Adding fold angle derived from velocity differences
# Changing velocity calculations. Previous idea was based on velocity difference between edges of sheet.
# New idea is to measure velocity difference between edge and closest atoms in other fold

# In[1]:


import math
import numpy as np
from functools import reduce

input = []
with open('INPUT','r') as file:
	input = file.readlines()
for j in input:
	input = j.split(" ")
	if (input[0].startswith('timestep')):
		timestep = float(input[1]) #fs
		timestep = timestep/1000 #ps



contents = []
with open('graphene_xyz_v.lammpstrj','r') as file:
	contents = file.readlines()
line_start = 1
line_num = 0
pi = math.pi
bond_length = 1.4
num_atoms = contents[3]
total_steps = contents.count(num_atoms)

num_atoms = int(num_atoms)
total_steps = int(total_steps)

nt = -1
t_switch = 0 #Switch for reading in timestep value
box_switch = -1 #Switch for reading in box dimensions
#index = numpy.zeros(num_atoms*total_steps, dtype=int)
x = np.zeros((total_steps,num_atoms))
y = np.zeros((total_steps,num_atoms))
z = np.zeros((total_steps,num_atoms))
vx = np.zeros((total_steps,num_atoms))
vy = np.zeros((total_steps,num_atoms))
vz = np.zeros((total_steps,num_atoms))
t = np.zeros(total_steps, dtype=int)

box_x = np.zeros((total_steps,2))
box_y = np.zeros((total_steps,2))
box_z = np.zeros((total_steps,2))

for i in contents:
	coords = i.split(' ')
	if t_switch == 1:
		t[nt] = int(coords[0])
		t_switch = 0
	if (len(coords) >1 ):
		if (coords[1].startswith('TIME')):
			nt += 1
			t_switch = 1
	if (len(coords) == 8):
		x[nt,int(coords[0])-1] = float(coords[2])#*box_dx
		y[nt,int(coords[0])-1] = float(coords[3])#*box_dy
		z[nt,int(coords[0])-1] = float(coords[4])#*box_dz
		vx[nt,int(coords[0])-1] = float(coords[5])
		vy[nt,int(coords[0])-1] = float(coords[6])
		vz[nt,int(coords[0])-1] = float(coords[7])




#Index of atoms initially on the top layer
zmax_value = np.amax(z[0])
zmax_index = np.where(z[0] == zmax_value)[0]

dx_init = np.amax(x[0,zmax_index]) -np.amin(x[0,zmax_index])
#print(dx_init)


#Index of atoms initially on the bottom layer
zmin_value = np.amin(z[0])
zmin_index_big = np.where(z[0] == zmin_value)[0]

xmax = np.amax(x[0,zmin_index_big])
zmin_index = np.where((z[0] == zmin_value) & (x[0] <=xmax) & (x[0] >= (xmax-dx_init)))[0]
print(len(zmax_index), len(zmin_index))
#Try to make sure the number of atoms selected in the bottom later are similar to the number in the top layer and span a similar range.

#S matrix is length along the folded sheet
S = np.zeros(num_atoms)



for atom in range(num_atoms-1, -1,-1):
	#print(atom)
	if (atom in zmin_index_big):
		S[atom] = x[0,atom]
	else:
		dr = math.sqrt((x[0,atom]-x[0,atom+1])**2 + (z[0,atom]-z[0,atom+1])**2)
		if (dr > 2*bond_length):
			if (atom in zmax_index): #on top layer
				#print((zmax_index > atom))
				ds_min = np.amin(abs(x[0,atom]-x[0,zmax_index[zmax_index > atom]]))
				ds_in = np.where(abs(x[0,atom]-x[0,atom+1:]) == ds_min)[0]
				ds_sign = np.sign(x[0,atom]-x[0,atom+ds_in])
				#print(ds_min,abs(x[0,atom]-x[0,atom+1+ds_in]), abs(x[0,atom]-x[0,ds_in+1]),ds_in,ds_sign)
			else:
				ds_min = np.amin(abs(z[0,atom]-z[0,atom+1:]))
				ds_in = np.where(abs(z[0,atom]-z[0,atom+1:]) == ds_min)[0]
				ds_sign = np.sign(z[0,atom]-z[0,atom+ds_in])
				#print(ds_min,abs(z[0,atom]-z[0,atom+ds_in]), abs(z[0,atom]-z[0,atom+ds_in+1]))

			ds_comp = S[atom+ds_in+1]
			ds = math.sqrt((x[0,atom] - x[0,atom+ds_in+1])**2 + (z[0,atom]-z[0,atom+ds_in+1])**2)
			S[atom] = ds_comp + ds_sign*ds
		else:
			dz = z[0,atom] - z[0,atom+1]
			dx = x[0,atom] - x[0,atom+1]
			if (dz == 0.0):#Checking direction of atoms in top layer
				if (dx > 0):
					ds_sign = 1.0
				else:
					ds_sign = -1.0
			elif (dz > 0.0): #Checking direction of atoms in the fold
				ds_sign = 1.0
			else:
				ds_sign = -1.0
			S[atom] = S[atom+1]-dr*ds_sign
			#print(ds_sign)
	#print(x[0,atom],z[0,atom])
	#print(S[atom],-0.1)
#	if (y[0,atom] < 2.0):
#		fid_S.write("%g \t %g \n" % (x[0,atom],z[0,atom]))
#		fid_S.write("%g \t %g \n" % (S[atom], -0.1))

Sm = np.amax(S)


filename_dist = 'Distance.dat'
filename_vel = 'velocity2.dat'
filename_angle = 'angle.dat'
fid_dist = open(filename_dist, 'w')
fid_vel = open(filename_vel, 'w')
fid_ang = open(filename_angle, 'w')

#fid_dist.write("#Time(ps) \t spacing 3SD \t spacing 2SD \t spacing 1SD \t SD \t  overhang 3SD \t overhang 2SD \t overhang 1SD \t SD\n")
fid_dist.write("#Time(ps) \t layer spacing \t overhang \n ")

fid_vel.write("#Time(ps) \t vx(A/ps) \t vy(A/ps) \t vz(A/ps) \t v(A/ps) \n")
fid_ang.write("#Time(ps) \t phi \t v*cos(phi)\n")
Sarray = np.zeros(max(len(zmax_index), len(zmin_index)))
rarray = np.zeros(max(len(zmax_index), len(zmin_index)))
#vxarray = np.zeros(max(len(zmax_index), len(zmin_index)))
#vyarray = np.zeros(max(len(zmax_index), len(zmin_index)))
#vzarray = np.zeros(max(len(zmax_index), len(zmin_index)))
vxtarray = np.zeros(max(len(zmax_index), len(zmin_index)))
vytarray = np.zeros(max(len(zmax_index), len(zmin_index)))
vztarray = np.zeros(max(len(zmax_index), len(zmin_index)))

vxbarray = np.zeros(max(len(zmax_index), len(zmin_index)))
vybarray = np.zeros(max(len(zmax_index), len(zmin_index)))
vzbarray = np.zeros(max(len(zmax_index), len(zmin_index)))



s_2SD_ave = 1
for i, dt in enumerate(t):
	if ( (i % 100)== 0):
		print("Time",i,dt)
	#Average positions of atoms in top layer
	x_top = np.mean(x[i,zmax_index])
	y_top = np.mean(y[i,zmax_index])
	z_top = np.mean(z[i,zmax_index])
	#Range of positions for atoms in the top layer
	dx_top = np.amax(x[i,zmax_index]) - np.amin(x[i,zmax_index])
	dy_top = np.amax(y[i,zmax_index]) - np.amin(y[i,zmax_index])
	dz_top = np.amax(z[i,zmax_index]) - np.amin(z[i,zmax_index])

	#Average positions of atoms in bottom layer
	x_bottom = np.mean(x[i,zmin_index])
	y_bottom = np.mean(y[i,zmin_index])
	z_bottom = np.mean(z[i,zmin_index])
	
	dx_sum =0.0
	dy_sum =0.0
	dz_sum =0.0
	S_sum =0.0
	r_sum =0.0
	if (s_2SD_ave > 0.0):
		#bottom layer is longer than top layer
		for jj, atom in enumerate(list(zmax_index)):
			r2 = (x[i,atom] - x[i,:])**2 + (y[i,atom] - y[i,:])**2 + (z[i,atom] - z[i,:])**2
			dS2 = (S - S[atom])**2
			#print(jj,atom,r2,dS2)
			r2 = np.where(dS2 <= 2*r2, 10000.0, r2)
			#print(r2,dS2)
			# Fixes so that atoms on the same plane are removed from the minimisation check
			r2min = np.amin(r2)
			if (r2min == 10000.0):
				rarray[jj] = 10000.0
				Sarray[jj] = 10000.0
				continue
				#Should be removed in standard deviation test
			r2min_index = np.where(r2 == r2min)[0]
			#print(dS2min,r2min, r2min_index)
			rmin = math.sqrt(r2min)
			dSmin = Sm - S[r2min_index] - S[atom]
			#print(i, S[r2min_index], Sm, S[atom],dSmin)
			#print(rmin,dSmin)
			#dx_sum += x[i,atom] - x[i,r2min_index]
			#dy_sum += y[i,atom] - y[i,r2min_index]
			#dz_sum += z[i,atom] - z[i,r2min_index]
			#S_sum += dSmin
			#r_sum += rmin
			#print(r2min,dSmin, jj, r2min_index)
			Sarray[jj] = dSmin
			rarray[jj] = rmin
			#vxarray[jj] = vx[i,r2min_index] - vx[i,atom]
			#vyarray[jj] = vy[i,r2min_index] - vy[i,atom]
			#vzarray[jj] = vz[i,r2min_index] - vz[i,atom]
			vxtarray[jj] = vx[i,atom]
			vxbarray[jj] = vx[i,r2min_index]
			vytarray[jj] = vy[i,atom]
			vybarray[jj] = vy[i,r2min_index]
			vztarray[jj] = vz[i,atom]
			vzbarray[jj] = vz[i,r2min_index]


	else:
		#top layer is longer than bottom layer
		for jj, atom in enumerate(list(zmin_index)):
			r2 = (x[i,atom] - x[i,:])**2 + (y[i,atom] - y[i,:])**2 + (z[i,atom] - z[i,:])**2
			dS2 = (S - S[atom])**2
			#print(jj,atom,r2,dS2)
			r2 = np.where(dS2 <= 2*r2, 10000.0, r2)
			#print(r2,dS2)
			# Fixes so that atoms on the same plane are removed from the minimisation check
			r2min = np.amin(r2)
			if (r2min == 10000.0):
				rarray[jj] = 10000.0
				Sarray[jj] = 10000.0
				continue
				#Should be removed in standard deviation test

			r2min_index = np.where(r2 == r2min)[0]
			dS2min = dS2[r2min_index]
			#print(dS2min,r2min, r2min_index)
			rmin = math.sqrt(r2min)
			dSmin = -(S[r2min_index][0] - Sm + S[atom])
			#print(i,S[r2min_index], Sm, S[atom],dSmin)
			#print(rmin,dSmin)
			#dx_sum += x[i,atom] - x[i,r2min_index]
			#dy_sum += y[i,atom] - y[i,r2min_index]
			#dz_sum += z[i,atom] - z[i,r2min_index]
			#S_sum += dSmin
			#r_sum += rmin
			Sarray[jj] = dSmin
			rarray[jj] = rmin
			#vxarray[jj] = vx[i,atom] - vx[i,r2min_index]
			#vyarray[jj] = vy[i,atom] - vy[i,r2min_index]
			#vzarray[jj] = vz[i,atom] - vz[i,r2min_index]
			vxbarray[jj] = vx[i,atom]
			vxtarray[jj] = vx[i,r2min_index]
			vybarray[jj] = vy[i,atom]
			vytarray[jj] = vy[i,r2min_index]
			vzbarray[jj] = vz[i,atom]
			vztarray[jj] = vz[i,r2min_index]


	r_ave1 = np.mean(rarray)
	r_SD = np.std(rarray)
	r_dist_from_mean = abs(rarray - r_ave1)
	#r_3SD_inc = r_dist_from_mean < 3*r_SD
	r_2SD_inc = r_dist_from_mean < 2*r_SD
	#r_1SD_inc = r_dist_from_mean < r_SD

	#r_3SD = rarray[r_3SD_inc]
	r_2SD = rarray[r_2SD_inc]
	#r_1SD = rarray[r_1SD_inc]

	#r_3SD_ave = np.mean(r_3SD)
	r_2SD_ave = np.mean(r_2SD)
	#r_1SD_ave = np.mean(r_1SD)
		
	#Without first set of outliers (2SD)
	#r_2SD_ave = np.mean(r_2SD)
	#r_SD = np.std(r_2SD)
	#r_dist_from_mean = abs(r_2SD - r_2SD_ave)
	#r_2SD_inc = r_dist_from_mean < 2*r_SD
	#r1_2SD = r_2SD[r_2SD_inc]
	#r1_ave = np.mean(r1_2SD)

	#Without second set of outliers (2SD)
	#r_2SD_ave = np.mean(r1_2SD)
	#r_SD = np.std(r1_2SD)
	#r_dist_from_mean = abs(r1_2SD - r_2SD_ave)
	#r_2SD_inc = r_dist_from_mean < 2*r_SD
	#r2_2SD = r1_2SD[r_2SD_inc]
	#r2_ave = np.mean(r2_2SD)


	s_ave1 = np.mean(Sarray)
	s_SD = np.std(Sarray)
	s_dist_from_mean = abs(Sarray - s_ave1)
	#s_3SD_inc = s_dist_from_mean < 3*s_SD
	s_2SD_inc = s_dist_from_mean < 2*s_SD
	#s_1SD_inc = s_dist_from_mean < s_SD

	#s_3SD = Sarray[s_3SD_inc]
	s_2SD = Sarray[s_2SD_inc]
	#s_1SD = Sarray[s_1SD_inc]

	#s_3SD_ave = np.mean(s_3SD)
	s_2SD_ave = np.mean(s_2SD)
	#s_1SD_ave = np.mean(s_1SD)

	#Without first set of outliers (2SD)
	#s_2SD_ave = np.mean(s_2SD)
	#s_SD = np.std(s_2SD)
	#s_dist_from_mean = abs(s_2SD - s_2SD_ave)
	#s_2SD_inc = s_dist_from_mean < 2*s_SD
	#s1_2SD = s_2SD[s_2SD_inc]
	#s1_ave = np.mean(s1_2SD)

	#Without second set of outliers (2SD)
	#s_2SD_ave = np.mean(s1_2SD)
	#s_SD = np.std(s1_2SD)
	#s_dist_from_mean = abs(s1_2SD - s_2SD_ave)
	#s_2SD_inc = s_dist_from_mean < 2*s_SD
	#s2_2SD = s1_2SD[s_2SD_inc]
	#s2_ave = np.mean(s2_2SD)



	#Velocity check
	vxarray = vxbarray - vxtarray
	vx_ave_raw = np.mean(vxarray)
	vx_SD = np.std(vxarray)
	vx_dist_from_mean = abs(vxarray - vx_ave_raw)
	vx_2SD_inc = vx_dist_from_mean < 2*vx_SD
	vx_2SD = vxarray[vx_2SD_inc]
	vx_ave = np.mean(vx_2SD)

	vyarray = vybarray - vytarray
	vy_ave_raw = np.mean(vyarray)
	vy_SD = np.std(vyarray)
	vy_dist_from_mean = abs(vyarray - vy_ave_raw)
	vy_2SD_inc = vy_dist_from_mean < 2*vy_SD
	vy_2SD = vyarray[vy_2SD_inc]
	vy_ave = np.mean(vy_2SD)

	vzarray = vzbarray-vztarray
	vz_ave_raw = np.mean(vzarray)
	vz_SD = np.std(vzarray)
	vz_dist_from_mean = abs(vzarray - vz_ave_raw)
	vz_2SD_inc = vz_dist_from_mean < 2*vz_SD
	vz_2SD = vzarray[vz_2SD_inc]
	vz_ave = np.mean(vz_2SD)

	vmag = math.sqrt(vx_ave**2 + vy_ave**2 +vz_ave**2)


	#vxt_array = vx[i,zmax_index]
	#vxt_ave = np.mean(vxt_array)
	#vxt_SD = np.std(vxt_array)
	#vxt_dfm = abs(vxt_array - vxt_ave)
	#vxt_2SD_inc = vxt_dfm < 2*vxt_SD
	#vxt_2SD = vxt_array[vxt_2SD_inc]
	#vx_top = np.mean(vxt_2SD)

	#vyt_array = vy[i,zmax_index]
	#vyt_ave = np.mean(vyt_array)
	#vyt_SD = np.std(vyt_array)
	#vyt_dfm = abs(vyt_array - vyt_ave)
	#vyt_2SD_inc = vyt_dfm < 2*vyt_SD
	#vyt_2SD = vyt_array[vyt_2SD_inc]
	#vy_top = np.mean(vyt_2SD)

	#vzt_array = vz[i,zmax_index]
	#vzt_ave = np.mean(vzt_array)
	#vzt_SD = np.std(vzt_array)
	#vzt_dfm = abs(vzt_array - vzt_ave)
	#vzt_2SD_inc = vzt_dfm < 2*vzt_SD
	#vzt_2SD = vzt_array[vzt_2SD_inc]
	#vz_top = np.mean(vzt_2SD)

	#vxb_array = vx[i,zmin_index]
	#vxb_ave = np.mean(vxb_array)
	#vxb_SD = np.std(vxb_array)
	#vxb_dfm = abs(vxb_array - vxb_ave)
	#vxb_2SD_inc = vxb_dfm < 2*vxb_SD
	#vxb_2SD = vxb_array[vxb_2SD_inc]
	#vx_bottom = np.mean(vxb_2SD)

	#vyb_array = vy[i,zmin_index]
	#vyb_ave = np.mean(vyb_array)
	#vyb_SD = np.std(vyb_array)
	#vyb_dfm = abs(vyb_array - vyb_ave)
	#vyb_2SD_inc = vyb_dfm < 2*vyb_SD
	#vyb_2SD = vyb_array[vyb_2SD_inc]
	#vy_bottom = np.mean(vyb_2SD)

	#vzb_array = vz[i,zmin_index]
	#vzb_ave = np.mean(vzb_array)
	#vzb_SD = np.std(vzb_array)
	#vzb_dfm = abs(vzb_array - vzb_ave)
	#vzb_2SD_inc = vzb_dfm < 2*vzb_SD
	#vzb_2SD = vzb_array[vzb_2SD_inc]
	#vz_bottom = np.mean(vzb_2SD)


	#vx_top = np.mean(vx[i,zmax_index])
	#vy_top = np.mean(vy[i,zmax_index])
	#vz_top = np.mean(vz[i,zmax_index])

	#vx_bottom = np.mean(vx[i,zmin_index])
	#vy_bottom = np.mean(vy[i,zmin_index])
	#vz_bottom = np.mean(vz[i,zmin_index])
    	

	#vx_ave = vx_top - vx_bottom
	#vy_ave = vy_top - vy_bottom
	#vz_ave = vz_top - vz_bottom
	vx_top = np.mean(vxtarray)
	vx_bottom = np.mean(vxbarray)
	vy_top = np.mean(vytarray)
	vy_bottom = np.mean(vybarray)
	vz_top = np.mean(vztarray)
	vz_bottom = np.mean(vzbarray)

	vdot = vx_top*vx_bottom + vy_top*vy_bottom + vz_top*vz_bottom
	vtmag = math.sqrt(vx_top**2 + vy_top**2 +vz_top**2)
	vbmag = math.sqrt(vx_bottom**2 + vy_bottom**2 +vz_bottom**2)

	ang = math.acos(-vdot/(vtmag*vbmag))
#	if ( (i % 100)== 0):
#		print(vdot, vmag, ang)

	vabs = abs(vx_ave)/vx_ave*vmag
#	print(vx_ave, vy_ave, vz_ave, vabs)
	#fid_dist.write("%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n"%(dt*timestep,r_3SD_ave, r_2SD_ave, r_1SD_ave, r_SD, s_3SD_ave, s_2SD_ave, s_1SD_ave, s_SD ))
	#fid_dist.write("%g \t %g \t %g \t %g \t %g \t %g \t %g \n"%(dt*timestep,r_ave1, r_2SD_ave, r1_ave, s_ave1, s_2SD_ave, s1_ave))
	fid_dist.write("%g \t %g \t %g\n"%(dt*timestep, r_2SD_ave, s_2SD_ave))
	fid_vel.write("%g \t %g \t %g \t %g \t %g \n"% (dt*timestep, vx_ave*1000,vy_ave*1000,vz_ave*1000,vabs*1000))
	fid_ang.write("%g \t %g \t %g \n" % (dt*timestep,ang, 1000*vabs*math.cos(ang)))
fid_dist.close()
fid_vel.close()
fid_ang.close()






