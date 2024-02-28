#!/usr/bin/env python
# coding: utf-8
#Extracts distance and velocity information from LAMMPS output files


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
#Extracts timestep and converts to ps


contents = []
with open('graphene_xyz_v.lammpstrj','r') as file:
	contents = file.readlines()
line_start = 1
line_num = 0
pi = math.pi
bond_length = 1.4
num_atoms = contents[3] #Total number of atoms
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
filename_vel = 'velocity.dat'
fid_dist = open(filename_dist, 'w')
fid_vel = open(filename_vel, 'w')

fid_dist.write("#Time(ps) \t layer spacing \t overhang \n ")
fid_vel.write("#Time(ps) \t vx(A/ps) \t vy(A/ps) \t vz(A/ps) \t v(A/ps) \n")

Sarray = np.zeros(max(len(zmax_index), len(zmin_index)))
rarray = np.zeros(max(len(zmax_index), len(zmin_index)))

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
			
			Sarray[jj] = dSmin
			rarray[jj] = rmin
			
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
			
			#S_sum += dSmin
			#r_sum += rmin
			Sarray[jj] = dSmin
			rarray[jj] = rmin
			
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

	vx_top = np.mean(vxtarray)
	vx_bottom = np.mean(vxbarray)
	vy_top = np.mean(vytarray)
	vy_bottom = np.mean(vybarray)
	vz_top = np.mean(vztarray)
	vz_bottom = np.mean(vzbarray)

	vabs = abs(vx_ave)/vx_ave*vmag

	fid_dist.write("%g \t %g \t %g\n"%(dt*timestep, r_2SD_ave, s_2SD_ave))
	fid_vel.write("%g \t %g \t %g \t %g \t %g \n"% (dt*timestep, vx_ave*1000,vy_ave*1000,vz_ave*1000,vabs*1000))
	
fid_dist.close()
fid_vel.close()
