#Python version
#Calls relevant folding scripts for 0, 30 and arbitrary angles
# INPUT PARAMETERS:
fileID = open("param_g",'r');
param = fileID.readlines()
ftheta = float(param[0])# Fold angle
if (ftheta < 0):
	print("Fold angle must be greater than 0")
	exit()
lx = float(param[1]) ;    # Length of sheet in the x direction
ly = float(param[2]);     # Length of sheet in the y direction
if (lx <= 0 and ly <= 0):
	print("dimensions of grid need to be positive")
	exit()
r1 = float(param[3]); #inner radius
if (r1 < 0):
	print("inner radius must be non-negative")
	exit()
r2 = float(param[4]); # length of top layer
if (r2 < 0):
	print(" Top layer must be non-negative")
	exit()
theta = float(param[5]); #Rotation of graphene sheet
if (theta < 0 or theta > 30):
	print("Angle of rotation for graphene sheet relative to scroll is between 0 and 30 degrees")
	exit()
periodic_y = int(param[6]) #Switch for periodic set-up
#if (ftheta > 0.0):
#    periodic_x = 0
    #Can't be periodic in the x directions
#else: 
periodic_x = int(param[7]) 
if (periodic_y == 0):
	y_extra = float(param[8])
else:
	y_extra = 0
#Added width for tearing simulations. Set to zero for periodic calculations
if (y_extra < 0):
	print("Tearing addition must be greater than zero")
	exit()
elif (y_extra ==0):
	print("No tearing")
s_lx = float(param[9])
s_ly = float(param[10])
s_lz = float(param[11])
if ((s_lx < 0) or (s_ly < 0) or (s_lz < 0)):
	print("Substrate dimensions must be non-negative")
	exit()
elif ((((s_lx < lx) and (s_lx>0)) or ((s_ly < ly) and (s_ly >0))) and (s_lz > 0)):
	print("Substrate should be larger than graphene sheet")
	exit()
elif ((s_lx*s_ly*s_lz) == 0):
	print("No substrate")
h2o = float(param[12])
if (h2o == 1):
	print("Water will be present")	
Hatoms = float(param[13])
if (Hatoms == 1):
    print("Edge H atoms included")
layer_num = int(param[14])
if (layer_num == 1):
    print("Single layer graphene")
elif (layer_num == 2):
    print("Bilayer graphene")
else:
    print("Too many layers")
    exit()
fileID.close()
filename = 'graphene.dat';
filename2 = 'graphene.xyz';
writeresults = True;
#showresults = True;


import folded_otherangles
