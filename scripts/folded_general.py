#Python version
#Calls relevant folding scripts for 0, 30 and arbitrary angles
# INPUT PARAMETERS:
fileID = open("param_g",'r');
param = fileID.readlines()
a = float(param[0])# interatomic distance
if (a <= 0):
	print("Interatomic distance between carbon atoms must be greater than zero")
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
periodic = int(param[6]) #Switch for periodic set-up
yextra = 0
if (len(param) > 7):
	yextra = int(param[7]) #length of unfolded edge pieces
	if (yextra < 0):
		print("Extra y must be non-negative")
		exit()
fileID.close()
filename = 'graphene.dat';
writeresults = True;
#showresults = True;


if ((theta == 0)| (theta == 30)): #Zigzag or armchair arrangements
	import folded_arm_zig
else:
	import folded_otherangles
