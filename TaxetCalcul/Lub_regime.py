################################################################################################
''' Developed by Kalenia Marquez Florez
    Aix Marseille Univ, CNRS, ISM, Marseille, France
    March 2023

    Script to calculate the film thickness of the elbow in extreme locomotion conditions.
'''
#

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import numpy as np
import math

# Type class for the definition of the parameters of contact
class Cas():
    def __init__(self, eta, E, nu, w, V, Rm):
        self.eta = eta
        self.E = E
        self.nu = nu
        self.w = w
        self.V = V
        self.Rm = Rm

# Function to calculate the dimensionless thickness
def AdHemin(eta, E, nu, w, Ra, Rb, V, Rm, k=50, Ru = 'mm'):
    Ef = [E/(1-nu**2), 0.0] # effective Young's Modulus PA
    Rf = [1/((1/Ra[0])-(1/Rb[0])), Ra[1]] # effective Radius mm
    Vp = [Ra[0]*V[0]/Rm[0], V[1]] # sliding speed at cartilage m/s
    if Ru == 'mm': # to use if radius in mm
        Rf[0] = Rf[0]*10**-3 # m
    # dimensionless velocity
    U = [(eta[0]*Vp[0])/(Ef[0]*Rf[0]), (eta[1]+Vp[1])-(Ef[1]+Rf[1])] 
    # dimensionless load
    W = [w[0]/(Ef[0]*Rf[0]**2),w[1]-(Ef[1]+2*Rf[1])]
    # dimensionless thickness
    Hemin = [7.43*(U[0]**0.65)*(W[0]**-0.21), 0.65*U[1]-0.21*W[1]]
    # *(1-0.85*math.exp(-0.31*k))
    return(Hemin, Rf)

# Function tu calculate the min film thickness
def dhemin(Hemin, Rf):
    # Hemin : dimensionless thickness
    # Rf: effective radius m
    # hemin:  min thickness um
    hemin = [Hemin[0] * Rf[0] * 10**6, Hemin[1]+Rf[1]] #um
    return(hemin)

# Function to calculate the radii with the cartilage and the opposing surface
def Rab(Rm, ci, cii):
    # ci: proportion for cartilage thickness
    # cii: proportion for the opposing surface radius
    # Ra: radius with cartilage mm
    # Rb: radius of the opposing surface mm
    Ra = [Rm[0]*(1+ci), Rm[1]]
    Rb = [Rm[0]*(1+ci+cii), Rm[1]]
    return(Ra, Rb)

# Function to calculate the min film thickness 
def Chemin(ci, cii, Data, val = False):
    # ci: proportion for cartilage thickness
    # cii: proportion for the opposing surface radius
    # Ra: radius with cartilage
    # Rb: radius of the opposing surface
    Ra, Rb = Rab(Data.Rm, ci, cii)
    # Hemin: dimensionless thickness
    # Rf: effective radius
    Hemin, Rf = AdHemin(Data.eta, Data.E, Data.nu, Data.w, Ra, Rb, Data.V, Data.Rm)
    # hemin:  min thickness um
    hemin = dhemin(Hemin, Rf)
    if val == False: # get both allometric coef. and exp.
        return(hemin)
    if val == True: # just get the allometric coef.
        return(hemin[0])

# Function to calculate Lambda for the classification of the regime
def Lambda(Ru, hmin):
    # Lambda: Lambda ratio for the classification of lub. regime
    # Ru: Cartilage rugosity um
    # hemin:  min thickness um
    Lambda = hmin / (Ru*2**0.5)
    return(Lambda)

# Function to assign a color according to lambda 
def Col(Lmb):
    # Classical classification
    if Lmb < 1: # boundary
        Colores = (0.247, 0.318, 0.710, 0.0) 
    if Lmb >= 1 and Lmb<3: # Mixed
        Colores= (0.957, 0.263, 0.212, 1.0)
    if Lmb>=3 and Lmb<5: #elasto
        Colores = (0.306, 0.784, 0.314, 1.0)
    if Lmb>=5: #Hydro
        Colores = (1.000, 0.757, 0.027, 1.0)
    ######################################
    # New classification
    # if Lmb < 0.1: # boundary
    #     Colores = (0.247, 0.318, 0.710, 0.0) 
    # if Lmb < 1 and Lmb<3: # Mixed
    #     Colores= (0.957, 0.263, 0.212, 1.0)
    # if Lmb>=1 and Lmb<5: #elasto
    #     Colores = (0.306, 0.784, 0.314, 1.0)
    # if Lmb>=5: #Hydro
    #     Colores = (1.000, 0.757, 0.027, 1.0)
    #
    return Colores

#3F51B5 - a deep blue color with RGBA code (0.247, 0.318, 0.710, 1.0)
#F44336 - a bright red color with RGBA code (0.957, 0.263, 0.212, 1.0)
#4CAF50 - a bright green color with RGBA code (0.306, 0.784, 0.314, 1.0)
#FFC107 - a bright yellow color with RGBA code (1.000, 0.757, 0.027, 1.0)

# Function to calculate Lambda for the classification of the regime
# With ci and cii as inputs
def LambdaCII(Ru, ci, cii, Data):
    # Ru: Cartilage rugosity um
    # ci: proportion for cartilage thickness
    # cii: proportion for the opposing surface radius
    # Data: object with the basic data of the problem 
    hemin = Chemin(ci, cii, Data, val = True) # min thickness um
    Lmb = Lambda(Ru, hemin) # Lambda ratio
    col = Col(Lmb) # color for classification
    return Lmb, col

# Function to calculate Lambda for the classification of the regime
# With ci and cii as inputs and changing the sliding velocity
def LambdaVts(Ru, ci, cii, Data, Vn):
    # Ru: Cartilage rugosity um
    # ci: proportion for cartilage thickness
    # cii: proportion for the opposing surface radius
    # Data: object with the basic data of the problem 
    # Vn: new sliding velocity
    Data2 = Cas(Data.eta, Data.E, Data.nu, Data.w, Vn, Data.Rm)
    # Data2: new object with the new information
    hemin = Chemin(ci, cii, Data2, val = True) # min thickness
    Lmb = Lambda(Ru, hemin) # Lambda ratio
    col = Col(Lmb) # color for classification
    return Lmb, col

# Function to calculate Lambda for the classification of the regime
# With ci and cii as inputs and changing the sliding velocity and the 
# force accordingly 
def LambdaVtsFs(Ra, ci, cii, Data, Vn):
    wn = [Data.w[0]*Vn[0]/Data.V[0], Data.w[1]]
    Data2 = Cas(Data.eta, Data.E, Data.nu, wn, Vn, Data.Rm)
    hemin = Chemin(ci, cii, Data2, val = True)
    Lmb = Lambda(Ra, hemin)
    col = Col(Lmb)
    return Lmb, col
################################################

################################################
# general variables
eta = [0.01, 0] # Viscosity Pa-s
nu = 0.4 # Poisson's ratio
E = 8.1*10**6 # Young's modulus Pa
w = [78.2, 0.78] # Muscle force N
V = [10**0.61, -0.07] # Sliding speed m/s
Rm = [(10**0.6)/2, 0.35] # Average radius mm 
Data = Cas(eta, E, nu, w, V, Rm) # Creation of the object Data

################################################
######### FILM THICKNESS EVALUATION ############

# TEST for a specific ci and cii
ci = 0.06
cii = 0.07
print(Chemin(ci, cii, Data))

# Print plot for ci and cii sensitivity
ci = np.linspace(0.0, 0.2, 1000)
cii = np.geomspace(0.01, 1, num=2000)
Ci, Cii = np.meshgrid(ci, cii)

hee = Chemin(Ci, Cii, Data, val = True)

# plot the surface
fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(111, projection='3d')
surf = ax1.plot_surface(Ci, Cii, hee, cmap='viridis',  linewidth=0, antialiased=False, alpha=0.8)

ax1.set_xlabel('ci')
ax1.set_ylabel('cii')
ax1.set_zlabel('hmin')

z_min = np.min(surf.get_array())
z_max = np.max(surf.get_array())

print('h_min the min:', z_min)
print('h_min the max:', z_max)

v = np.linspace(z_min, z_max, 10)
fig.colorbar(surf, ticks = v)

plt.show()

##############################################
############ LAMBDA EVALUATION ###############
# To evaluate just one value 
print(LambdaCII(1.42, 0.06, 0.07, Data))

# Plot parameters
Az = -70
elv = 27
ctm = 1/2.54  # centimeters in inches
szx = 20
szy = 20

## Uncomment for evaluation in a wide range 
# c_i = np.linspace(0.0, 1.0, 20 ) #5
# c_ii = np.linspace(0.01, 10, 20) #10
# ra = np.linspace(0.68, 80, 20)#

## Uncomment for evaluation in a range near elbow param 
c_i = np.linspace(0.0, 0.1, 25 ) #5
c_ii = np.linspace(0.01, 0.2, 25) #10
ru = np.linspace(0.68, 6, 25)#

# Create vectors
Ci = []
Cii = []
Ru = []
Lmbd = []
Colors = []
for Ici in c_i:
    for Icii in c_ii:
        for Iru in ru:
            Ci.append(Ici)
            Cii.append(Icii)
            Ru.append(Iru)
            Lmb, col = LambdaCII(Iru, Ici, Icii, Data)
            Lmbd.append(Lmb)
            Colors.append(col)

# create figure
fig3 = plt.figure(figsize=(szx*ctm, szy*ctm))
ax3 = plt.axes(projection='3d')

# plot the cube of points
p = ax3.scatter(Ru, Ci, Cii, 
                c = Colors, 
                cmap= 'viridis', 
                marker = "o", 
                alpha= 0.5,
                edgecolors= None,
                linewidths = 0,
                )

# boundaries of the cube
ax3.plot([min(Ru), max(Ru)], [min(Ci), min(Ci)], [min(Cii), min(Cii)], color='grey') # plot the line
ax3.plot([max(Ru), max(Ru)], [min(Ci), max(Ci)], [min(Cii), min(Cii)], color='grey') # plot the line
ax3.plot([max(Ru), min(Ru)], [max(Ci), max(Ci)], [min(Cii), min(Cii)], color='grey') # plot the line
ax3.plot([min(Ru), min(Ru)], [max(Ci), min(Ci)], [min(Cii), min(Cii)], color='grey') # plot the line

ax3.plot([min(Ru), max(Ru)], [min(Ci), min(Ci)], [max(Cii), max(Cii)], color='grey') # plot the line
ax3.plot([max(Ru), max(Ru)], [min(Ci), max(Ci)], [max(Cii), max(Cii)], color='grey') # plot the line
ax3.plot([max(Ru), min(Ru)], [max(Ci), max(Ci)], [max(Cii), max(Cii)], color='grey') # plot the line
ax3.plot([min(Ru), min(Ru)], [max(Ci), min(Ci)], [max(Cii), max(Cii)], color='grey') # plot the line

ax3.plot([max(Ru), max(Ru)], [min(Ci), min(Ci)], [min(Cii), max(Cii)], color='grey') # plot the line
ax3.plot([max(Ru), max(Ru)], [max(Ci), max(Ci)], [min(Cii), max(Cii)], color='grey') # plot the line
ax3.plot([min(Ru), min(Ru)], [max(Ci), max(Ci)], [min(Cii), max(Cii)], color='grey') # plot the line
ax3.plot([min(Ru), min(Ru)], [min(Ci), min(Ci)], [min(Cii), max(Cii)], color='grey') # plot the line

Xmin = min(Ru)
Xmax = max(Ru)
Ymin = min(Ci)
Ymax = max(Ci)
Zmin = min(Cii)
Zmax = max(Cii)

ax3.set_xlim([Xmin, Xmax])
ax3.set_ylim([Ymin, Ymax])
ax3.set_zlim([Zmin, Zmax])

ax3.set_xlabel('Ru')
ax3.set_ylabel('Ci')
ax3.set_zlabel('Cii')

# orientation
ax3.view_init(azim= Az, elev=elv)
plt.show()

#### REMOVAL OF THE REGIMES 

# Uncomment if more density is needed 
# c_i = np.linspace(0.0, 1.0, 70 ) #5
# c_ii = np.linspace(0.01, 10, 70) #10
# ru = np.linspace(0.68, 80, 70)#(0.68, 6.5, 5) #10
# Ci = []
# Cii = []
# Ru = []
# Lmbd = []
# Colors = []
# for Ici in c_i:
#     for Icii in c_ii:
#         for Iru in ru:
#             Ci.append(Ici)
#             Cii.append(Icii)
#             Ru.append(Iru)
#             Lmb, col = LambdaCII(Iru, Ici, Icii, Data)
#             Lmbd.append(Lmb)
#             Colors.append(col)

# Create the figure
fig3 = plt.figure(figsize=(szx*ctm, szy*ctm))
ax3 = plt.axes(projection='3d')

# BOUNDARY
# Removal of the values of BOUNDARY
Lmbd_sans_bd =[]
Ci_sans_bd = [] 
Cii_sans_bd = []
Ru_sans_bd = []
Colors_sans_bd = []
for ILmb in range(len(Lmbd)):
    if Lmbd[ILmb] >= 1:
        Lmbd_sans_bd.append(Lmbd[ILmb])
        Ci_sans_bd.append(Ci[ILmb])
        Cii_sans_bd.append(Cii[ILmb])
        Ru_sans_bd.append(Ru[ILmb])
        Colors_sans_bd.append(Colors[ILmb])

# plot
p = ax3.scatter(Ru_sans_bd, Ci_sans_bd, Cii_sans_bd, 
                c =  Colors_sans_bd, 
                cmap= 'viridis', 
                marker = ".", 
                #s=Lmbd
                alpha= 0.5,
                edgecolors= None,
                linewidths = 0,
                )

ax3.set_xlim([Xmin, Xmax])
ax3.set_ylim([Ymin, Ymax])
ax3.set_zlim([Zmin, Zmax])

ax3.set_xlabel('Ru')
ax3.set_ylabel('Ci')
ax3.set_zlabel('Cii')

# orientation
ax3.view_init(azim= Az, elev=elv)
plt.show()

# Create the figure
fig3 = plt.figure(figsize=(szx*ctm, szy*ctm))
ax3 = plt.axes(projection='3d')

# MIXED
# Removal of the values of MIXED
Lmbd_sans_bdM =[]
Ci_sans_bdM = [] 
Cii_sans_bdM = []
Ru_sans_bdM = []
Colors_sans_bdM = []
for ILmb in range(len(Lmbd_sans_bd)):
    if Lmbd_sans_bd[ILmb] >= 3:
        Lmbd_sans_bdM.append(Lmbd_sans_bd[ILmb])
        Ci_sans_bdM.append(Ci_sans_bd[ILmb])
        Cii_sans_bdM.append(Cii_sans_bd[ILmb])
        Ru_sans_bdM.append(Ru_sans_bd[ILmb])
        Colors_sans_bdM.append(Colors_sans_bd[ILmb])

# plot
p = ax3.scatter(Ru_sans_bdM, Ci_sans_bdM, Cii_sans_bdM, 
                c =  Colors_sans_bdM, 
                cmap= 'viridis', 
                marker = ".", 
                #s=Lmbd
                alpha= 0.5,
                edgecolors= None,
                linewidths = 0,
                )

ax3.set_xlim([Xmin, Xmax])
ax3.set_ylim([Ymin, Ymax])
ax3.set_zlim([Zmin, Zmax])

ax3.set_xlabel('Ru')
ax3.set_ylabel('Ci')
ax3.set_zlabel('Cii')

# orientation
ax3.view_init(azim= Az, elev=elv)
plt.show()

# Create the figure
fig3 = plt.figure(figsize=(szx*ctm, szy*ctm))
ax3 = plt.axes(projection='3d')

# ELASTOHYDRODYNAMIC
# Removal of the values of ELASTOHYDRODYNAMIC
Lmbd_sans_bdME =[]
Ci_sans_bdME = [] 
Cii_sans_bdME = []
Ru_sans_bdME = []
Colors_sans_bdME = []
for ILmb in range(len(Lmbd_sans_bdM)):
    if Lmbd_sans_bdM[ILmb] >= 5:
        Lmbd_sans_bdME.append(Lmbd_sans_bdM[ILmb])
        Ci_sans_bdME.append(Ci_sans_bdM[ILmb])
        Cii_sans_bdME.append(Cii_sans_bdM[ILmb])
        Ru_sans_bdME.append(Ru_sans_bdM[ILmb])
        Colors_sans_bdME.append(Colors_sans_bdM[ILmb])

# plot
p = ax3.scatter(Ru_sans_bdME, Ci_sans_bdME, Cii_sans_bdME, 
                c =  Colors_sans_bdME, 
                cmap= 'viridis', 
                marker = ".", 
                #s=Lmbd
                alpha= 0.5,
                edgecolors= None,
                linewidths = 0,
                )

ax3.set_xlim([Xmin, Xmax])
ax3.set_ylim([Ymin, Ymax])
ax3.set_zlim([Zmin, Zmax])

ax3.set_xlabel('Ru')
ax3.set_ylabel('Ci')
ax3.set_zlabel('Cii')

# orientation
ax3.view_init(azim= Az, elev=elv)
plt.show()

##############################################
############ SPEED EVALUATION ###############
## Uncomment for evaluation in a range near elbow param 
c_i = 0.06 #
c_ii = 0.07 #
ru = np.linspace(0.68, 6, 50)#
vn = np.linspace(0.001, 5, 50)#

# Create vectors
Vn = []
Ru = []
Lmbd = []
Colors = []
for Ivn in vn:
    for Iru in ru:
        Vn.append(Ivn)
        Ru.append(Iru)
        Lmb, col = LambdaVts(Iru, c_i, c_ii, Data, [Ivn, -0.07])
        Lmbd.append(Lmb)
        Colors.append(col)

# create figure
fig3 = plt.figure(figsize=(szx*ctm, szy*ctm))
ax3 = plt.axes(projection='3d')

# plot the cube of points
p = ax3.scatter(Ru, Vn, Lmbd, 
                c = Colors, 
                cmap= 'viridis', 
                marker = "o", 
                alpha= 0.5,
                edgecolors= None,
                linewidths = 0,
                )

Xmin = min(Ru)
Xmax = max(Ru)
Ymin = min(Vn)
Ymax = max(Vn)
Zmin = min(Lmbd)
Zmax = max(Lmbd)

ax3.set_xlim([Xmin, Xmax])
ax3.set_ylim([Ymin, Ymax])
ax3.set_zlim([Zmin, Zmax])

ax3.set_xlabel('Ru')
ax3.set_ylabel('Vn')
ax3.set_zlabel('Lamb')

# orientation
ax3.view_init(azim= Az, elev=elv)
plt.show()

##############################################
############ SPEED EVALUATION ###############
## Uncomment for evaluation in a range near elbow param 
c_i = 0.06 #
c_ii = 0.07 #
# ra = np.linspace(0.68, 6, 100)#
# vn = np.linspace(0.001, 5, 100)#

ra = np.geomspace(0.68, 6, num = 100)
vn = np.geomspace(0.001, 5, num = 100)#

# Create vectors
Vn = []
Ru = []
Lmbd = []
Colors = []
for Ivn in vn:
    for Iru in ra:
        Vn.append(Ivn)
        Ru.append(Iru)
        Lmb, col = LambdaVtsFs(Iru, c_i, c_ii, Data, [Ivn, -0.07])
        Lmbd.append(Lmb)
        Colors.append(col)

# create figure
fig3 = plt.figure(figsize=(szx*ctm, szy*ctm))
ax3 = plt.axes(projection='3d')

# plot the surf of points
p = ax3.scatter(Ru, Vn, Lmbd, 
                c = Colors, 
                cmap= 'viridis', 
                marker = "o", 
                alpha= 0.3,
                edgecolors= None,
                linewidths = 0,
                )

Xmin = min(Ru)
Xmax = max(Ru)
Ymin = min(Vn)
Ymax = max(Vn)
Zmin = min(Lmbd)
Zmax = max(Lmbd)

ax3.set_xlim([Xmin, Xmax])
ax3.set_ylim([Ymin, Ymax])
ax3.set_zlim([Zmin, Zmax])

ax3.set_xlabel('Ru')
ax3.set_ylabel('Vn')
ax3.set_zlabel('Lamb')

# orientation
ax3.view_init(azim= Az, elev=elv)
plt.show()


# v = np.linspace(0, 120, 13)
# t = fig3.colorbar(p, ax = ax3)

# plt.xticks((-1,0,1))
# plt.yticks((-1,0,1))
# ax3.set_zticks((-1,0,1))

# for plne in faces:
#     face = mp3d.art3d.Poly3DCollection([plne], alpha=0.02, linewidth=1)
#     alpha1 = 0.05
#     alpha2 = 1.0
#     face.set_facecolor((0, 0, 0, alpha1))
#     face.set_edgecolor((0, 0, 0, alpha2))
#     ax3.add_collection3d(face)

# solo -1 0 1




# plt.


# # Planes
# xy1 = [ (-1, -1, 0),
#         (-1,  0, 0),
#         ( 0,  0, 0),
#         ( 0, -1, 0),]

# xy2 = [ ( 1, -1, 0),
#         ( 1,  0, 0),
#         ( 0,  0, 0),
#         ( 0, -1, 0),]

# xy3 = [ ( 1,  1, 0),
#         ( 1,  0, 0),
#         ( 0,  0, 0),
#         ( 0,  1, 0),]

# xy4 = [ (-1,  1, 0),
#         (-1,  0, 0),
#         ( 0,  0, 0),
#         ( 0,  1, 0),]

# xy5 = [ (-1, -1, 1),
#         (-1,  0, 1),
#         ( 0,  0, 1),
#         ( 0, -1, 1),]

# xy6 = [ ( 1, -1, 1),
#         ( 1,  0, 1),
#         ( 0,  0, 1),
#         ( 0, -1, 1),]

# xy7 = [ ( 1,  1, 1),
#         ( 1,  0, 1),
#         ( 0,  0, 1),
#         ( 0,  1, 1),]

# xy8 = [ (-1,  1, 1),
#         (-1,  0, 1),
#         ( 0,  0, 1),
#         ( 0,  1, 1),]

# xy9 = [ (-1, -1, -1),
#         (-1,  0, -1),
#         ( 0,  0, -1),
#         ( 0, -1, -1),]

# xy10 = [ ( 1, -1, -1),
#          ( 1,  0, -1),
#          ( 0,  0, -1),
#          ( 0, -1, -1),]

# xy11 = [ ( 1,  1, -1),
#          ( 1,  0, -1),
#          ( 0,  0, -1),
#          ( 0,  1, -1),]

# xy12 = [ (-1,  1, -1),
#          (-1,  0, -1),
#          ( 0,  0, -1),
#          ( 0,  1, -1),]

# yz1 = [ (0, -1, -1),
#         (0, -1,  0),
#         (0,  0,  0),
#         (0,  0, -1),]

# yz2 = [ (0,  1, -1),
#         (0,  1,  0),
#         (0,  0,  0),
#         (0,  0, -1),]

# yz3 = [ (0,  1,  1),
#         (0,  1,  0),
#         (0,  0,  0),
#         (0,  0,  1),]

# yz4 = [ (0, -1,  1),
#         (0, -1,  0),
#         (0,  0,  0),
#         (0,  0,  1),]

# yz5 = [ (1, -1, -1),
#         (1, -1,  0),
#         (1,  0,  0),
#         (1,  0, -1),]

# yz6 = [ (1,  1, -1),
#         (1,  1,  0),
#         (1,  0,  0),
#         (1,  0, -1),]

# yz7 = [ (1,  1,  1),
#         (1,  1,  0),
#         (1,  0,  0),
#         (1,  0,  1),]

# yz8 = [ (1, -1,  1),
#         (1, -1,  0),
#         (1,  0,  0),
#         (1,  0,  1),]

# yz9 = [ (-1, -1, -1),
#         (-1, -1,  0),
#         (-1,  0,  0),
#         (-1,  0, -1),]

# yz10 = [(-1,  1, -1),
#         (-1,  1,  0),
#         (-1,  0,  0),
#         (-1,  0, -1),]

# yz11 = [(-1,  1,  1),
#         (-1,  1,  0),
#         (-1,  0,  0),
#         (-1,  0,  1),]

# yz12 = [(-1, -1,  1),
#         (-1, -1,  0),
#         (-1,  0,  0),
#         (-1,  0,  1),]

# xz1 = [ (-1, 0, -1),
#         (-1, 0,  0),
#         ( 0, 0,  0),
#         ( 0, 0, -1),]

# xz2 = [ ( 1, 0, -1),
#         ( 1, 0,  0),
#         ( 0, 0,  0),
#         ( 0, 0, -1),]

# xz3 = [ ( 1, 0,  1),
#         ( 1, 0,  0),
#         ( 0, 0,  0),
#         ( 0, 0,  1),]

# xz4 = [ (-1, 0,  1),
#         (-1, 0,  0),
#         ( 0, 0,  0),
#         ( 0, 0,  1),]

# xz5 = [ (-1, 1, -1),
#         (-1, 1,  0),
#         ( 0, 1,  0),
#         ( 0, 1, -1),]

# xz6 = [ ( 1, 1, -1),
#         ( 1, 1,  0),
#         ( 0, 1,  0),
#         ( 0, 1, -1),]

# xz7 = [ ( 1, 1,  1),
#         ( 1, 1,  0),
#         ( 0, 1,  0),
#         ( 0, 1,  1),]

# xz8 = [ (-1, 1,  1),
#         (-1, 1,  0),
#         ( 0, 1,  0),
#         ( 0, 1,  1),]

# xz9 = [ (-1, -1, -1),
#         (-1, -1,  0),
#         ( 0, -1,  0),
#         ( 0, -1, -1),]

# xz10 = [( 1, -1, -1),
#         ( 1, -1,  0),
#         ( 0, -1,  0),
#         ( 0, -1, -1),]

# xz11 = [( 1, -1,  1),
#         ( 1, -1,  0),
#         ( 0, -1,  0),
#         ( 0, -1,  1),]

# xz12 = [(-1, -1,  1),
#         (-1, -1,  0),
#         ( 0, -1,  0),
#         ( 0, -1,  1),]

# faces = [xy1,xy2,xy3,xy4,xy5,xy6,xy7,xy8,xy9,xy10,xy11,xy12,yz1,yz2,yz3,yz4,yz5,yz6,yz7,yz8,yz9,yz10,yz11,yz12,xz1,xz2,xz3,xz4,xz5,xz6,xz7,xz8,xz9,xz10,xz11,xz12]



# # Fix ci
# c_i = 0.06 # de Medley
# # cii = np.linspace(0.01, 0.05, 500)
# ra = np.linspace(0.086, 6.5, 500)

# Cii, RA = np.meshgrid(cii, ra)

# Lmbd = LambdaCII(RA, c_i, Cii, Data)

# # plot the surface
# ax2 = fig.add_subplot(132, projection='3d')
# surf = ax2.plot_surface(Cii, RA, Lmbd, cmap='viridis')

# ax2.set_xlabel('cii')
# ax2.set_ylabel('ra')
# ax2.set_zlabel('lmb')

# z_min = np.min(surf.get_array())
# z_max = np.max(surf.get_array())

# print('lmb the min:', z_min)
# print('lmb the max:', z_max)

# # # Fix cii
# c_ii = 5000 # de Medley
# # cii = np.linspace(0.01, 0.05, 500)
# ra = np.linspace(0.086, 6.5, 500)

# Ci, RA = np.meshgrid(ci, ra)

# Lmbd = LambdaCII(RA, Ci, c_ii, Data)

# # plot the surface
# ax2 = fig.add_subplot(133, projection='3d')
# surf = ax2.plot_surface(Ci, RA, Lmbd, cmap='viridis')

# ax2.set_xlabel('ci')
# ax2.set_ylabel('ra')
# ax2.set_zlabel('lmb')

# z_min = np.min(surf.get_array())
# z_max = np.max(surf.get_array())

# print('lmb the min_ci:', z_min)
# print('lmb the max_ci:', z_max)

# fig.subplots_adjust(top=0.985,
# bottom=0.015,
# left=0.0,
# right=0.97,
# hspace=0.195,
# wspace=0.073)