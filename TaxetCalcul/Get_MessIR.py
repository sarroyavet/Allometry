# import meshio
from cmath import pi
from vtk import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from cylinder_fitting import *
import numpy as np
from scipy import interpolate
# from scipy import find_peaks
import os
from matplotlib.pyplot import cm
import peakdetect

delta_theta = 90 # number of angles
knots = 2
imgpath = '/home/kale/Documents/Allometry/DATA/DistalImg'

# class surface: for each distal humerus surface
class Surface():
    def __init__(self, ID, file, scale, morph_type):
        self.ID = ID
        self.file = file
        self.scale = scale
        self.morph_type = morph_type
    def GetPoints(self):
        self.nbp, self.X, self.Y, self.Z = read_mesh(self.file, self.scale)
    def Get_cyl(self):
        self.w_fit, self.C_fit, self.r_fit, self.RMSD = Axis_Fitting(self.X, self.Y, self.Z)
    def PlotSur(self):
        Plot3Dscatter(self.ID, self.X, self.Y, self.Z)
    def Origin(self):
        self.X_W, self.Y_W, self.Z_W = GoToOrigin(self.X, self.Y, self.Z, self.w_fit, self.C_fit, self.r_fit)
    def MessRL(self):
        self.Rmax, self.Rmin, self.Rmed, self.L = Get_measurements(self.X_W, self.Y_W, self.Z_W, np.array([1,0,0]))
    def Get_profiles(self):
        self.profiles = get_Prof_points(self.X_W, self.Y_W,self.Z_W)
    def print_profiles(self):
        self.xfitl, self.yfitl = print_prof2(self.profiles, self.morph_type, self.ID)
        # self.xfitl, self.yfitl, self.pts, self.slopes = print_prof(self.profiles, self.morph_type, self.ID)


# Read the mesh in .vtk any of the two functions below works fine
def read_mesh(file, scale):
    mesh = vtkDataSetReader()
    mesh.SetFileName(file)
    # mesh.ReadAllVectorsOn()
    # mesh.ReadAllScalarsOn()
    mesh.Update()
    data = mesh.GetOutput()
    X = [0]*data.GetNumberOfPoints()
    Y = [0]*data.GetNumberOfPoints()
    Z = [0]*data.GetNumberOfPoints()

    for i in range(data.GetNumberOfPoints()):
        x,y, z = data.GetPoint(i)
        X[i] = x*scale
        Y[i] = y*scale
        Z[i] = z*scale

    return (data.GetNumberOfPoints(),X,Y,Z)

def openVTK(vtk_file):
    f = open(vtk_file, 'r')
    line = f.readline()
    X = []
    Y = []
    Z = []
    while line:
        if 'POINTS' in line:
            line2 = line.split()
            NoPoints = int(line2[1])
            print(NoPoints)
            for l in range(int(NoPoints/3)):
                line = f.readline()
                if 'POLYGONS' in line:
                    print('aqui')
                    break
                else:
                    line2 = line.split()
                    npts = int(len(line2)/3 )
                    for i in range(npts):
                        X.append(float(line2[i*3]))
                        Y.append(float(line2[i*3+1]))
                        Z.append(float(line2[i*3+2]))
        else:
            line = f.readline()
    return(NoPoints, X, Y, Z)

# Plot the surface points in 3D
def Plot3Dscatter(ID,X,Y,Z, w_fit=None, C_fit=None, r_fit=None, plot_cylinder = False):
    plt.ion
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot(X,Y,Z, marker='o', ls="", markersize = 2)
    plt.title(ID)

    if np.any(w_fit != None):
        ax.quiver(C_fit[0], C_fit[1], C_fit[2], r_fit*w_fit[0], r_fit*w_fit[1], r_fit*w_fit[2], color = 'red')

    # plt.show()

# Fit a cylinder to the points
def Axis_Fitting(X, Y, Z):
    pdata = []
    for i in range(len(X)):
        pdata.append(np.array([X[i],Y[i],Z[i]]))

    w_fit, C_fit, r_fit, fit_err = fit(pdata)

    RMSD = fitting_rmsd(w_fit, C_fit, r_fit, pdata) # root-mean-square deviation (RMSD)

    # show_fit(w_fit, C_fit, r_fit, pdata)

    # show_G_distribution(pdata)
    # print(w_fit, C_fit, r_fit, fit_err, RMSD)

    return(w_fit, C_fit, r_fit, RMSD)

def get_distance_to_axis(X, Y, Z, w):
    D =[]
    for i in range(len(X)):
        D.append((np.linalg.norm(np.cross(np.array([X[i], Y[i], Z[i]]), w))))
    return(D)

# Transform coordinates
def GoToOrigin(X, Y, Z, w_fit, C_fit, r_fit):
    C_new = [0,0,0]
    desplz = C_new - C_fit

    X_new = X+desplz[0]
    Y_new = Y+desplz[1]
    Z_new = Z+desplz[2]

    D = get_distance_to_axis(X_new, Y_new, Z_new, w_fit)

    # Further point P
    P = np.array([X_new[D.index(min(D))], Y_new[D.index(min(D))], Z_new[D.index(min(D))]])
    Pw = w_fit*np.dot(P,w_fit) # in the direction of the axis
    Pv = P - Pw # in the direction of perpendicular to the axis

    v_fit = Pv/np.linalg.norm(Pv)
    q_fit = np.cross(w_fit,v_fit)/np.linalg.norm(np.cross(w_fit,v_fit))

    # Rotation matrix
    RotMat = np.linalg.inv(np.transpose(np.array([w_fit,v_fit,q_fit])))

    # Point rotation
    X_W = [0]*len(X_new)
    Y_W = [0]*len(X_new)
    Z_W = [0]*len(X_new)
    for i in range(len(X_new)):
        CC = np.matmul(RotMat, np.array([X_new[i], Y_new[i], Z_new[i]]))
        X_W[i] = CC[0]
        Y_W[i] = CC[1]
        Z_W[i] = CC[2]

    # Move in X
    desplz_X = 0-min(X_W)
    X_WW = X_W + desplz_X
    C_g = np.array([desplz_X, 0, 0]) # to visualize

    # Visualization
    pdata = []
    for i in range(len(X)):
        pdata.append(np.array([X_WW[i],Y_W[i],Z_W[i]]))

    #### SHOW FIT
    # show_fit(np.array([1,0,0]), C_g, r_fit, pdata)

    return(X_WW, Y_W, Z_W)

def Get_measurements(X, Y, Z, W):
    D = get_distance_to_axis(X, Y, Z, W)
    Rmax = max(D)
    Rmin = min(D)
    Rmed = sum(D)/len(D)
    L = max(X)-min(X)

    return(Rmax, Rmin, Rmed, L)

def get_Prof_points(X, Y, Z):
    prof = {}
    angles = np.linspace(0, 2*np.pi, num = delta_theta)
    for th in angles:
        Xs = []
        Ys = []
        Zs = []
        Rs = []
        for i in range(len(X)):
            rho = (Z[i]**2+Y[i]**2)**0.5
            theta = np.arctan2(Y[i],Z[i]) # rad
            if abs(theta-th) <= 10e-2:
                Xs.append(X[i])
                Ys.append(Y[i])
                Zs.append(Z[i])
                Rs.append(rho)
        points = {'Xs': Xs}
        points.update({'Ys': Ys})
        points.update({'Zs': Zs})
        points.update({'Rs': Rs})
        if Rs:
            dd = max(Xs)-min(Xs)
        else:
            dd = 0.0
        prof.update({'th%5.2f'%th: {'points': points, 'dist': dd, 'angle': th}})
    return(prof)

# def print_prof(profile_points, morph_type, ID):
    fig, axs = plt.subplots(4, figsize=(30, 20))
    axs[0].axis('equal')
    axs[1].axis('equal')
    axs[2].axis('equal')
    axs[3].axis('equal')
    color = iter(cm.nipy_spectral(np.linspace(0, 1, delta_theta)))
    ids = []
    clr = []
    for angle in profile_points:
        rho = profile_points[angle]['points']['Rs']
        Xs =  profile_points[angle]['points']['Xs']
        if rho:
            if len(rho) > 5:
                Xs, rho = (list(t) for t in zip(*sorted(zip(Xs, rho))))
                yfit, xfit, _, _, _, _ = spline(Xs, rho)
                c = next(color)
                # axs[0].plot(Xs, rho,'o', c=c) # to plot points
                axs[0].plot(xfit, yfit, alpha = 0.5, c=c)
                clr.append(c)
                ids.append(profile_points[angle]['angle'])
                Xfitx, Yfitx = min_x(xfit.tolist(), yfit.tolist())
                axs[1].plot(Xfitx, Yfitx, alpha = 0.5, c = c)

    # get the longest opposed to where there is 0 points
    deltatheta= abs(max(ids) - min(ids))
    angleav = (max(ids) - min(ids))/2
    if deltatheta < np.pi:
        Angle = angleav+np.pi
        if Angle >= 2*np.pi:
            Angle = Angle-2*np.pi
    else:
        Angle = angleav
    
    #angle for profile
    ang = min(ids, key=lambda x:abs(x-Angle))
    # print(ang)
    rhol = profile_points['th%5.2f'%ang]['points']['Rs']
    Xsl =  profile_points['th%5.2f'%ang]['points']['Xs']
    Xsl, rhol = (list(t) for t in zip(*sorted(zip(Xsl, rhol))))
    yfitl, xfitl, pts, slopes, yder_3, yder_4 = spline(Xsl, rhol)
    axs[2].plot(xfitl, yder_3, c='black')
    axs[2].plot(xfitl, yder_4, c='red')
    axs[3].plot(xfitl, yfitl, alpha=0.5, c=clr[ids.index(ang)])
    if len(slopes)>morph_type*2:
        remove = len(slopes)-(morph_type*2)
        U = [abs(x) for x in slopes]
        for r in range(remove):
            pts.pop(U.index(min(U)))
            slopes.pop(U.index(min(U)))
            U.pop(U.index(min(U)))
    tt = 0
    for pt in pts:
        plt.axline(xy1= pt, slope=slopes[tt], c='grey', alpha=0.5)
        tt+=1
    # plt.show()
    xfitl = xfitl.tolist()
    yfitl = yfitl.tolist()
    fig.savefig(imgpath+'/'+ID+'.svg', format='svg')
    plt.clf()
    return(xfitl, yfitl, pts, slopes)

def print_prof2(profile_points, morph_type, ID):
    fig, axs = plt.subplots(2,figsize=(30, 20))
    axs[0].axis('equal')
    axs[1].axis('equal')
    color = iter(cm.nipy_spectral(np.linspace(0, 1, delta_theta)))
    ids = []
    clr = []
    for angle in profile_points:
        rho = profile_points[angle]['points']['Rs']
        Xs =  profile_points[angle]['points']['Xs']
        if rho:
            if len(rho) > 5:
                Xs, rho = (list(t) for t in zip(*sorted(zip(Xs, rho))))
                yfit, xfit, _, _, _, _ = spline(Xs, rho)
                c = next(color)
                axs[0].plot(xfit, yfit, alpha = 0.5, c=c)
                clr.append(c)
                ids.append(profile_points[angle]['angle'])

    # get the longest opposed to where there is 0 points
    deltatheta= abs(max(ids) - min(ids))
    angleav = (max(ids) - min(ids))/2
    if deltatheta < np.pi:
        Angle = angleav+np.pi
        if Angle >= 2*np.pi:
            Angle = Angle-2*np.pi
    else:
        Angle = angleav
    
    #angle for profile
    ang = min(ids, key=lambda x:abs(x-Angle))
    # print(ang)
    rhol = profile_points['th%5.2f'%ang]['points']['Rs']
    Xsl =  profile_points['th%5.2f'%ang]['points']['Xs']
    Xsl, rhol = (list(t) for t in zip(*sorted(zip(Xsl, rhol))))
    yfitl, xfitl, pts, slopes, yder_3, yder_4 = spline(Xsl, rhol)
    axs[1].plot(xfitl, yfitl, alpha=0.5, c=clr[ids.index(ang)])
    # plt.show()
    xfitl = xfitl.tolist()
    yfitl = yfitl.tolist()
    fig.savefig(imgpath+'/'+ID+'.svg', format='svg')
    # plt.cla()
    plt.close('all')
    return(xfitl, yfitl)

def scale_profile(x,y):
    Xmax = max(x)
    Xmin = min(x)
    dd = Xmax-Xmin
    xs = []
    ys = []
    for xx in x:
        xs.append((xx-Xmin)/dd)
    for rr in y:
        ys.append(rr/dd)
    return(xs, ys, dd)

def min_x(x,y):
    ymin = min(y)
    xmin = x[y.index(ymin)]
    xs = []
    for xx in x:
        xs.append(xx-xmin)
    return(xs, y)

def spline(x, y):
    # Original
    s=len(x)-np.sqrt(2*len(x))
    t, c, kk = interpolate.splrep(x, y, k = 5, s=s)
    xnew = np.linspace(min(x), max(x))
    yfit = interpolate.BSpline(t,c, kk)(xnew)
    # First derivate
    s=len(xnew)-np.sqrt(2*len(xnew))
    yder_4 = interpolate.splev(xnew, (t,c, kk), der=1)
    t4, c4, kk4 = interpolate.splrep(xnew, yder_4, k = 4, s=s)
    # Second derivate
    yder_3 = interpolate.splev(xnew, (t,c, kk), der=2)
    t3, c3, kk3 = interpolate.splrep(xnew, yder_3, k = 3, s=s)
    #Roots
    Xroots = interpolate.sproot((t3, c3, kk3),)
    pts = []
    slopes = []
    for xx in Xroots:
        yy = interpolate.BSpline(t,c, kk)(xx)
        yyy = yy.item(0)
        pts.append((xx,yyy))
        slp = interpolate.BSpline(t4, c4, kk4)(xx)
        slpp = slp.item(0)
        slopes.append(slpp)
    return yfit, xnew, pts, slopes, yder_3, yder_4

def ReadFile(ID, filename, scale, morph_type):
    sf = Surface(ID, filename, scale, morph_type)
    if os.path.exists(filename):
        sf.GetPoints()
        sf.PlotSur()
        # pointsok = input('Lo hago o paila?\n')
        pointsok = 's'
        if pointsok == 's':
            sf.Get_cyl()
            sf.Origin()
            sf.MessRL()
            sf.Get_profiles()
            sf.print_profiles()
            Exists = True
        else:
            Exists = False
    else:
        Exists = False
    return(sf, Exists)