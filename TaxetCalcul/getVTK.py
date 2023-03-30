################################################################################################
''' Developed by Kalenia Marquez Florez
    Aix Marseille Univ, CNRS, ISM, Marseille, France
    March 2023

'''
#
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
            f.readline()
            if 'POLYGONS' in line:
                break
            else:
                line2 = line.split()
                npts = int(len(line2)/3 )
                for i in range(npts):
                    X.append(line2[i*3])
                    Y.append(line2[i*3+1])
                    Z.append(line2[i*3+2])
                line = f.readline()
    
    return(NoPoints, X, Y, Z)




