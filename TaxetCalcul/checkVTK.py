import json
from Get_MessIR import read_mesh
from Get_MessIR import Plot3Dscatter
import os

fileDir = '/home/kale/Documents/Allometry/DATA/'
jsonfile = '/home/kale/Documents/Allometry/DATA/fls.json'
with open(jsonfile,'r') as outfile:
        JS = json.load(outfile)
        items = [item for item in JS]

        for esp in items:
            for bn in JS[esp]:
                file = fileDir+'Distal/'+bn+'.vtk'
                if os.path.exists(file):
                    nbp, X, Y, Z = read_mesh(file,1)
                    print(bn)
                    Plot3Dscatter(bn,X,Y,Z)
