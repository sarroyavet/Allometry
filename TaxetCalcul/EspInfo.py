#Read the CSV file and create the objects of the species with the info
from inspect import trace
import json
from pyexpat import features
import numpy as np
import matplotlib.pyplot as plt
from numpy import NaN
import scipy.optimize as opt
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
from ete3 import NCBITaxa
from ete3 import Tree, TreeStyle
from sklearn.metrics import r2_score
from GetInfoII import *
from Func_Regre import *

# General variables
poptHirt = PopHirt()

ncbi = NCBITaxa()
# ncbi.update_taxonomy_database()

# Get info for the tree
def getNCBIforTree(items, ChrH, ChrV, JS):
    # Get the information that will be displayed in tre graph
    V = [] # vertical axis
    V2 = [] # for calculated data
    H = [] # horizontal axis
    H2 = [] # horizontal axis for calculated data
    NCBI = []
    sci_name = []
    for item in items:
        if (JS[item][ChrV] != "NaN"):
            if (JS[item][ChrH] != "NaN"):
                V.append(float(JS[item][ChrV]))
                H.append(float(JS[item][ChrH]))
                NCBI.append(item)
                sci_name.append(JS[item]['Sci_name'])
        if ChrV == "Maxspeed":
            if (JS[item][ChrV] != "NaN"):
                if (JS[item][ChrH] != "NaN"):
                    V2.append(float(JS[item][ChrV]))
                    H2.append(float(JS[item][ChrH]))
                    # NCBI2.append(item)
                    # sci_name2.append(JS[item]['Sci_name'])
            else:
                if (JS[item][ChrH] != "NaN"):
                    V2.append(float(JS[item]['speedHirt']))
                    H2.append(float(JS[item][ChrH]))
        if ChrH == "Maxspeed":
            if (JS[item][ChrH] != "NaN"):
                if (JS[item][ChrV] != "NaN"):
                    V2.append(float(JS[item][ChrV]))
                    H2.append(float(JS[item][ChrH]))
                    # NCBI2.append(item)
                    # sci_name2.append(JS[item]['Sci_name'])
            else:
                if (JS[item][ChrH] != "NaN"):
                    V2.append(float(JS[item][ChrV]))
                    H2.append(float(JS[item]['speedHirt']))
    return(NCBI, sci_name, H, V, H2, V2)

# Plot the tree
def PlotTree(NCBI, filesBN, Plot_tree):
    if Plot_tree == True:
        tree = ncbi.get_topology(NCBI, intermediate_nodes=False, )
        for leaf in tree:
            if leaf.name in filesBN:
                leaf.add_features(filesID = ', '.join(filesBN[leaf.name]), files = str(len(filesBN[leaf.name])))
            else:
                leaf.add_features(filesID = " ", files = '0')
        circular_style = TreeStyle()
        circular_style.mode = "c" # draw tree in circular mode
        tree.show(tree_style=circular_style,)

####################################################################
fileDir = '/home/kale/Documents/Allometry/DATA/'
espfile = fileDir+'New tax.csv'
Spout = fileDir+'Esp_out'

D3fileDB = fileDir+'3DfilesDb.csv'
Bonesout = fileDir+'Bone_out'

Fls = fileDir+'fls.json'

Get_Info(D3fileDB, espfile, Bonesout, Spout, Fls, fileDir, New=False)
# GetInfoSpecies(espfile, espjason, espcsv)
# GetInfoBones(D3fileDB, Bonesjson, Bonescsv, fileDir)
# UpdateSpecies(D3fileDB, espjason, Fls, espcsv)
# graphloglog("Maxspeed","massAvg", espjason, Fls,)

# # Linear regression
# def linRegression(H,V, func):
#     popt, pcov = opt.curve_fit(func, H, V)
#     perr = np.sqrt(np.diag(pcov))
#     residuals = V- func(H, *popt)
#     ss_res = np.sum(residuals**2)
#     ss_tot = np.sum((V-np.mean(V))**2)
#     r_squared = 1 - (ss_res / ss_tot)
#     y_pred = func(H, *popt)
#     r2 = r2_score(V, y_pred)
#     return(popt, pcov, r_squared)


# def graphloglog(ChrV, ChrH, jsonfile, Fls, Plot_tree = False, Mass=None, stance=None, tax=None, speed =None, stridef=None, PCSA= None, Diet = None, Colorclf = None):
#     #################################################################
#     #   Mass = [liminf, limsup]
#     #   stance = [list of the wanted stance]
#     #   tax = [list of the wanted taxes]
#     #   speed = [liminf, limsup]
#     #   stride = [liminf, limsup]
#     #   PCSA = [liminf, limsup]
#     #   Diet = [list of the wanted food habits]
#     #   Colorclf = [list of color according to what]

#     C = [] # color

#     # NCBI codes to ID codes
#     with open(Fls,'r') as outfile:
#         filesBN = json.load(outfile)
    
#     # Open the info of the species and get specific data
#     with open(jsonfile,'r') as outfile:
#         JS = json.load(outfile)
#         items = [item for item in JS]
#         for item in JS:
#             if Mass != None:
#                 if Mass[0] >= float(JS[item]["massAvg"]) >= Mass[-1]:
#                     items.remove(item)
#             if stance != None:
#                 if JS[item]["Stance"] not in stance:
#                     items.remove(item)
#             if tax != None:
#                 if any([itt in JS[item]["lineage_names"] for itt in tax])== False:
#                     items.remove(item)
#             if speed != None:
#                 if speed[0] >= float(JS[item]["Maxspeed"]) >= speed[-1]:
#                         items.remove(item)
#             if stridef != None:
#                 if stridef[0] >= float(JS[item]["StrideFreq"]) >= stridef[-1]:
#                     items.remove(item)
#             if PCSA != None:
#                 if PCSA[0] >= float(JS[item]["PCSA"]) >= PCSA[-1]:
#                     items.remove(item)
#             if Diet != None:
#                 if JS[item]["Diet"] not in Diet:
#                     items.remove(item)
        
#         NCBI, sci_name,  H, V, H2, V2 = getNCBIforTree(items, ChrH, ChrV, JS)
    
#     PlotTree(NCBI, filesBN, Plot_tree)

#     # LOG FUNCTION
#     # optimization function
#     popt, pcov, perr = linRegression(H,V, regfunct)
#     print('normal', perr)

#     # Plot of the function
#     newX = np.logspace(-2, 4, base=10)
#     fig, (ax1) = plt.subplots(1,1)
#     ax1.scatter(H,V,s=60, alpha=0.7, label='data')
#     ax1.plot(newX, regfunct(newX, *popt), 'r-', label='('+str(popt[0])+'*x^'+ str(popt[1])+') r^2 =' + str(perr))

#     # Plot extras
#     if ChrH or ChrV == "Maxspeed":
#         poptH, pcovH, perrH = linRegression(H2, V2, regfunct)
#         print(poptH)
#         print('con extra data', perrH)
#         ax1.plot(newX, regfunct(newX, *poptH), label='('+str(poptH[0])+'*x^'+ str(poptH[1])+') r^2 =' + str(perrH))
#         ax1.plot(newX, SpeedHirt(newX, *poptHirt), label="Hirt et al.")
    
#     # plot details
#     ax1.set_yscale("log")
#     ax1.set_xscale("log")
#     ax1.grid(b='on')
#     plt.legend(loc='lower right')
#     plt.show()

    
#     # # To see details in Potly
#     # df =  pd.DataFrame({'NCBI': NCBI, ChrH: H, ChrV: V, 'Sci_name': sci_name})
#     # fig = px.scatter(df,
#                     #  x = df[ChrH], 
#                     #  y = df[ChrV],
#                     #  hover_name= sci_name)
#     # # fig.add_trace(px.line())
#     # fig.update_xaxes(type="log")
#     # fig.update_yaxes(type="log")
#     # fig.show()