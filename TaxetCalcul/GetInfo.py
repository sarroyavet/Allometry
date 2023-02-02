from cmath import e
import csv
import json
from ete3 import NCBITaxa
from Get_MessIR import ReadFile
from Func_Regre import *
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import math

poptHirt = PopHirt()

ncbi = NCBITaxa()
# ncbi.update_taxonomy_database()

# Species
class Specie:
    def __init__(self, NCBI, massI, massS, massAvg, refmass, Maxspeed, refspeed, spdMissing, Stance, StrideFreq, refStride, PCSA, refPCSA, Diet):
        self.NCBI = NCBI
        self.massI = massI
        self.massS = massS
        self.massAvg = massAvg
        self.refmass = refmass
        self.Maxspeed = Maxspeed
        self.refspeed = refspeed
        self.speedHirt = SpeedHirt(self.massAvg,*poptHirt)
        self.spdMissing = spdMissing
        self.Stance = Stance
        self.StrideFreq = StrideFreq
        self.refStride = refStride
        self.PCSA = PCSA
        self.refPCSA = refPCSA
        self.Diet = Diet
    def bones(self):
        self.D3files = {}
        self.nofiles = 0
    def taxa(self):
        self.Sci_name = ncbi.get_taxid_translator([self.NCBI])[self.NCBI]
        try:
            self.common_name = ncbi.get_common_names([self.NCBI])[self.NCBI]
        except:
            self.common_name = 'NaN'
        self.lineage = ncbi.get_lineage(self.NCBI)
        self.lineage_names = [ncbi.get_taxid_translator([id])[id] for id in self.lineage]

# bones
class bone:
    def __init__(self, ID, NCBI, side, Long): # latter more properties will be added
        self.ID = ID
        self.NCBI = NCBI
        self.side = side
        self.Long = Long
    def taxa(self):
        self.Sci_name = ncbi.get_taxid_translator([self.NCBI])[self.NCBI]
        try:
            self.common_name = ncbi.get_common_names([self.NCBI])[self.NCBI]
        except:
            self.common_name = 'NaN'
        self.lineage = ncbi.get_lineage(self.NCBI)
        self.lineage_names = [ncbi.get_taxid_translator([id])[id] for id in self.lineage]
    def DataJoint(self, file_dir, scale):
        self.Rmax, self.Rmin, self.Rmed, self.La = MeasureBones(self.ID, file_dir,scale)
    def mass(self):
        self.massCalcul = CalcMass(float(self.Long), self.lineage_names)
        if self.massCalcul > 0e-10:
            self.speedHirt = SpeedHirt(self.massCalcul,*poptHirt)
        else:
            self.speedHirt = float('nan')
    def speed(self):
        self.Max_speed = speed(self.NCBI)

def CalcMass(Long, lineage_names):
    if 'Carnivora' in lineage_names: #christiansen
        a = 56.157
        b = 0.3451
    elif 'Artiodactyla' in lineage_names:
        a = 53.315
        b = 0.3045
    else: 
        a = 56.969
        b = 0.3109
    if Long > 10e-12:
        return((Long/a)**(1/b))
    else:
        return(float("nan"))

# Measure bones
def MeasureBones(ID, file_dir, scale):
    fileDB = file_dir+'Distal/'+ID+'.vtk'
    print(ID)
    SF, Exists = ReadFile(ID,fileDB,scale)
    if Exists:
        return(SF.Rmax, SF.Rmin, SF.Rmed, SF.L)
    else: 
        return(float("nan"), float("nan"), float("nan"), float("nan"))

# Read bones csv file make the object and the json file 
def GetInfoBones(D3fileDB, Bonesjson, Bonescsv, fileDir):
    HuesosIDs = {}
    with open(D3fileDB, 'r') as csvfile: 
        case = csv.DictReader(csvfile)
        for item in case:
            ID = item["ID"]
            NCBI = int(item["NCBI"])
            side = item["side"]
            long = float(item["Long(mm)"])*float(item["Mult"])
            scale = float(item["Mult"])
            bn = bone(ID,NCBI, side, long)
            bn.taxa()
            # bn.DataJoint(fileDir, scale)
            bn.mass()
            HuesosIDs.update({str(item["ID"]) : bn.__dict__})
    
    GetDF(HuesosIDs, Bonescsv)

    with open(Bonesjson,'w') as outfile:
        json.dump(HuesosIDs, outfile,indent=4)

# Read species csv file make the object and the json file
def GetInfoSpecies(espfile, espjason, espcsv):
    Especies = {}
    EspeciesClss = {}
    with open(espfile, 'r') as csvfile: 
        case = csv.DictReader(csvfile)
        for item in case:
            NCBI = int(item["NCBI"])
            if item["mass inf kg"]=='':
                massI = float(item["Avg mass kg"])
                massS = float(item["Avg mass kg"])
                massAvg = float(item["Avg mass kg"])
            else:
                massI = float(item["mass inf kg"])
                massS = float(item["mass sup kg"])
                massAvg = (massI+massS)/2.0
            refmass = item["refMass"]
            if item["Speed max km/hr"] == '':
                Maxspeed = float("nan")
                refspeed = float("nan")
                spdMissing = SpeedHirt(massAvg,*poptHirt)
            else:
                Maxspeed = float(item["Speed max km/hr"])
                refspeed = item["refSpeed"]
                spdMissing = float(item["Speed max km/hr"])
            Stance = item["stance"]
            if item["stride freq"] == '':
                StrideFreq = ''
                refStride  = ''
            else:
                StrideFreq = float(item["stride freq"])
                refStride  = item["refStride"]
            if item["physiological cross-sectional area (PCSA)"] == '':
                PCSA = ''
                refPCSA = ''
            else:
                PCSA = float(item["physiological cross-sectional area (PCSA)"])
                refPCSA = item["refPCSA"]
            Diet = item["Diet"]

            sp = Specie(NCBI, massI, massS, massAvg, refmass, Maxspeed, refspeed, spdMissing, Stance, StrideFreq, refStride, PCSA, refPCSA, Diet)
            sp.bones()
            sp.taxa()

            Especies.update({str(item["NCBI"]) : sp.__dict__})
            EspeciesClss.update({str(item["NCBI"]) : sp})
    
    with open(espjason,'w') as outfile:
        json.dump(Especies, outfile,indent=4)

# Update dictionary of species with the bones    
def UpdateSpecies(D3fileDB, jsonfile, Fls, espcsv):
    with open(jsonfile,'r') as outfile:
        Especies = json.load(outfile)
    
    Filesdic = {}
    with open(D3fileDB, 'r') as csvfile: 
        case = csv.DictReader(csvfile)
        for item in case:
            NCBI = item["NCBI"]
            if NCBI in Especies:
                Especies[item["NCBI"]]["D3files"].update({item["ID"]: {"side": item["side"]}})
                Especies[item["NCBI"]]["nofiles"] = len(Especies[item["NCBI"]]["D3files"])
                if NCBI in Filesdic:
                    Filesdic[NCBI].append(item["ID"])
                else:
                    Filesdic.update({NCBI: [item["ID"]]})
  
    GetDF(Especies, espcsv)

    with open(jsonfile,'w') as outfile:
        json.dump(Especies, outfile,indent=4)

    with open(Fls,'w') as outfile:
        json.dump(Filesdic, outfile,indent=4)

# # get the averages 
# def Measure_Average():




# get df in pandas 
def GetDF(infodict, csvfile):
    dfinf = pd.DataFrame.from_dict(infodict, orient='index') 
    dfinf.to_csv(csvfile+'.csv')
    dfinf.to_pickle(csvfile+'.pkl')