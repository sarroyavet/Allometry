from cmath import e
import csv
import json
from ete3 import NCBITaxa
from Get_MessIR import ReadFile
from Func_Regre import *
import pandas as pd
import math

# for speed based on Hirt equation: 
# A general scaling law reveals why the largest animals are not the fastest 2017 
poptHirt = PopHirt()

# To get the names order and families
ncbi = NCBITaxa()
# ncbi.update_taxonomy_database()

# Mammals Orders:
Orders= ['Monotremata', 'Didelphimorphia', 'Microbiotheria', 'Paucituberculata', 'Dasyuromorphia', 'Peramelemorphia', 'Notoryctemorphia', 'Diprotodontia', 'Afrosoricida', 'Macroscelidea', 'Tubulidentata', 'Hyracoidea', 'Proboscidea', 'Sirenia', 'Cingulata', 'Pilosa', 'Scandentia', 'Dermoptera', 'Primates', 'Lagomorpha', 'Erinaceomorpha', 'Soricomorpha', 'Chiroptera', 'Pholidota', 'Carnivora', 'Perissodactyla', 'Artiodactyla', 'Cetacea', 'Rodentia']

# Mammals Families:
Families=['Ornithorhynchidae', 'Tachyglossidae', 'Didelphidae', 'Microbiotheria', 'Paucituberculata', 'Dasyuridae', 'Myrmecobiidae', 'Thylacinidae', 'Malleodectidae', 'Peramelidae', 'Chaeropodidae', 'Peroryctidae', 'Notoryctidae', 'Phascolarctidae', 'Vombatidae', 'Phalangeridae', 'Potoroidae', 'Macropodidae', 'Burramyidae', 'Pseudocheiridae', 'Petauridae', 'Tarsipedidae', 'Acrobatidae', 'Hypsiprymnodontidae', 'Chrysochloridae', 'Tenrecidae', 'Macroscelididae', 'Orycteropodidae', 'Procaviidae', 'Pliohyracidae', 'Elephantidae', 'Dugongidae', 'Trichechidae', 'Prorastomidae', 'Chlamyphoridae', 'Dasypodidae', 'Cyclopedidae', 'Myrmecophagidae', 'Choloepodidae', 'Bradypodidae', 'Megalonychidae',  'Tupaiidae', 'Ptilocercidae',  'Cynocephalidae', 'Cheirogaleidae', 'Lemuridae', 
'Lepilemuridae', 'Indriidae', 'Daubentoniidae' , 'Lorisidae', 'Galagidae', 'Tarsiidae','Cebidae', 'Aotidae', 'Pitheciidae', 'Atelidae', 'Cercopithecidae', 'Hylobatidae', 'Hominidae', 'Leporidae', 'Ochotonidae', 'Prolagidae', 'Erinaceidae', 'Solenodontidae', 'Soricidae', 'Talpidae', 
'Vespertilionidae', 'Thyropteridae', 'Rhinopomatidae', 'Rhinolophidae', 'Pteropodidae', 'Phyllostomidae', 'Nycteridae', 'Noctilionidae', 'Natalidae', 'Myzopodidae', 'Mystacinidae', 'Mormoopidae', 'Molossidae', 'Megadermatidae', 'Hipposideridae', 'Furipteridae', 'Emballonuridae', 'Craseonycteridae', 'Manidae', 'Felidae', 'Viverridae', 'Eupleridae', 'Nandiniidae', 'Herpestidae', 'Hyaenidae', 'Canidae', 'Ursidae', 'Odobenidae' , 'Otariidae', 'Phocidae', 'Mustelidae', 'Procyonidae', 'Ailuridae', 'Equidae', 'Tapiridae', 'Rhinocerotidae', 'Suidae', 'Tayassuidae', 'Hippopotamidae', 'Camelidae', 'Tragulidae', 'Moschidae', 'Cervidae', 'Antilocapridae', 'Giraffidae', 'Bovidae', 'Balaenidae', 'Balaenopteridae', 'Eschrichtiidae', 'Cetotheriidae', 'Delphinidae', 'Monodontidae', 'Phocoenidae', 'Physeteridae', 'Platanistidae', 'Iniidae', 'Ziphiidae', 'Aplodontiidae', 'Sciuridae', 'Gliridae', 'Castoridae', 'Heteromyidae', 'Geomyidae', 'Dipodidae', 'Platacanthomyidae', 'Spalacidae', 'Calomyscidae', 'Nesomyidae', 'Cricetidae', 'Muridae', 'Anomaluridae', 'Pedetidae', 'Ctenodactylidae', 'Bathyergidae', 'Hystricidae', 'Petromuridae', 'Thryonomyidae', 'Erethizontidae', 'Chinchillidae', 'Dinomyidae', 'Caviidae', 'Dasyproctidae', 'Cuniculidae', 'Ctenomyidae', 'Octodontidae', 'Abrocomidae', 'Echimyidae', 'Myocastoridae', 'Capromyidae', 'Heptaxodontidae', 'Hydrochaeridae']

# BONES
# class bones
class bone:
    def __init__(self, ID, NCBI, side, Long, morph_type, Morphotype): # latter more properties will be added
        self.ID = ID
        self.NCBI = NCBI
        self.side = side
        self.Long = Long
        self.morph_type = morph_type
        self.Morphotype = Morphotype
        self.Max_speed = float('nan')
    def taxa(self):
        self.Sci_name = ncbi.get_taxid_translator([self.NCBI])[self.NCBI]
        try:
            self.common_name = ncbi.get_common_names([self.NCBI])[self.NCBI]
        except:
            self.common_name = 'NaN'
        self.lineage = ncbi.get_lineage(self.NCBI)
        self.lineage_names = [ncbi.get_taxid_translator([id])[id] for id in self.lineage]
        self.order = get_order(self.lineage_names)
    def DataJoint(self, file_dir, scale):
        # self.Rmax, self.Rmin, self.Rmed, self.La, self.RMSD, self.xfitl, self.yfitl, self.pts, self.slopes = MeasureBones(self.ID, file_dir,scale,self.morph_type)
        self.Rmax, self.Rmin, self.Rmed, self.La, self.RMSD, self.xfitl, self.yfitl = MeasureBones(self.ID, file_dir,scale,self.morph_type)
    def mass(self):
        self.massCalcul = CalcMass(float(self.Long), self.lineage_names)
        if self.massCalcul > 0e-10:
            self.speedHirt = SpeedHirt(self.massCalcul,*poptHirt)
        else:
            self.speedHirt = float('nan')

# BONES SPECIES
# Class species 
class Specie:
    def __init__(self, NCBI, massI, massS, massAvg, refmass, Maxspeed, refspeed, spdMissing, Stance, Diet):
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
        self.Diet = Diet
    def bones(self, ListObj, Fls):
        self.D3files, self.nofiles, self.Long, self.Rmax, self.Rmin, self.Rmed, self.La, self.LrL, self.LrLa, self.LrRmin, self.LrRmax, self.LrRmed, self.Morphotype = GetBones(self.NCBI, ListObj, Fls)
        # self.D3files, self.nofiles, self.Long, self.Rmax, self.Rmin, self.Rmed, self.La, self.LrL, self.LrLa, self.LrRmin, self.LrRmax, self.LrRmed, self.Lrpts, self.Lrslopes, self.Morphotype = GetBones(self.NCBI, ListObj, Fls)
    def taxa(self):
        self.Sci_name = ncbi.get_taxid_translator([self.NCBI])[self.NCBI]
        try:
            self.common_name = ncbi.get_common_names([self.NCBI])[self.NCBI]
        except:
            self.common_name = 'Na'
        self.lineage = ncbi.get_lineage(self.NCBI)
        self.lineage_names = [ncbi.get_taxid_translator([id])[id] for id in self.lineage]
        self.order = get_order(self.lineage_names)
        self.family = get_family(self.lineage_names)
    def massChristiansen(self):
        self.mass_Ch = CalcMass(self.Long, self.lineage_names)
        if self.mass_Ch > 0e-10:
            self.speedHirtChris = SpeedHirt(self.mass_Ch,*poptHirt)
        else:
            self.speedHirtChris = float('nan')

#############################################################################################
def Get_Info(D3fileDB, espfile, Bonesout, Spout, Fls, fileDir, New = False):
    if New == False:
        with open(Bonesout+'.json', 'r') as jsonfile:
            DicBones = json.load(jsonfile)
    BonesList =[]
    with open(D3fileDB, 'r') as csvfile: 
        case = csv.DictReader(csvfile)
        for item in case:
            ID = item["ID"]
            if New == False:
                if ID in DicBones:
                    NCBI = DicBones[ID]["NCBI"]
                    morph_type = DicBones[ID]["morph_type"]
                    Morphotype = DicBones[ID]["Morphotype"]
                    long = DicBones[ID]["Long"]
                    side = DicBones[ID]["side"]
                    BonesList.append(bone(ID,NCBI, side, long, morph_type, Morphotype))
                    BonesList[-1].taxa()
                    BonesList[-1].Rmax = DicBones[ID]["Rmax"]
                    BonesList[-1].Rmin = DicBones[ID]["Rmin"]
                    BonesList[-1].Rmed = DicBones[ID]["Rmed"]
                    BonesList[-1].La = DicBones[ID]["La"]
                    BonesList[-1].RMSD = DicBones[ID]["RMSD"]
                    BonesList[-1].xfitl = DicBones[ID]["xfitl"]
                    BonesList[-1].yfitl = DicBones[ID]["yfitl"]
                    # BonesList[-1].pts = DicBones[ID]["pts"]
                    # BonesList[-1].slopes = DicBones[ID]["slopes"]
                    # BonesList[-1].mass()
            else:
                NCBI = int(item["NCBI"])
                side = item["side"]
                morph_type = int(item["# trochlea"])
                Morphotype = item["Morph_type"]
                if float(item["Long(mm)"])>10e-10:
                    long = float(item["Long(mm)"])*float(item["Mult"])
                else:
                    long = float('nan')
                scale = float(item["Mult"])
                BonesList.append(bone(ID,NCBI, side, long, morph_type, Morphotype))
                BonesList[-1].taxa()
                BonesList[-1].DataJoint(fileDir, scale)
                # BonesList[-1].mass()

            # Save the info for R (o otra cosa)
            HuesosIDs = {}
            for bn in BonesList:
                HuesosIDs.update({str(bn.ID) : bn.__dict__})
            GetDF(HuesosIDs, Bonesout)

    SpeciesList =[]
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
                refspeed = ''
                spdMissing = SpeedHirt(massAvg,*poptHirt)
            else:
                Maxspeed = float(item["Speed max km/hr"])
                refspeed = item["refSpeed"]
                spdMissing = float(item["Speed max km/hr"])
            Stance = item["stance"]
            Diet = item["Diet"]

            SpeciesList.append(Specie(NCBI, massI, massS, massAvg, refmass, Maxspeed, refspeed, spdMissing, Stance, Diet))
            SpeciesList[-1].bones(BonesList, Fls)
            SpeciesList[-1].taxa()
    
    # update bones speed
    for bn in BonesList:
        for sp in SpeciesList:
            if bn.NCBI == sp.NCBI: 
                bn.Maxspeed = sp.Maxspeed

    SpeciesIDs = {}
    for sp in SpeciesList:
        SpeciesIDs.update({str(sp.NCBI) : sp.__dict__})
    GetDF(SpeciesIDs, Spout)

###############################################################################
# Measure bones
def MeasureBones(ID, file_dir, scale, morph_type):
    fileDB = file_dir+'Distal/'+ID+'.vtk'
    print(ID)
    SF, Exists = ReadFile(ID,fileDB,scale, morph_type)
    if Exists:
        return(SF.Rmax, SF.Rmin, SF.Rmed, SF.L, SF.RMSD, SF.xfitl, SF.yfitl)
    else: 
        return(float("nan"), float("nan"), float("nan"), float("nan"), float("nan"), float("nan"), float("nan"))

###############################################################################
def GetBones(NCBI, ListObj, Fls):
    lbn = []
    Ll = []
    Lal = []
    Rminl = []
    Rmaxl = []
    Rmedl = []
    xfitl = []
    yfitl = []
    pts = []
    slopes = []
    Mtpy = []

    with open(Fls,'r') as outfile:
        Filesdic = json.load(outfile)

    for bn in range(len(ListObj)):
        if ListObj[bn].NCBI == NCBI:
            lbn.append(ListObj[bn].ID)
            Ll.append(ListObj[bn].Long)
            Rmaxl.append(ListObj[bn].Rmax)
            Rminl.append(ListObj[bn].Rmin)
            Rmedl.append(ListObj[bn].Rmed)
            Lal.append(ListObj[bn].La)
            xfitl.append(ListObj[bn].xfitl)
            yfitl.append(ListObj[bn].yfitl)
            Mtpy.append(ListObj[bn].morph_type)
            if NCBI in Filesdic:
                Filesdic[NCBI].append(ListObj[bn].ID)
            else:
                Filesdic.update({NCBI: [ListObj[bn].ID]})
    size = 0
    L= 0
    La= 0
    Rmin= 0
    Rmax= 0
    Rmed = 0

    # for the longest
    LrL = 0
    LrLa= 0
    LrRmin= 0
    LrRmax= 0
    LrRmed = 0
    # Lrpts = 0
    # Lrslopes = 0
    Morphotype = 0


    if lbn:
        size =  len(lbn)
        L= sum(Ll)/size
        La= sum(Lal)/size
        Rmin= sum(Rminl)/size
        Rmax= sum(Rmaxl)/size
        Rmed = sum(Rmedl)/size
    
        # for the longest
        LrL = max(Ll)
        LrLa= Lal[Ll.index(max(Ll))]
        LrRmin= Rminl[Ll.index(max(Ll))]
        LrRmax= Rmaxl[Ll.index(max(Ll))]
        LrRmed = Rmedl[Ll.index(max(Ll))]
        # Lrpts = pts[Ll.index(max(Ll))]
        # Lrslopes = slopes[Ll.index(max(Ll))]
        Morphotype = Mtpy[Ll.index(max(Ll))]

    with open(Fls,'w') as outfile:
        json.dump(Filesdic, outfile,indent=4)

    return(lbn, size, L, Rmax, Rmin, Rmed, La, LrL, LrLa, LrRmin, LrRmax, LrRmed, Morphotype)

###############################################################################
# get df in pandas 
def GetDF(infodict, outfile):
    dfinf = pd.DataFrame.from_dict(infodict, orient='index') 
    dfinf.to_csv(outfile+'.csv')
    dfinf.to_pickle(outfile+'.pkl')
    with open(outfile+'.json', 'w') as outfile:
        json.dump(infodict, outfile,indent=4)

# Get order
def get_order(lineage_names):
    ord = 'nota'
    for lin in lineage_names:
        if lin in Orders:
            ord = lin
    return(ord)

# Get Family
def get_family(lineage_names):
    fam = 'nota'
    for lin in lineage_names:
        if lin in Families:
            fam = lin
    return(fam)

# Calculate mass accordin to Christiansen
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