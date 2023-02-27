#Read the CSV file and create the objects of the species with the info
from GetInfoII import *

####################################################################
# File directory
fileDir = '/home/kale/Documents/Allometry/DATA/'
# Species file
espfile = fileDir+'New tax.csv'
Spout = fileDir+'Esp_out'

# Bones file
D3fileDB = fileDir+'3DfilesDb.csv'
Bonesout = fileDir+'Bone_out'

Fls = fileDir+'fls.json'

# call the get info function
Get_Info(D3fileDB, espfile, Bonesout, Spout, Fls, fileDir, New=False)