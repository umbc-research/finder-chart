import pandas as pd
from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import SkyCoord

simbad = Simbad()
simbad.add_votable_fields('ids')
simbad.add_votable_fields('otype')
simbad.add_votable_fields('flux(V)')

#varStar = "HIP 1213"

# Returns the Variable Star's Coordinates and Catalog Number
def getStarData(star=None):

     

    # Ensures the star exists in SIMBAD, returns None, None if star is not found
     if star == None:
            print(f"Source '{star}' not found in SIMBAD.")
            return None, None, None, None, None # Return placeholders if the first source is not found
     else:
          resultTableStar=simbad.query_object(star)

          # Finds the RA and Dec of the Variable Star in h:m:s and deg:':'', respectively
          ra=resultTableStar["RA"][0]
          dec=resultTableStar["DEC"][0]

          ids = resultTableStar['IDS'][0]#.decode('utf-8')
          catalog_ids = [id.strip() for id in ids.split('|')]
          catalogNum = str(selectCatalog(catalog_ids))
          if catalogNum == None:
               catalogNum = resultTableStar["MAIN_ID"]
          
          catalogNumConv = catalogNum.replace(" ", "_")

          # Converts RA and Dec into a SkyCoord object
          #skyCoordStar = SkyCoord(ra, dec, frame="icrs", unit=(u.hourangle, u.deg))

          magnitudeV = resultTableStar["FLUX_V"][0]
          
          return ra, dec, catalogNum, catalogNumConv, magnitudeV # Returns RA and Dec in decimal degrees

# Uses simbad.query_region() to find stars that could be used as references within the same frame (+-1 magnitude ideally)
def getPotentialReferences(coordinates):
          
     potentialRef=simbad.query_region(coordinates, radius="0d6m")
     
     star_results = potentialRef[potentialRef['OTYPE'] == 'Star']

     return star_results

# Selects a catalog number from a list of catalogs an object belongs to
# Hierarchy:
# HD -> HIP -> HIC -> 2MASS -> User-Given ID

def selectCatalog(catalogs):
    
    # Loops through list of catalogs and searches for each one in the hierarchy
    for entry in catalogs:
         if(entry.startswith("HD")):
              return entry
    for entry in catalogs:
         if(entry.startswith("HIP")):
              return entry
    for entry in catalogs:
         if(entry.startswith("HIC")):
              return entry
    for entry in catalogs:
         if(entry.startswith("2MASS")):
              return entry
         
    return

# Separates hh mm and ss or deg ' and '' from one string into three strings
def extractAngle(coordinates):
     coordSeparated = [ang.strip() for ang in coordinates.split(' ')]
     # Removes "+" from Dec degrees (Prevents equations in Google Sheets)
     decSoln = coordSeparated[0].replace("+", "")
     print(coordSeparated)
     return decSoln, coordSeparated[1], coordSeparated[2]


# List of targets that need ID's and Reference Stars
targetList = "C:\\Users\\conno\\Downloads\\VariableStarsTargets_TargetCommons.csv"

inputDF = pd.read_csv(targetList)
#print(commonVarNames)
commonVarNames = []
for i in inputDF.itertuples():
     commonVarNames.append(i[1])

# Variable Star IDs
varIDs = []
# Variable Stars IDs with underscores instead of whitespaces
varIDsConv = []
# Reference Star IDs
refIDs = []
# Reference Star IDs with underscores instead of whitespaces
refIDsConv = []
# Reference RA hours
refRAh = []
# Reference RA minutes
refRAm = []
# Reference RA seconds
refRAs = []
# Reference Dec degrees
refDecd = []
# Reference Dec arcminutes
refDecam = []
# Reference Dec arcseconds
refDecas = []
# Reference Epoch
refEpochs = []
# Reference Star V Magnitudes
refMagsV = []



for target in commonVarNames:

     # Finds the coordinates of a given variable star

     #varRA, varDec, varID, varIDConv, varMagV  = getStarData(varStar if varStar else None)
     varRA, varDec, varID, varIDConv, varMagV  = getStarData(target if target else None)

     if(varRA != None and varDec != None):

          originCoords = SkyCoord(varRA, varDec, frame="icrs", unit=(u.hourangle, u.deg))

          if(originCoords!=None and varID!=None):

     #    print(varRA)
     #    print(varDec)
     #    print(varID)
     #    print(varIDConv)
     #    print(varMagV)

               refComp = 100
               references = getPotentialReferences(originCoords)
               if references is not None:
                    #for row in references:
                         #print(f"MAIN_ID: {row['MAIN_ID']}, FLUX_V: {row['FLUX_V']}")        
                         for row in references:
                              if abs(varMagV-row['FLUX_V']) < refComp:
                                   refComp = abs(varMagV-row['FLUX_V'])
                                   refName=row['MAIN_ID']

               refRA, refDec, refID, refIDConv, refMagV = getStarData(refName if refName else None)
               print(target)
               # Splits RA and Dec into hms and deg arcmin arcsec
               refSplitRA = extractAngle(refRA)
               refSplitDec = extractAngle(refDec)

               # Populates the current row with current target's info
               varIDs.append(varID)
               varIDsConv.append(varIDConv)
               refIDs.append(refID)
               refIDsConv.append(refIDConv)
               refRAh.append(refSplitRA[0])
               refRAm.append(refSplitRA[1])
               refRAs.append(refSplitRA[2])
               refDecd.append(refSplitDec[0])
               refDecam.append(refSplitDec[1])
               refDecas.append(refSplitDec[2])
               refEpochs.append("J2000")
               refMagsV.append(refMagV)
          
          else:
               # Populates the current row with current target's info
               varIDs.append(None)
               varIDsConv.append(None)
               refIDs.append(None)
               refIDsConv.append(None)
               refRAh.append(None)
               refRAm.append(None)
               refRAs.append(None)
               refDecd.append(None)
               refDecam.append(None)
               refDecas.append(None)
               refEpochs.append(None)
               refMagsV.append(None)


df = pd.DataFrame({"Object Name": varIDsConv, "Catalog No.": varIDs ,"Reference Star Object Name": refIDsConv, "Reference Star Catalog No.": refIDs, "hh": refRAh, "mm": refRAm, "ss": refRAs, "deg": refDecd, "'": refDecam, "''": refDecas, "Reference Epoch": refEpochs, "Reference Star Magnitude (V)": refMagsV})
df.to_csv("VariablesExtraInfo.csv", index=False)