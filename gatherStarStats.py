import pandas as pd
from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import SkyCoord

# List of targets that need ID's and Reference Stars
targetList = "C:\\Users\\conno\\Downloads\\VariableStarsTargets_TargetCommons.csv"

df=pd.read_csv(targetList)
print(df)

simbad = Simbad()
simbad.add_votable_fields('ids')
simbad.add_votable_fields('otype')
simbad.add_votable_fields('flux(V)')




varStar = "HIP 1213"

# Returns the Variable Star's Coordinates and Catalog Number
def getStarData(star=None):

     resultTableStar=simbad.query_object(star)

    # Ensures the star exists in SIMBAD, returns None, None if star is not found
     if resultTableStar is None:
            print(f"Source '{star}' not found in SIMBAD.")
            return None, None # Return placeholders if the first source is not found
     else:
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
          skyCoordStar = SkyCoord(ra, dec, frame="icrs", unit=(u.hourangle, u.deg))

          magnitudeV = resultTableStar["FLUX_V"][0]
          
          return skyCoordStar, catalogNum, catalogNumConv, magnitudeV # Returns RA and Dec in decimal degrees

# Uses simbad.query_region() to find stars that could be used as references within the same frame (+-1 magnitude ideally)
def getPotentialReferences(coordinates):
          
     potentialRef=simbad.query_region(coordinates, radius="0d8m")
     
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

# Finds the coordinates of a given variable star
originCoords, varID, varIDConv, varMagV,  = getStarData(varStar if varStar else None)

if(originCoords!=None and varID!=None):

    varRA = originCoords.ra.deg
    varDec = originCoords.dec.deg

    print(varRA)
    print(varDec)
    print(varID)
    print(varIDConv)
    print(varMagV)

    refComp = 100
    references = getPotentialReferences(originCoords)
    if references is not None:
         #for row in references:
              #print(f"MAIN_ID: {row['MAIN_ID']}, FLUX_V: {row['FLUX_V']}")        
          for row in references:
              if abs(varMagV-row['FLUX_V'] < refComp):
                   refName=row['MAIN_ID']

    refCoords, refID, refIDConv, refMagV = getStarData(refName if refName else None)
    refRA = refCoords.ra.deg
    refDec = refCoords.dec.deg

    print(refRA)
    print(refDec)
    print(refID)
    print(refIDConv)
    print(refMagV)


