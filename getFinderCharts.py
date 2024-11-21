import pandas as pd
from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np
import math

simbad=Simbad()
simbad.add_votable_fields("flux(V)")

def getCenterCoords(objTable=None):
    objRA, objDec = objTable["RA"], objTable["DEC"]
    center = SkyCoord(objRA, objDec, frame="icrs", unit=(u.hourangle, u.deg))
    return center

def extractAngle(coordinates):
     coordSeparated = [ang.strip() for ang in coordinates.split(' ')]

     # Removes "+" from Dec degrees (Prevents equations in Google Sheets)
     decSoln = coordSeparated[0].replace("+", "")

     if len(coordSeparated) == 3:
          return decSoln, coordSeparated[1], coordSeparated[2]
     else:
          return None, None, None

star = "U Aql"
angFoVstr = "0d4m"
angFoV = 4 # arcminutes

starTable = simbad.query_object(star)

centerCoords = getCenterCoords(starTable)
#print(centerCoords.ra.deg)
surroundings = simbad.query_region(centerCoords, angFoVstr)

surroundingsTable = []
surroundingsNames = []
surroundingsRA = []
surroundingsDec = []
surroundingsMag = []

for i in surroundings:
    #print(i["MAIN_ID"])
        starLoc=getCenterCoords(i)
        surroundingsNames.append(i["MAIN_ID"])
        surroundingsRA.append(starLoc.ra.deg)
        surroundingsDec.append(starLoc.dec.deg)

        if isinstance(i["FLUX_V"], np.float32):
            adjMag = 10/math.log(i["FLUX_V"])     
        else:
            adjMag = 1
        
        surroundingsMag.append(pow(adjMag,2))

raArray = np.array(surroundingsRA)

decArray = np.array(surroundingsDec)
sizeArray = np.array(surroundingsMag)

plt.scatter(raArray, decArray, s=sizeArray, color="black")



plt.xlim(centerCoords.ra.deg+(angFoV/60),centerCoords.ra.deg-(angFoV/60))
plt.ylim(centerCoords.dec.deg-(angFoV/60), centerCoords.dec.deg+(angFoV/60))

plt.show()