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

star = "eta Aql"
angFoVstr = "0d6m"
angFoV = 6 # arcminutes

starTable = simbad.query_object(star)

centerCoords = getCenterCoords(starTable)
#print(centerCoords.ra.deg)
surroundings = simbad.query_region(centerCoords, angFoVstr)

surroundingsTable = []
surroundingsNames = []
surroundingsRA = []
surroundingsDec = []
surroundingsMag = []

#print(surroundings["RA"], surroundings["DEC"])

for i in surroundings:
    #print(i["MAIN_ID"])
        starLoc=getCenterCoords(i)
        surroundingsNames.append(i["MAIN_ID"])
        surroundingsRA.append(starLoc.ra.deg)
        surroundingsDec.append(starLoc.dec.deg)
        print(i["FLUX_V"])
        adjMag = 10/math.log(i["FLUX_V"])
        surroundingsMag.append(pow(adjMag,2))
#print(surroundings)
#print(surroundingsNames)

#print(f"\nMagnitude After Log:\n{surroundingsMag}")

raArray = np.array(surroundingsRA)

decArray = np.array(surroundingsDec)
sizeArray = np.array(surroundingsMag)

plt.scatter(raArray, decArray, s=sizeArray, color="black")

plt.xlim(centerCoords.ra.deg+(angFoV/60),centerCoords.ra.deg-(angFoV/60))
plt.ylim(centerCoords.dec.deg-(angFoV/60), centerCoords.dec.deg+(angFoV/60))

plt.show()