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

star = "SY Cas"
reference = "2MASS J00144120+5827154"
angFoVstr = "0d8m"
angFoV = 8 # arcminutes
textOffset = angFoV/300

starTable = simbad.query_object(star)
referenceTable = simbad.query_object(reference)

centerCoords = getCenterCoords(starTable)
referenceCoords = getCenterCoords(referenceTable)

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

        if isinstance(i["FLUX_V"], np.float32):
            adjMag = 5/math.log(i["FLUX_V"])     
        else:
            adjMag = 1
    
        surroundingsNames.append(i["MAIN_ID"])
        surroundingsRA.append(starLoc.ra.deg)
        surroundingsDec.append(starLoc.dec.deg)
        surroundingsMag.append(pow(adjMag,3))

raArray = np.array(surroundingsRA)

decArray = np.array(surroundingsDec)
sizeArray = np.array(surroundingsMag)

plt.scatter(raArray, decArray, s=sizeArray, color="black")

plt.annotate("Variable", [centerCoords.ra.deg, centerCoords.dec.deg], [centerCoords.ra.deg, centerCoords.dec.deg+textOffset], horizontalalignment='center', arrowprops=dict(arrowstyle='-'))
plt.annotate("Reference", [referenceCoords.ra.deg, referenceCoords.dec.deg], [referenceCoords.ra.deg, referenceCoords.dec.deg+textOffset], horizontalalignment='center', arrowprops=dict(arrowstyle='-'))

plt.title(star)

plt.xlabel("RA (deg)")
plt.xlim(centerCoords.ra.deg+(angFoV/60),centerCoords.ra.deg-(angFoV/60))

plt.ylabel("Dec (deg)")
plt.ylim(centerCoords.dec.deg-(angFoV/60), centerCoords.dec.deg+(angFoV/60))

plt.show()