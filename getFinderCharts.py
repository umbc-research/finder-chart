import pandas as pd
from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import numpy as np
import math

simbad=Simbad()
simbad.add_votable_fields("flux(V)")

# Input a Simbad object from simbad.query_object
# Output a SkyCoord object
def getCenterCoords(objTable=None):
    objRA, objDec = objTable["RA"], objTable["DEC"]
    center = SkyCoord(objRA, objDec, frame="icrs", unit=(u.hourangle, u.deg))
    return center


# Input coordinates in string of format (+deg ' '')
# Output each individual coordinate as its own variable
def extractAngle(coordinates):
     coordSeparated = [ang.strip() for ang in coordinates.split(' ')]

     # Removes "+" from Dec degrees (Prevents equations in Google Sheets)
     decSoln = coordSeparated[0].replace("+", "")

     if len(coordSeparated) == 3:
          return decSoln, coordSeparated[1], coordSeparated[2]
     else:
          return None, None, None

# Manually input star and reference here

star = "SY Cas"
reference = "2MASS J00144120+5827154"

# Manually input FoV here

angFoVstr = "0d8m"
angFoV = 8 # arcminutes

# Set up vertical offset for labelling stars
textOffset = angFoV/300

starTable = simbad.query_object(star)
referenceTable = simbad.query_object(reference)

centerCoords = getCenterCoords(starTable)
referenceCoords = getCenterCoords(referenceTable)

# Find stars in a given region

surroundings = simbad.query_region(centerCoords, angFoVstr)

 # Initialize star info lists

surroundingsTable = []
surroundingsNames = []
surroundingsRA = []
surroundingsDec = []
surroundingsMag = []

# Organize info of all stars into lists

for i in surroundings:
    starLoc=getCenterCoords(i)

    # If a star magnitude is known, scale it,
    # Else, make the size of the point equal to 1
        
    if isinstance(i["FLUX_V"], np.float32):
        adjMag = 5/math.log(i["FLUX_V"])     
    else:
        adjMag = 1
    
    # Initializes lists of info
    surroundingsNames.append(i["MAIN_ID"])
    surroundingsRA.append(starLoc.ra.deg)
    surroundingsDec.append(starLoc.dec.deg)
    surroundingsMag.append(pow(adjMag,3))

# Translates lists of star info into arrays to be plotted

raArray = np.array(surroundingsRA)
decArray = np.array(surroundingsDec)
sizeArray = np.array(surroundingsMag)

# Plots the known stars in the field

plt.scatter(raArray, decArray, s=sizeArray, color="black")

# Labels both the variable and the reference star

plt.annotate("Variable", [centerCoords.ra.deg, centerCoords.dec.deg], [centerCoords.ra.deg, centerCoords.dec.deg+textOffset], horizontalalignment='center', arrowprops=dict(arrowstyle='-'))
plt.annotate("Reference", [referenceCoords.ra.deg, referenceCoords.dec.deg], [referenceCoords.ra.deg, referenceCoords.dec.deg+textOffset], horizontalalignment='center', arrowprops=dict(arrowstyle='-'))

# Sets the figure Title to the name of the Variable Star 

plt.title(star)

# Sets the axis labels and ranges

plt.xlabel("RA (deg)")
plt.xlim(centerCoords.ra.deg+(angFoV/60),centerCoords.ra.deg-(angFoV/60))

plt.ylabel("Dec (deg)")
plt.ylim(centerCoords.dec.deg-(angFoV/60), centerCoords.dec.deg+(angFoV/60))

# Adds a scalebar to the figure

ax = plt.gca()
scalebar = AnchoredSizeBar(ax.transData, 1/60, "1 minute", loc='lower left', pad=1, frameon=False)

ax.add_artist(scalebar)

# Either display or save image here (Comment out to exclude)

#plt.show()

path = "FinderChartPNGs\\finder_"
fileName = path+star.replace(" ", "_")

plt.savefig(fileName)
