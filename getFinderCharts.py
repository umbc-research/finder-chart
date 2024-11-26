import pandas as pd
from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
import matplotlib.pyplot as plt
import matplotlib.axis as Axis
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

# Input the RA component of a SkyCoord Object
# Output the RA as a string in the hms format
def deci2RA(skyCoordsRA):
    formatCoordsRA = str(skyCoordsRA.hms.h)+"h"+str(skyCoordsRA.hms.m)+"m"+str(round(skyCoordsRA.hms.s, 2))+"s"
    return formatCoordsRA

# Input the Dec component of a SkyCoord Object
# Output the Dec as a string in the dms format
def deci2Dec(skyCoordsDec):
    formatCoordsDec = str(skyCoordsDec.dms.d)+"deg"+str(skyCoordsDec.dms.m)+"'"+str(round(skyCoordsDec.dms.s, 1))+"\""
    return formatCoordsDec

# Manually input star and reference here

star = "HD 236542"
reference = "2MASS_J00500726+6009235"

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
    #print(starLoc.ra.to_string(u.hourangle))
    #print(starLoc.dec.to_string(u.degree))

# Translates lists of star info into arrays to be plotted

raArray = np.array(surroundingsRA)
decArray = np.array(surroundingsDec)
sizeArray = np.array(surroundingsMag)


# Plots the known stars in the field
fig, ax = plt.subplots()
plt.scatter(raArray, decArray, s=sizeArray, color="black")

# Labels both the variable and the reference star

plt.annotate("Variable", [centerCoords.ra.deg, centerCoords.dec.deg], [centerCoords.ra.deg, centerCoords.dec.deg+textOffset], horizontalalignment='center', arrowprops=dict(arrowstyle='-'))
plt.annotate("Reference", [referenceCoords.ra.deg, referenceCoords.dec.deg], [referenceCoords.ra.deg, referenceCoords.dec.deg+textOffset], horizontalalignment='center', arrowprops=dict(arrowstyle='-'))

# Sets the figure Title to the name of the Variable Star 

plt.title(star)

# Sets the axis labels and ranges

plt.xlabel("RA (hms)")
plt.xlim(centerCoords.ra.deg+(angFoV/60),centerCoords.ra.deg-(angFoV/60))

plt.ylabel("Dec (dms)")
plt.ylim(centerCoords.dec.deg-(angFoV/60), centerCoords.dec.deg+(angFoV/60))


# Adds a scalebar to the figure

#ax = plt.gca()
scalebar = AnchoredSizeBar(ax.transData, 1/60, "1'", loc='lower left', pad=1, frameon=False)

ax.add_artist(scalebar)

# Reformats the values of axis ticks

newXLabels = []
xlocations, xlabels = plt.xticks()
#print(xlocations)

labels = [item.get_text() for item in ax.get_xticklabels()]
labels = ax.get_xticks()
for i in labels:
    i=(SkyCoord(ra=i, dec=0, frame="icrs", unit=(u.deg))).ra.to_string(u.hourangle)
    newXLabels.append(i)
ax.set_xticks(xlocations, newXLabels, rotation='vertical', font=dict(size=8))

newYLabels = []
ylocations, ylabels = plt.yticks()

labels = [item.get_text() for item in ax.get_yticklabels()]
labels = ax.get_yticks()
for i in labels:
    i=(SkyCoord(ra=0, dec=i, frame="icrs", unit=(u.hourangle, u.deg))).dec.to_string(u.degree)
    newYLabels.append(i)

ax.set_yticks(ylocations, newYLabels, font=dict(size=8))

# Either display or save image here (Comment out to exclude)

plt.tight_layout()
plt.show()

#path = "FinderChartPNGs\\finder_"
#fileName = path+star.replace(" ", "_")

#plt.savefig(fileName)
