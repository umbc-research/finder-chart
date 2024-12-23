import pandas as pd
from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
import matplotlib.pyplot as plt
import matplotlib.axes as Axis
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import numpy as np
import math

from sys import exit

from finderChartHelpers import *

simbad=Simbad()
simbad.add_votable_fields("flux(V)")

# List of targets that need ID's and Reference Stars
targetList = "FinderChartList.csv"

inputDF = pd.read_csv(targetList)
#print(commonVarNames)
targList = []

for i in inputDF.itertuples():
    targList.append(i[1:])

for targ in targList:
    # Manually input star and reference here

    star = targ[0]
    reference = targ[1]

    # Manually input FoV here

    angFoVstr = "0d8m"
    angFoV = 8 # arcminutes

    # Set up vertical offset for labelling stars

    textOffset = angFoV/300

    # Gather coordinates for variable and reference stars

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
            adjMag = 4/pow(math.log(i["FLUX_V"]), 0.5)     
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

    fig, ax = plt.subplots(figsize=(8,8))
    
    plt.scatter(raArray, decArray, s=sizeArray, color="black")

    # Labels both the variable and the reference star

    plt.annotate("Variable", [centerCoords.ra.deg, centerCoords.dec.deg], [centerCoords.ra.deg, centerCoords.dec.deg+textOffset], horizontalalignment='center', arrowprops=dict(arrowstyle='-'))
    plt.annotate("Reference", [referenceCoords.ra.deg, referenceCoords.dec.deg], [referenceCoords.ra.deg, referenceCoords.dec.deg+textOffset], horizontalalignment='center', arrowprops=dict(arrowstyle='-'))

    # Sets the figure Title to the name of the Variable Star 

    plt.title(star)
    
    # Adds a scalebar to the figure

    scalebar = AnchoredSizeBar(ax.transData, 1/60, "1'", loc='lower left', pad=1, frameon=False)

    ax.add_artist(scalebar)

    plt.legend([], [], loc='lower right', frameon=False, fontsize=12, title=f'FoV: {angFoV*2}\'x{angFoV*2}\'')

    # Reformats the values of axis ticks

    halfFoVSkyCoordX=SkyCoord(ra=angFoV/60, dec=0, frame="icrs", unit=u.deg)
    startTickX = SkyCoord(ra=centerCoords.ra.deg-halfFoVSkyCoordX.ra.deg, dec=0, frame="icrs", unit=u.deg)
    
    # Reformats the tick labels of the x-axis (RA)
    xLocs = []
    xLocsStr = []
    for tick in getLocations(angFoV, startTickX, "RA"):
        strTick = deci2RA(SkyCoord(ra=tick, dec=0, frame="icrs", unit=u.deg).ra)
        xLocs.append(tick[0])
        xLocsStr.append(strTick)

    xlocations, xlabels = plt.xticks()

    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels = ax.get_xticks()
    
    ax.set_xticks(xLocs, xLocsStr, rotation=60, font=dict(size=10))

    # Reformats the tick labels of the y-axis (Dec)

    halfFoVSkyCoordY=SkyCoord(ra=0, dec=angFoV/60, frame="icrs", unit=u.deg)
    startTickY = SkyCoord(ra=0, dec=centerCoords.dec.deg-halfFoVSkyCoordY.dec.deg, frame="icrs", unit=u.deg)
    
    yLocs = []
    yLocsStr = []
    for tick in getLocations(angFoV, startTickY, "DEC"):
        strTick = deci2Dec(SkyCoord(ra=0, dec=tick, frame="icrs", unit=u.deg).dec)
        yLocs.append(tick[0])
        yLocsStr.append(strTick)

    ylocations, ylabels = plt.yticks()

    labels = [item.get_text() for item in ax.get_yticklabels()]
    labels = ax.get_yticks()
    
    ax.set_yticks(yLocs, yLocsStr, font=dict(size=10))
    
    # Sets the axis labels and ranges

    plt.xlabel("RA (hms)")
    plt.xlim(centerCoords.ra.deg+(angFoV/60),centerCoords.ra.deg-(angFoV/60))

    plt.ylabel("Dec (dms)")
    plt.ylim(centerCoords.dec.deg-(angFoV/60), centerCoords.dec.deg+(angFoV/60))
    
    # Set the aspect ratio to 1:1
    plt.gca().set_aspect(1)

    plt.grid(1, alpha=0.7)
    plt.tight_layout()
    
    # Either display or save image here (Comment out to exclude)
    
    #plt.show() # USE THIS TO SEE EACH PLOT

    path = "FinderChartPNGs\\finder_"
    fileName = path+star.replace(" ", "_")

    plt.savefig(fileName) # USE THIS TO SEE EACH PLOT
    plt.close()
