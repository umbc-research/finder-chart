from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import SkyCoord

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

def getLocations(halfFoVmin, startingTick, coordsSystem):
    halfFoVdeg = halfFoVmin/60
    #print(f"Half the FoV is {halfFoVdeg} degrees, or {halfFoVmin} minutes")
    locations = []
    #startFormat = SkyCoord(ra=0, dec=startingTick, frame="icrs", unit=(u.deg))
    #print(f"Starting Tick = {startingTick.dec.dms}")

    counter = 0
    if(coordsSystem == "RA"):
        while counter <= halfFoVmin:
            locations.append(startingTick.ra.deg + 1/30)
            counter = counter + 1

    if(coordsSystem == "DEC"):
            while counter < halfFoVmin:
                locations.append(startingTick.dec.deg + counter/30)
                counter = counter + 1

    return locations
