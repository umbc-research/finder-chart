import pandas as pd
from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import SkyCoord

simbad=Simbad()

def getObjTable(Obj=None):
    if Obj == None:
        print("Cannot find provided star")
    else:
        return simbad.query_object(Obj)

star = "HD 143454"

starTable = getObjTable(star)

print(starTable)