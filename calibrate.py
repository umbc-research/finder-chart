# Author: Connor stole from Roy
# Before Use: Provide all locations of the files in lines 36 - 46,
#             and change the output location to the desired directory


import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

from scipy.ndimage import gaussian_filter
from scipy.signal import resample

from glob import glob

from photutils.detection import DAOStarFinder

def extractRadialData(subFrame, xC, yC):
    #Get matrix of integer indices associated with subFrame
    y, x = np.indices((subFrame.shape))

    #Generate matrix of radius values
    r = np.sqrt((x - xC)**2 + (y - yC)**2)

    #Force integer values (np.sqrt gives floats)
    r = r.astype(int)

    #Generate a histogram of radius bin, each weighed by corresponding counts
    weightedRadiusHistogram = np.bincount(r.ravel(), weights=subFrame.ravel())
    unweightedRadiusHistogram = np.bincount(r.ravel())

    #Get average for each radius bin
    averageCountsPerRadiusBin = weightedRadiusHistogram / unweightedRadiusHistogram
    return averageCountsPerRadiusBin




fileLocation = "C:\\Users\\conno\\Documents\\UMBCObservatory\\Calibration\\rawData\\Comet"

directories = glob(fileLocation+"\\*")

darks = fileLocation+"\\Darks\\*\\*.fits"

flatdarks = fileLocation+"\\FlatDarks\\*\\*.fits"

flats = fileLocation+"\\Flats\\*\\*.fits"

lights = fileLocation+"\\Light\\*\\*.fits"
print(lights)
# Import all calibration frames and light frames into their own arrays

darkData = [fits.open(d)[0].data for d in glob(darks)]
flatDarkData = [fits.open(d)[0].data for d in glob(flatdarks)]
flatData = [fits.open(d)[0].data for d in glob(flats)]
lightData = [fits.open(d)[0].data for d in glob(lights)]

# Show/Save example of each type
#plt.imshow(lightData[0])
#plt.show()

# Perform Calibration

## Read in all flat darks
#DONE

## Average all flat darks into master_flatdark
master_flatDark = np.mean(flatDarkData, axis=0)

## Read in all flats
#DONE

## subtract master_faltdark from each flat
## Average result into master_flat * Ct
master_flat = np.mean(flatData-master_flatDark, axis=0)

## Identify the constant Ct
Ct = np.mean(master_flat)

## Divide by Ct to get master_flat
master_flat /= Ct

## Read in all darks
#DONE

## average all darks into master_dark
master_dark = np.mean(darkData, axis=0)

## read in all lights
#DONE

## subtract master dark from each light and divide result by master flat
## average results into Science Frame

scienceFrame = np.mean((lightData-master_dark)/master_flat, axis=0)
#scienceFrame = np.sum((lightData-master_dark)/master_flat, axis=0)

plt.imshow(scienceFrame, cmap='gray')
plt.show()

hdu = fits.PrimaryHDU(data=scienceFrame)
hdul = fits.HDUList([hdu])
hdul.writeto(".\\calibratedData\\Comet\\comet_cal.fits")
