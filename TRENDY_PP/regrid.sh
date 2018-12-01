#!/bin/bash

path2data="../MODEL_DATA/JULES-ES.1p0.vn5.1.500.CRUJRA1.365.S3.Monthly.2D/"

# shift longitudes such that they run from -179.0625 to 179.0625
#cdo sellonlatbox,-180,180,90,-90 JULES-ES.1p0.vn5.1.500.CRUJRA1.365.S3.Monthly.2D.ILAMB.nc JULES-ES.1p0.vn5.1.500.CRUJRA1.365.S3.Monthly.2D.ILAMB_updated.nc

# reverse latitudes so positive latitudes are first
#cdo invertlat JULES-ES.1p0.vn5.1.500.CRUJRA1.365.S3.Monthly.2D.ILAMB_updated.nc JULES-ES.1p0.vn5.1.500.CRUJRA1.365.S3.Monthly.2D.ILAMB_invlat.nc

# regrid from 1.25x1.875 to 1x1
cdo -O remapcon,mygrid_0.5x0.5 $path2data/JULES-ES.1p0.vn5.1.500.CRUJRA1.365.S3.Monthly.2D.ILAMB_invlat.nc $path2data/JULES-ES.1p0.vn5.1.500.CRUJRA1.365.S3.Monthly.2D.ILAMB_05x05.nc 









































