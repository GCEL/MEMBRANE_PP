#!/bin/bash

path2data="INLAND_aline/S3/"

#cd $path2data

#cdo mergetime inland-monthly-{2000..2016}.nc inland-monthly-2000-2016.nc

# regrid from 0.5x0.5 to 1x1
cdo -O remapcon,mygrid_1x1_inland INLAND_aline/S3-mon-cor.nc S3-mon-cor_1x1.nc



















