#!/bin/bash

# add yearly into 1 file
cdo mergetime INLAND_aline/S3/inland-yearly-{2001..2010}.nc inland-yearly-2001_2010.nc

# extract rootbio, woodbio, leafbio, npptot
cdo select,name=rootbio,woodbio,leafbio,npptot,totfall inland-yearly-2001_2010.nc inland-yearly-2001_2010_tmp.nc

# calculate cveg
cdo expr,'cveg=rootbio+woodbio+leafbio;npptot=npptot;totfall=totfall' inland-yearly-2001_2010_tmp.nc inland-yearly-2001_2010_tmp1.nc

# extract tile 3 data
ncks -v cveg,npptot,totfall -d tile,2 inland-yearly-2001_2010_tmp1.nc inland-yearly-2001_2010_tmp2.nc

# delete tile dimension
ncwa -a tile inland-yearly-2001_2010_tmp2.nc inland-yearly-2001_2010_tmp3.nc

# delete tile variable
ncks -C -O -x -v tile  inland-yearly-2001_2010_tmp3.nc inland-yearly-2001_2010_tmp4.nc

# sum over all years
ncwa -a time -y sum  inland-yearly-2001_2010_tmp4.nc inland-yearly-2001_2010_tmp5.nc

# delete time variable
ncks -C -O -x -v time inland-yearly-2001_2010_tmp5.nc inland-yearly-2001_2010_tmp6.nc

# calculate residence times
cdo expr,'cveg=cveg;npptot=npptot;totfall=totfall;rt_npp=cveg/npptot;rt_lit=cveg/totfall' inland-yearly-2001_2010_tmp6.nc inland-yearly-2001_2010_cveg.nc

# add units
ncatted -O -a units,cveg,o,c,"kg.m-2" inland-yearly-2001_2010_cveg.nc
ncatted -O -a units,totfall,o,c,"kg.m-2" inland-yearly-2001_2010_cveg.nc
ncatted -O -a units,npptot,o,c,"kg.m-2.yr-1" inland-yearly-2001_2010_cveg.nc
ncatted -O -a units,rt_npp,a,c,"yrs" inland-yearly-2001_2010_cveg.nc
ncatted -O -a units,rt_lit,a,c,"yrs" inland-yearly-2001_2010_cveg.nc

# regrid from 0.5x0.5 to 1x1
cdo -O remapcon,mygrid_1x1_inland inland-yearly-2001_2010_cveg.nc inland-yearly-2001_2010_cveg_1x1.nc

# tidy up
rm inland-yearly-2001_2010.nc inland-yearly-2001_2010_tmp.nc inland-yearly-2001_2010_tmp1.nc inland-yearly-2001_2010_tmp2.nc inland-yearly-2001_2010_tmp3.nc inland-yearly-2001_2010_tmp4.nc inland-yearly-2001_2010_tmp5.nc inland-yearly-2001_2010_tmp6.nc









