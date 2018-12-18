ncks -A -v lsmask WFD_EI_global_landmask_2d_1x1.nc CARDAMOM_2001_2010_RT_WOO.nc
ncatted -a _FillValue,Median,o,f,-9999 CARDAMOM_2001_2010_RT_WOO.nc
ncap2 -s 'where(lsmask != 1) Median=Median@_FillValue' CARDAMOM_2001_2010_RT_WOO.nc CARDAMOM_2001_2010_RT_WOO_MASKED.nc
ncap2 -s 'defdim("time",1);time[time]=58075.0;time@long_name="Time";' -O CARDAMOM_2001_2010_RT_WOO_MASKED.nc CARDAMOM_2001_2010_RT_WOO_MASKED_TIME.nc
ncatted -a units,time,a,c,"days since 1850-01-01 00:00:00" CARDAMOM_2001_2010_RT_WOO_MASKED_TIME.nc
ncatted -a long_name,time,a,c,time CARDAMOM_2001_2010_RT_WOO_MASKED_TIME.nc
ncatted -a calendar,time,a,c,noleap CARDAMOM_2001_2010_RT_WOO_MASKED_TIME.nc
ncrename -v Median,rt_woo CARDAMOM_2001_2010_RT_WOO_MASKED_TIME.nc
ncap2 -s 'rt_woo[$time, $lat, $lon]=rt_woo' CARDAMOM_2001_2010_RT_WOO_MASKED_TIME.nc CARDAMOM_2001_2010_RT_WOO_MASKED_TIME_3DIM.nc
ncks --mk_rec_dmn time CARDAMOM_2001_2010_RT_WOO_MASKED_TIME_3DIM.nc CARDAMOM_2001_2010_RT_WOO_MASKED_TIME_3DIM_REC.nc
#no tbnds step?
ncatted -a _FillValue,rt_woo,a,c,-9999.0 CARDAMOM_2001_2010_RT_WOO_MASKED_TIME_3DIM_REC.nc
ncatted -a _FillValue,rt_woo,o,f,-9999.0 CARDAMOM_2001_2010_RT_WOO_MASKED_TIME_3DIM_REC.nc
ncatted -a bounds,time,a,c,time_bounds CARDAMOM_2001_2010_RT_WOO_MASKED_TIME_3DIM_REC.nc
ncap2 -O -s 'defdim("tbnds",2) ; time_bnds[time,tbnds]={54788, 58075}' CARDAMOM_2001_2010_RT_WOO_MASKED_TIME_3DIM_REC.nc CARDAMOM_2001_2010_RT_WOO_MASKED_TIME_3DIM_REC_TBNDS.nc
