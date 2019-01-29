import netCDF4
import numpy as np

#Open the original INPE landmask as a readable dataset object

old_filename = "mask-AM.nc"
new_filename = "ilamb-mask-AM.nc"

#orig_landmask = netCDF4.Dataset(filename, "r")

with netCDF4.Dataset(old_filename) as src, netCDF4.Dataset(new_filename) as dst:
    # Copy across the original variables and dimensions
    for name, dimension in src.dimensions.iteritems():
        dst.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    for name, variable in src.variables.iteritems():
        #processing
        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst.variables[x][:] = src.variables[x][i]

    # Or should we do it manually
    # Create netCDF dimensions
    dset = Dataset("regions.nc",mode="w")
    dset.createDimension("lat" ,size=lat.size)
    dset.createDimension("lon" ,size=lon.size)
    dset.createDimension("nb"  ,size=2       )
    dset.createDimension("n"   ,size=lbl.size)

    # Create netCDF variables
    X  = dset.createVariable("lat"        ,lat.dtype,("lat"      ))
    XB = dset.createVariable("lat_bounds" ,lat.dtype,("lat","nb" ))
    Y  = dset.createVariable("lon"        ,lon.dtype,("lon"      ))
    YB = dset.createVariable("lon_bounds" ,lon.dtype,("lon","nb" ))
    I  = dset.createVariable("ids"        ,ids.dtype,("lat","lon"))
    L  = dset.createVariable("labels"     ,lbl.dtype,("n"        ))

    # Load data and encode attributes
    X [...] = lat
    X.units = "degrees_north"
    XB[...] = latbnd

    Y [...] = lon
    Y.units = "degrees_east"
    YB[...] = lonbnd

    I[...]  = ids
    I.labels= "labels"

    L[...]  = lbl

    dset.close()
    







