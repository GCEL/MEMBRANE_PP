import netCDF4
import numpy as np

#Open the original INPE landmask as a readable dataset object

mask_dir = "/exports/csce/datastore/geos/groups/gcel/MEMBRANE_database/ANCILLARY_NC_FILES/"

old_filename = mask_dir + "mask-CERRADO.nc"
new_filename = mask_dir + "ilamb-mask-CERRADO.nc"

print(old_filename)
print(new_filename)

#orig_landmask = netCDF4.Dataset(filename, "r")

with netCDF4.Dataset(old_filename) as src, netCDF4.Dataset(new_filename, "w") as dst:
    # Copy across the original variables and dimensions
    """
    for name, dimension in src.dimensions.iteritems():
        dst.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    for name, variable in src.variables.iteritems():
        #processing
        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst.variables[x][:] = src.variables[x][i]
    """

    # Or should we do it manually
    # Create netCDF dimensions
    #dset = Dataset("regions.nc",mode="w")

    # We can copy over exisitng values for lat lon from the existing 

    # Create the lat/lon values using numpy
    # We later write them into the nc file
    lon = src.variables["longitude"][:] 
    lat = src.variables["latitude"][:]

    #res    = 0.5
    #latbnd = np.asarray([np.arange(- 90    , 90     ,res),
    #                     np.arange(- 90+res, 90+0.01,res)]).T
    #lonbnd = np.asarray([np.arange(-180    ,180     ,res),
    #                     np.arange(-180+res,180+0.01,res)]).T
    #lat    = latbnd.mean(axis=1)
    #lon    = lonbnd.mean(axis=1)

    # Create the number array and initialise to missing value
    miss = -9999
    # Take the mask from the original data (Check the name first)
    ids = src.variables["Band1"][:]  # warning: check incoming datatype!
    print("The datatype of incoming landmask is: ", src.variables["Band1"].dtype)

    src_missing_val = src.variables["Band1"].getncattr("missing_value")

    print("The incoming missing_value is: ", src_missing_val)
    print(np.min(ids), np.max(ids))

    # Make sure the mask values match the ILAMB numbering system
    # i.e. 0, 1, 2, ...
    #ids[ids==1] = 0  # set the ones to zero    # doesn't work because comparing floats/ints

    # This is a bit hacky, but deals with the above issue
    ids[np.where(np.logical_and(ids < 2, ids > -1))] = 0

    ids = np.ma.masked_values(ids,src_missing_val)
    ids = ids.astype(np.int64)

    # Create the array of labels
    lbl = np.asarray(["Cerrado"])

    # Create the relevant dimensions
    dst.createDimension("lat" ,size=lat.size)
    dst.createDimension("lon" ,size=lon.size)
    dst.createDimension("nb"  ,size=2       )
    dst.createDimension("n"   ,size=lbl.size)

    # Create netCDF variables
    X  = dst.createVariable("lat"        ,lat.dtype,("lat"      ))
    #XB = dst.createVariable("lat_bounds" ,lat.dtype,("lat","nb" ))
    Y  = dst.createVariable("lon"        ,lon.dtype,("lon"      ))
    #YB = dst.createVariable("lon_bounds" ,lon.dtype,("lon","nb" ))
    I  = dst.createVariable("ids"        ,ids.dtype,("lat","lon"))
    L  = dst.createVariable("labels"     ,lbl.dtype,("n"        ))

    # Load data and encode attributes
    X [...] = lat
    X.units = "degrees_north"
    #XB[...] = latbnd

    Y [...] = lon
    Y.units = "degrees_east"
    #YB[...] = lonbnd

    I[...]  = ids
    I.labels= "labels"

    L[...]  = lbl

    #dst.close()
    







