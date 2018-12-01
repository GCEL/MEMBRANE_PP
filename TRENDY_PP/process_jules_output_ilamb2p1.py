"""
process jules data for use by ILAMB vn2.1
"""
import sys
import jules
from iris.time import PartialDateTime
import subprocess
import iris
import glob
import cf_units
import numpy.ma as ma
import numpy as np
from copy import deepcopy
import os

forceExtract = False # True always use ncrcat to get variables from raw JULES output
                     # False skip this step if it has been done before

ignoreFinalFile = True #True remove final file from list of inputs
                       #False use all file in "filesIn"
                       #For use when JULES has produced an empty file with zero times in
extend_latitudes = True #True add an extra row of points to the south edge of gridded data
                        #False do nothing extra
                        #To prevent ILAMB regriding interpolating data from South America
                        # to the South Pole, or Greenland northward

#DEFAULT OPTIONAL INPUTS
varListDefault = ["cv","gpp_gb","cs_gb","le_flux","runoff","ftl_gb","fqw_gb",
                  "sw_net","sw_down","lw_up","lw_net","rad_net",
                  "t1p5m_gb","precip","lw_down","nep","nee","burntArea","lai",
                  "nbp","reco","EvapFrac","twsa","rhums","albedo"]

def define_derived_vars():
    """
    this is a function at the top of the code to allow new derived variables to be added easily
    """
    #SOME ILAMB VARIABLES CAN NOT BE DIRESTLY OUTPUT FROM JULES
    #dictionary = ilamb var: [[jules vars], conversion function]
    global derivedVar
    derivedVar = {"nep"      : [["npp_n_gb","resp_s_to_atmos_gb"], nep_func],
                  "nee"      : [["npp_n_gb","resp_s_to_atmos_gb"], nee_func],
                  "burntArea": [["burnt_area_gb"], burntArea_func],
                  "lai"      : [["lai","frac"], lai_func],
                  "nbp"      : [["npp_n_gb","WP_fast_out","WP_med_out","WP_slow_out","resp_s_to_atmos_gb","harvest_gb"], nbp_func],
                  "reco"     : [["resp_p_gb","resp_s_to_atmos_gb"], reco_func],
                  "EvapFrac" : [["latent_heat","ftl_gb"], EvapFrac_func],
                  "twsa"     : [["zw","canopy_gb","snow_grnd_gb","smc_tot"], twsa_func],
                  "rhums"    : [["q1p5m_gb","t1p5m_gb","pstar"], rhums_func],
                  "albedo"   : [["sw_down","sw_net"], albedo_func], 
                  "le_flux"  : [["latent_heat"], le_flux_func]
                  }

def get_test_inputs():
    """
    It is expected that this program will be called from a rose suite 
    which will provide the necesary inputs.
    However, if you want to run it on its own, this function can be used
    to provide the inputs. And from_shell should be set to false in the 
    call to "main" at the bottom of this file.
    """
    filesIn = "/hpc/data/d05/eroberts/jules_output/u-at275/at275_ilambTest2.S3.ilamb.*.nc"
    procDir = "/scratch/eroberts/ILAMB_2.1/PROC"
    outDir = "/scratch/eroberts/ILAMB_2.1/MODEL_DATA"
    outName = "at275_ilambTest2.S3"
    monStart  = 1     # =None: use all times
    yearStart = 1980  # =None: use all times
    monEnd    = 12    # =None: use all times
    yearEnd   = 2013  # =None: use all times
    varListIn = ["sw_net"] # =None: use varListDefault
    return filesIn,procDir,outDir,outName,monStart,yearStart,monEnd,yearEnd,varListIn

def main(from_shell=True):
    if from_shell: #read inputs from shell script
##should use this method
 #   ncycle = os.environ['NCYCLE']
#    seasons = os.environ['SEASONS'].split()
##
        filesIn   = os.environ['filesIn']
        #account for filesIn being passed in with speachmarks
        #useful to stop rose expanding the list of file names
        if filesIn[0] in ['"',"'"]: filesIn = filesIn[1:-1]
        procDir   = os.environ['procDir']
        outDir    = os.environ['ILAMB_MODEL_DIR']
        outName   = os.environ['outName']
        monStart  = int(os.getenv('ILAMB_START_MON',0))
        yearStart = int(os.getenv('ILAMB_START_YEAR',0))
        monEnd    = int(os.getenv('ILAMB_END_MON',0))
        yearEnd   = int(os.getenv('ILAMB_END_YEAR',0))
        varListIn,shellError = subprocess.Popen("ilamb_setup.sh",
                                     shell=True, 
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()
        varListIn = varListIn.split()
        print filesIn
        print procDir
        print outDir
        print outName
        print monStart
        print yearStart
        print monEnd
        print yearEnd
        print varListIn
    else: #read inputs from above function
        filesIn,procDir,outDir,outName,monStart,yearStart,monEnd,yearEnd,varListIn = get_test_inputs()
    #Call processing script
    if varListIn:
        print "using user provided var list"
        proc_jules(filesIn,procDir,outDir,outName,monStart,yearStart,monEnd,yearEnd,varList=varListIn)
    else:
        print "using default var list"
        proc_jules(filesIn,procDir,outDir,outName,monStart,yearStart,monEnd,yearEnd)

def proc_jules(filesIn,procDir,outDir,outName,monStart,yearStart,monEnd,yearEnd,varList=varListDefault,makeNewDir=False):
    if makeNewDir:
        print "WARNING: ILAMB may not find the output if makNewDir=True"
    define_derived_vars()
    procDir = procDir+"/"+outName
    if not os.path.exists(procDir):
        os.makedirs(procDir)
    outDir0 = deepcopy(outDir+"/"+outName)
    i = 1
    outDir = deepcopy(outDir0)
    while i > 0: 
        if os.path.exists(outDir) and makeNewDir:
            outDir = outDir0+"_"+str(i)
            i = i + 1
        elif not os.path.exists(outDir):
            os.makedirs(outDir)
            break
        else:
            break
    if ignoreFinalFile: filesIn = ' '.join(sorted(glob.glob(filesIn))[:-1]) #JULES sometimes outputs empty files at end of run
    if any(item == 0 for item in [monStart,yearStart,monEnd,yearEnd]):
        print "process whole timeseries"
        constrainTime = False
        timeName = "."
    else:
        constrainTime = True
        timeName = "."+str(yearStart)+"{:0>2}".format(monStart)+"-"+str(yearEnd)+"{:0>2}".format(monEnd)+"."
        timeConstraint = iris.Constraint(time=lambda cell:                               \
                                         PartialDateTime(year=yearStart, month=monStart) \
                                         <= cell <=                                      \
                                         PartialDateTime(year=yearEnd, month=monEnd))
    print PartialDateTime(year=yearStart, month=monStart), PartialDateTime(year=yearEnd, month=monEnd)
    outNameTime = outName + timeName
    print outNameTime
    firstSubVar = True #is this the first variable of the loop?
    for varIn in varList:
        print varIn
        if varIn in derivedVar.keys():
            subVars = derivedVar[varIn][0]
            derivedFunc = derivedVar[varIn][1]
            derivedCubeList = iris.cube.CubeList([])
        else:
            subVars = [varIn]
            derivedFunc = None
        for var in subVars:
            procVarDir = procDir+"/"+var
            if not os.path.exists(procVarDir):
                os.makedirs(procVarDir)
            varFile = procVarDir+"/"+outName+"."+var+".nc"
            if not glob.glob(varFile) or forceExtract:
                retcode = subprocess.call("ncrcat -O -v "+var+",latitude,longitude,time_bounds "+filesIn+" "+varFile,shell=True)
                if retcode != 0: 
                    print "ncrcat -O -v "+var+",latitude,longitude,time_bounds "+filesIn+" "+varFile
                    exit()
            varConstraint = iris.Constraint(cube_func=lambda x: x.var_name == var)
            cube = jules.load_cube(varFile, varConstraint, missingdata=ma.masked, callback=add_time_middle)
            cube = check_units(cube)
            if constrainTime:
                cube.coord("time").bounds = None  
                with iris.FUTURE.context(cell_datetime_objects=True):
                    cube = cube.extract(timeConstraint) 
                cube.coord('time').guess_bounds()
            if extend_latitudes:
                cube = add_lats(cube)
            if derivedFunc != None:
                derivedCubeList.append(cube)
            else:
                outVarDir = outDir+"/"+var
                if not os.path.exists(outVarDir):
                    os.makedirs(outVarDir)
                varFileOut = outVarDir+"/"+outNameTime+var+".nc"
                iris.save(cube,varFileOut)
                print "SAVED: ",varFileOut
            firstSubVar = False
        if derivedFunc != None:
            derivedCube = derivedFunc(derivedCubeList)
            derivedCube.var_name = varIn
            outVarDir = outDir+"/"+varIn
            if not os.path.exists(outVarDir):
                os.makedirs(outVarDir)
            varFileOut = outVarDir+"/"+outNameTime+varIn+".nc"
            iris.save(derivedCube,varFileOut)
            print "SAVED: ",varFileOut

def albedo_func(cubeList):
    """
    should check that the input cubes are the variables expected.
    albedo = up/down = (down-net)/down
    """
    return (cubeList[0]-cubeList[1])/cubeList[0]

def le_flux_func(cubeList):
    """
    just need to rename variable
    because ILAMB can't differentiate between "lat" and "latent_heat"
    """
    return cubeList[0]

def rhums_func(cubeList):
    """
    calculate relative humidity from specific humidity
    inputs are 1.5m q, 1.5m T and p*
    """
    q = cubeList[0]
    T = cubeList[1]
    p = cubeList[2]
    #saturated vapor pressure(?) = es
    es = 6.1078 * np.exp((17.26938818 * (T.data - 273.15)) / (237.3 + (T.data - 273.15)))
    #calculate relative humidity
    hurs = p.copy(q.data * p.data / (es * (0.622 - (1.0 - 0.622) * q.data)))
    hurs.units = cf_units.Unit("%")
    return hurs

def twsa_func(cubeList):
    """
    variation in total water mass
    (-1 * depth to water table) * water density +
    canopy water + snow + soil moisture
    should check that the input cubes are the variables expected.
    """
    water_density = 1000.0 #kg m-3
    return cubeList[1].copy(-1.0 * cubeList[0].data) * water_density + cubeList[1] + cubeList[2] + cubeList[3]
    
def EvapFrac_func(cubeList):
    """
    should check that the input cubes are the variables expected.
    """
    cubeList[0] = cubeList[0] / (cubeList[0] + cubeList[1])
    cubeList[0].data[cubeList[0].data > 1.0] = np.nan
    cubeList[0].data[cubeList[0].data < 0.0] = np.nan
    return cubeList[0]

def nbp_func(cubeList):
    """
    should check that the input cubes are the variables expected.
    assumes first cube is npp and all others are loss terms
    """
    ncube = len(cubeList)
    print cubeList[0].var_name, cubeList[0].units
    for i in range(1,ncube): 
        print cubeList[i].var_name, cubeList[i].units
        cubeList[0] = cubeList[0] - cubeList[i]
    return cubeList[0]

def nee_func(cubeList):
    """
    should check that the input cubes are the variables expected.
    """
    return -1.0 * nbp_func(cubeList)

def reco_func(cubeList):
    """
    should check that the input cubes are the variables expected.
    """
    return cubeList[0] + cubeList[1]

def nep_func(cubeList):
    """
    should check that the input cubes are the variables expected.
    """
    return nbp_func(cubeList)

def burntArea_func(cubeList):
    """
    converts units from "fraction of land per second"
    to "% of land per year"
    """
    cubeList[0].data = cubeList[0].data*365.0*86400.0*100.0
    cubeList[0].units = cf_units.Unit("%")
    return cubeList[0]
    
def lai_func(cubeList):
    """
    assumes JULES's cubes have "tile" as second dimension:
    time, tile, latitude, longitude
    """
    npft = np.shape(cubeList[0])[1]
    return cubeList[0].collapsed("pft",
                                 iris.analysis.SUM,
                                 weights=cubeList[1].data[:,:npft])

def check_units(cube):
    print cube.var_name
    units = str(cube.units)
    print units
    if "360" in units:
        cube.data = cube.data/(360.0*86400.0)
        units = units.replace("days","")
        if "per" in units:
            units = units.replace("per","")
            units = units.replace("360","s-1")
        else:
            units = units.replace("360","s")
    units = units.replace("C","")
    units = units.replace("N","")
    cube.units = cf_units.Unit(units)
    print units
    return cube

def add_lats(cube):
    """
    Add a row of masked data to the south 
    and a row to the north
    """
    ndims = len(np.shape(cube.data))
    if ndims == 2: #assume dims are [lat,long]
        south = cube[0:1]
        north = cube[-1:]
    elif ndims == 3:
        south = cube[:,0:1]
        north = cube[:,-1:]
    elif ndims == 4:
        south = cube[:,:,0:1]
        north = cube[:,:,-1:]
    south.data.data[:] = np.nan
    north.data.data[:] = np.nan
    south.data.mask[:] = 1
    north.data.mask[:] = 1
    dlat = cube.coord("latitude").points[1] - cube.coord("latitude").points[0]
    south.coord("latitude").points = south.coord("latitude").points[0] - dlat
    north.coord("latitude").points = north.coord("latitude").points[0] + dlat
    if south.coord("latitude").bounds is not None:
        south.coord("latitude").bounds -= dlat
        north.coord("latitude").bounds += dlat
    cube = iris.cube.CubeList([cube,south,north]).concatenate_cube()
    return cube

def make_grid_cube(inCube):
    """
    make an example cube with an extra row of NAN data to the south 
    and an extra row to the north
    """
    cube = deepcopy(inCube[0]) # assumes inCube has 3 time dimension
    south = cube[0:1]          # assumes inCube has no z dimension
    north = cube[-1:]        # assumes inCube has no z dimension
    dlat = cube.coord("latitude").points[1] - cube.coord("latitude").points[0]
    south.coord("latitude").points = south.coord("latitude").points[0] - dlat
    north.coord("latitude").points = north.coord("latitude").points[0] + dlat
    if south.coord("latitude").bounds is not None:
        south.coord("latitude").bounds -= dlat
        north.coord("latitude").bounds += dlat
    cube = iris.cube.CubeList([cube,south,north]).concatenate_cube()
    return cube
                
def add_time_middle(cube, field, filename):
    ''' Add a new coord, with points that are the mid values of the bounds on the coord given as argument'''
    cube.coord("time").points = (cube.coord("time").bounds[:,0] + cube.coord("time").bounds[:,1])/2
        
if __name__ == '__main__':
#    retcode = subprocess.call("ilamb_setup.sh",shell=True)
    main()
