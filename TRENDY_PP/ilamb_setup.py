import os
def main():
    varConf = {'cv'         :'Biomass',
               'gpp_gb'     :'GrossPrimaryProductivity',
               'cs_gb'      :'SoilCarbon',
               'le_flux'    :'LatentHeat',
               'runoff'     :'Runoff',
               'ftl_gb'     :'SensibleHeat',
               'fqw_gb'     :'Evapotranspiration',
               'sw_net'     :'SurfaceNetSWRadiation',
               'sw_down'    :'SurfaceDownwardSWRadiation',
               'lw_up'      :'SurfaceUpwardLWRadiation',
               'lw_net'     :'SurfaceNetLWRadiation',
               'rad_net'    :'SurfaceNetRadiation',
               't1p5m_gb'   :'SurfaceAirTemperature',
               'precip'     :'Precipitation',
               'lw_down'    :'SurfaceDownwardLWRadiation',
               'nee'        :'NetEcosystemExchange',
               'burntArea'  :'BurnedArea',
               'lai'        :'LeafAreaIndex',
               'nbp'        :'GlobalNetEcosystemCarbonBalance',
               'reco'       :'EcosystemRespiration',
               'EvapFrac'   :'EvaporativeFraction',
               'twsa'       :'TerrestrialWaterStorageAnomaly',
               'rhums'      :'SurfaceRelativeHumidity',
               'albedo'     :'Albedo',
               'pfarea'     :'PermafrostExtent',
               'swe'        :'SnowWaterEquivalent',
               'sw_up'      :'SurfaceUpwardSWRadiation'
               }
    if os.environ['PROCESS_ALL'] == "True":
        confrontations = ''
    else:
        confrontations = ''
        for var,conf in varConf.items():
            if os.environ[var] == "True":
                confrontations = confrontations+" "+conf
    print confrontations

if __name__ == '__main__':
    main()

