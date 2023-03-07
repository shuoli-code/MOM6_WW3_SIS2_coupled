import netCDF4 as NC

with NC.Dataset("grid_spec.nc") as src, NC.Dataset("grid_spec_waves.nc","w") as dst:
    dst.setncatts(src.__dict__)
    for name, dimension in src.dimensions.items():
        dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else none))
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst[name][:] = src[name][:]
        dst[name].setncatts(src[name].__dict__)

    dst.createDimension('nfile_aXw',1)

#    a=dst.createVariable('wav_mosaic_dir','c',('string'))
#    STRING = './'
#    dst.variables['wav_mosaic_dir'][:len(STRING)] = STRING
#    a.setncattr('standard_name','directory_storying_wave_mosaic')

#    a=dst.createVariable('wav_mosaic_file','c',('string'))
#    STRING = 'wave_mosaic.nc'
#    dst.variables['wav_mosaic_file'][:len(STRING)] = STRING
#    a.setncattr('standard_name','wave_mosaic_file_name')

#    a=dst.createVariable('wav_mosaic','c',('string'))
#    STRING = 'wave_mosaic'
#    dst.variables['wav_mosaic'][:len(STRING)] = STRING
#    a.setncattr('standard_name','wave_mosaic_name')

    a=dst.createVariable('aXw_file','c',('nfile_aXw','string'))
    STRING = 'atmos_mosaic_tile1Xwave_mosaic_tile1.nc'
    dst.variables['aXw_file'][0,:len(STRING)] = STRING
    a.setncattr('standard_name','atmXwav_exchange_grid_file')
