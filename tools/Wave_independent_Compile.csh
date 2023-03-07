#!/bin/csh
#

setenv site `perl -T -e "use Net::Domain(hostdomain); print hostdomain" | sed 's/\.$//'`

#Set 1 to build FMS libraries (needed for coupled model)
setenv BuildFMS 0
#Set 1 to build ww3_grid exec (needed for creating wave model grid)
setenv BuildWW3grid 0
#Set 1 to build ww3_prnc exec (optional, needed to run wave model from external forcing)
setenv BuildWW3prnc 1
#Set 1 to build ww3_multi exec (optional, needed to run wave model in stand-alone from multi driver)
setenv BuildWW3multi 1
#Set 1 to build ww3_shel exec (optional, needed to run wave model in stand-alone from shel driver)
setenv BuildWW3shel 1
#Set 1 to build WW3 library for coupled model
setenv BuildWW3lib 0
#Set 1 to build ww3_ounf exec (optional, needed to process WW3 output to NetCDF)
setenv BuildWW3ounf 0
#Set 1 to build the actual coupled model
setenv BuildMOM6 0

setenv HEADDIR `pwd`
mkdir -p build/intel

# These options correctly set the compiler templates and modules for Gaea or the GFDL workstations.
# Use these as a guide to set-up your own environment, and feel free to add so others may use!
if ( "$site" == "gfdl.noaa.gov" ) then
    setenv TEMPLATE '/home/bgr/Custom_Files/linux-intel-OMPI.mk'
    setenv TEMPLATE_WW3 '/home/bgr/Custom_Files/linux-intel-OMPI-WW3.mk'
    cat <<EOF > build/intel/env
    module load netcdf/4.2
    module swap intel_compilers/11.1.073 intel_compilers/18.0.3
EOF
else if ( "$site" == "ncrc.gov" ) then
    # For intel
    setenv TEMPLATE '../../../../src/mkmf/templates/ncrc-intel.mk'
    setenv TEMPLATE '../../../../src/mkmf/templates/ncrc-intel-WW3.mk'
    cat <<EOF > build/intel/env
    module unload PrgEnv-pgi
    module unload PrgEnv-pathscale
    module unload PrgEnv-intel
    module unload PrgEnv-gnu
    module unload PrgEnv-cray

    module load PrgEnv-intel
    module swap intel intel/16.0.3.210
    module unload netcdf
    module load cray-netcdf
    module load cray-hdf5
EOF
#shuo
else
    setenv TEMPLATE "../../../../src/mkmf/templates/linux-intel.mk"
    setenv TEMPLATE_WW3 "../../../../src/mkmf/templates/linux-intel-WW3.mk"
    cat <<EOF > build/intel/env
    module load openmpi
    module load netcdf
    module load intel-compiler
EOF
#
endif
echo '2'
if ($BuildFMS == 1) then
    echo ''
    echo ''
    echo 'Building FMS'
    echo ''
    echo ''
    #FMS
    mkdir -p build/intel/FMSlib/repro/
    (cd build/intel/FMSlib/repro/; rm -f path_names;\
    ../../../../src/mkmf/bin/list_paths -l ../../../../src/FMS; \
    ../../../../src/mkmf/bin/mkmf -t $TEMPLATE -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DSPMD " path_names)
    (cd build/intel/FMSlib/repro/; source ../../env; make NETCDF=4 REPRO=1 libfms.a -j )
endif

if ($BuildWW3grid == 1) then
    echo ''
    echo ''
    echo 'Building ww3_grid'
    echo ''
    echo ''
    #WW3_grid
    mkdir -p build/intel/wave_ice_ocean/ww3_grid/
    (cd build/intel/wave_ice_ocean/ww3_grid/; rm -f path_names; \
    ../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/WW3/model/ww3_grid)
    (cd build/intel/wave_ice_ocean/ww3_grid/; \
    ../../../../src/mkmf/bin/mkmf -t $TEMPLATE_WW3 -p ww3_grid path_names)
    (cd build/intel/wave_ice_ocean/ww3_grid/; source ../../env; make NETCDF=3 REPRO=1 ww3_grid)
endif

if ($BuildWW3prnc == 1) then
    echo ''
    echo ''
    echo 'Building ww3_prnc'
    echo ''
    echo ''
    #WW3_prnc
    mkdir -p build/intel/wave_ice_ocean/ww3_prnc/
    (cd build/intel/wave_ice_ocean/ww3_prnc/; rm -f path_names;\
    ../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/WW3/model/ww3_prnc)
    (cd build/intel/wave_ice_ocean/ww3_prnc/; \
    ../../../../src/mkmf/bin/mkmf -t $TEMPLATE_WW3 -p ww3_prnc path_names)
    (cd build/intel/wave_ice_ocean/ww3_prnc/; source ../../env; make NETCDF=3 REPRO=1 ww3_prnc)
endif

if ($BuildWW3multi == 1) then
    echo ''
    echo ''
    echo 'Building ww3_multi'
    echo ''
    echo ''
    #Wave
    #WW3_multi
    mkdir -p build/intel/wave_ice_ocean/ww3_multi/
    (cd build/intel/wave_ice_ocean/ww3_multi/; rm -f path_names;\
    ../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/WW3/model/ww3_multi)
    (cd build/intel/wave_ice_ocean/ww3_multi/; \
    ../../../../src/mkmf/bin/mkmf -t $TEMPLATE_WW3 -p ww3_multi path_names)
    (cd build/intel/wave_ice_ocean/ww3_multi/; source ../../env; make NETCDF=3 REPRO=1 ww3_multi)
endif

if ($BuildWW3shel == 1) then
    echo ''
    echo ''
    echo 'Building ww3_multi'
    echo ''
    echo ''
    #Wave
    #WW3_multi
    mkdir -p build/intel/wave_ice_ocean/ww3_shel/
    (cd build/intel/wave_ice_ocean/ww3_shel/; rm -f path_names;\
    ../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/WW3/model/ww3_shel)
    (cd build/intel/wave_ice_ocean/ww3_shel/; \
    ../../../../src/mkmf/bin/mkmf -t $TEMPLATE_WW3 -p ww3_shel path_names)
    (cd build/intel/wave_ice_ocean/ww3_shel/; source ../../env; make NETCDF=3 REPRO=1 ww3_shel)
endif

if ($BuildWW3lib == 1) then
    echo ''
    echo ''
    echo 'Building WW3lib'
    echo ''
    echo ''
    mkdir -p build/intel/WW3lib/repro
    (cd build/intel/WW3lib/repro/; rm -f path_names; \
    ../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/WW3/model/ww3_multi)
    (cd build/intel/WW3lib/repro/; \
    ../../../../src/mkmf/bin/mkmf -t $TEMPLATE_WW3 -p libww3.a -c "-Duse_libMPI -Duse_netCDF -DSPMD " path_names)
    (cd build/intel/WW3lib/repro/; source ../../env; make NETCDF=4 REPRO=1 libww3.a -j)
endif

if ($BuildMOM6 == 1) then
    echo ''
    echo ''
    echo 'Building MOM6'
    echo ''
    echo ''
    #Wave-ice-ocean
    mkdir -p build/intel/wave_ice_ocean/repro/
    (cd build/intel/wave_ice_ocean/repro/; rm -f path_names; \
		../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/MOM6/config_src/{dynamic,coupled_driver} ../../../../src/MOM6/src/{*,*/*}/ ../../../../src/{atmos_null,coupler,land_null,ice_param,icebergs,SIS2,FMS/coupler,FMS/include,WW3/model/CPL}/)
    (cd build/intel/wave_ice_ocean/repro/; \
	../../../../src/mkmf/bin/mkmf -t $TEMPLATE -o '-I../../FMSlib/repro -I../../WW3lib/repro' -p MOM6 -l '-L../../FMSlib/repro -lfms -L../../WW3lib/repro -lww3' -c '-Duse_libMPI -Duse_netCDF -DSPMD -Duse_AM3_physics -D_USE_LEGACY_LAND_ ' path_names )
    (cd build/intel/wave_ice_ocean/repro/; source ../../env; make NETCDF=3 REPRO=1 MOM6 -j)
endif

if ($BuildWW3ounf == 1) then
    echo ''
    echo ''
    echo 'Building ww3_ounf'
    echo ''
    echo ''
    #WW3_ounf
    mkdir -p build/intel/wave_ice_ocean/ww3_ounf/
    (cd build/intel/wave_ice_ocean/ww3_ounf/; rm -f path_names; \
    ../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/WW3/model/ww3_ounf)
    (cd build/intel/wave_ice_ocean/ww3_ounf/; \
    ../../../../src/mkmf/bin/mkmf -t $TEMPLATE_WW3 -p ww3_ounf path_names)
    (cd build/intel/wave_ice_ocean/ww3_ounf/; source ../../env; make NETCDF=3 REPRO=1 ww3_ounf)
endif
