#!/bin/csh
#

setenv HEADDIR `pwd`
cd src/WW3/model/bin
cp -u ~/Custom_Files/wwatch3.env-workstation wwatch3.env
./w3_setup -q -c Intel -s NCEP_st4 -t `pwd`/../tmp ../.
cp -u ~/Custom_Files/comp-workstation comp
cp -u ~/Custom_Files/link-workstation link
cp -u ~/Custom_Files/switch-workstation switch
foreach sfile (ww3_grid ww3_multi ww3_prnc ww3_ounf)
    if ( ! -d ../$sfile ) then
	w3_source $sfile
	mkdir -p ../$sfile
	cd ../$sfile
	cp -u ../work/{$sfile}.tar.gz .
	gunzip -f {$sfile}.tar.gz
	tar -xf {$sfile}.tar
	mv constants.F90 constants_ww3.F90
	cd ../bin
    endif
end
cd $HEADDIR

# For intel
mkdir -p build/intel
cat <<EOFA > build/intel/env
module load netcdf/4.2
module swap intel_compilers/11.1.073 intel_compilers/18.0.3
EOFA

echo ''
echo ''
echo 'Building FMS'
echo ''
echo ''
#FMS
mkdir -p build/intel/FMSlib/repro/
(cd build/intel/FMSlib/repro/; rm -f path_names;\
../../../../src/mkmf/bin/list_paths -l ../../../../src/FMS; \
../../../../src/mkmf/bin/mkmf -t /home/bgr/Custom_Files/linux-intel-OMPI.mk -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DSPMD " path_names)
(cd build/intel/FMSlib/repro/; make NETCDF=4 REPRO=1 libfms.a -j)

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
../../../../src/mkmf/bin/mkmf -t /home/bgr/Custom_Files/linux-intel-OMPI-WW3.mk -o -p ww3_grid path_names)
(cd build/intel/wave_ice_ocean/ww3_grid/; source ../../env; make NETCDF=3 REPRO=1 ww3_grid -j)

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
../../../../src/mkmf/bin/mkmf -t /home/bgr/Custom_Files/linux-intel-OMPI-WW3.mk -o -p ww3_prnc path_names)
(cd build/intel/wave_ice_ocean/ww3_prnc/; source ../../env; make NETCDF=3 REPRO=1 ww3_prnc -j)

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
../../../../src/mkmf/bin/mkmf -t /home/bgr/Custom_Files/linux-intel-OMPI-WW3.mk -o -p ww3_multi path_names)
(cd build/intel/wave_ice_ocean/ww3_multi/; source ../../env; make NETCDF=3 REPRO=1 ww3_multi -j)
echo ''
echo ''
echo 'Building libww3'
echo ''
echo ''
mkdir -p build/intel/WW3lib/repro
(cd build/intel/WW3lib/repro/; rm -f path_names; \
../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/WW3/model/ww3_multi)
(cd build/intel/WW3lib/repro/; \
../../../../src/mkmf/bin/mkmf -t /home/bgr/Custom_Files/linux-intel-OMPI-WW3.mk -p libww3.a -c "-Duse_libMPI -Duse_netCDF -DSPMD " path_names)
(cd build/intel/WW3lib/repro/; make NETCDF=4 REPRO=1 libww3.a -j)
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
../../../../src/mkmf/bin/mkmf -t /home/bgr/Custom_Files/linux-intel-OMPI.mk -o '-I../../FMSlib/repro -I../../WW3lib/repro' -p MOM6 -l '-L../../FMSlib/repro -lfms -L../../WW3lib/repro -lww3' -c '-Duse_libMPI -Duse_netCDF -DSPMD -Duse_AM3_physics -D_USE_LEGACY_LAND_ ' path_names )
(cd build/intel/wave_ice_ocean/repro/; source ../../env; make NETCDF=3 REPRO=1 MOM6 -j)
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
../../../../src/mkmf/bin/mkmf -t /home/bgr/Custom_Files/linux-intel-OMPI-WW3.mk -o -p ww3_ounf path_names)
(cd build/intel/wave_ice_ocean/ww3_ounf/; source ../../env; make NETCDF=3 REPRO=1 ww3_ounf -j)
