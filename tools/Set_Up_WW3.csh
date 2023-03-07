#!/bin/csh
#

# This needs to reflect your WW3 compiler set.  Look in the WW3/model/bin directory for further details.
# We won't compile the wavemodel with WW3 compile scripts, but we still need to direct the w3_source
#  to valid compilers.
# setenv comp 'cray_xc.Intel'
# setenv switch 'NCEP_st4'
setenv comp 'Intel'
# setenv switch 'UoM_nl1'
# setenv switch 'NCEP_st4'
setenv switch 'UoM_nl1_qx'
#--------------------
#The following generates, unpacks, and sets up folders for wave model executables and libraries needed for coupling to MOM6.

setenv HEADDIR `pwd`
cd src/WW3/model/
setenv WW3DIR `pwd`

cd bin
./w3_setup -q -c $comp -s $switch -t $WW3DIR/tmp $WW3DIR
foreach sfile (ww3_grid ww3_multi ww3_prnc ww3_ounf ww3_shel)
    ./w3_source $sfile
    if ( ! -d ../$sfile ) then
	mkdir -p ../$sfile
    endif
    cd ../$sfile
    cp -u ../work/{$sfile}.tar.gz .
    gunzip -f {$sfile}.tar.gz
    tar -xf {$sfile}.tar
    mv constants.F90 constants_ww3.F90
    cd ../bin
end
cd $HEADDIR
