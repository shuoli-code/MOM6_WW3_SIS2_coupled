import numpy as np
import netCDF4 as NC


grid_handle = NC.Dataset ('ocean_hgrid.nc',mode = 'r' )
topog_handle = NC.Dataset ('ocean_topog.nc',mode = 'r' )

# 2. Read in Super Grid
# 3. Read in Topog

Lon = grid_handle.variables['x'][:,:]
Lat = grid_handle.variables['y'][:,:]
H_DPT = topog_handle.variables['depth'][:,:]

print(np.shape(H_DPT))

H_MSK6=np.ones(np.shape(H_DPT))

H_LAT=np.zeros(np.shape(H_DPT))
H_LON=np.zeros(np.shape(H_DPT))

H_LAT = Lat[1::2,1::2]
H_LON = Lon[1::2,1::2]

LLAT = len(H_LAT[:,0])
LLON = len(H_LAT[0,:])

GRIDNAME = 'OM4_025'

f6=open(GRIDNAME+'.Mask','w')
f7=open(GRIDNAME+'.Dpt','w')
f9=open(GRIDNAME+'.Obstr','w')
f10=open(GRIDNAME+'.Lat','w')
f11=open(GRIDNAME+'.Lon','w')
for ii in np.arange(0,LLAT):
    for jj in np.arange(0,LLON):
        f6.write(str(H_MSK6[ii,jj])+' ')
        f7.write(str(-abs(H_DPT[ii,jj]))+' ')
        f9.write(str(0.)+' ')
        f10.write(str(H_LAT[ii,jj])+' ')
        f11.write(str(H_LON[ii,jj])+' ')
    f6.write('\n')
    f7.write('\n')
    f9.write('\n')
    f10.write('\n')
    f11.write('\n')

for ii in np.arange(0,LLAT):
    for jj in np.arange(0,LLON):
        f9.write(str(0.)+' ')
    f9.write('\n')

f6.close()
f9.close()
f10.close()
f11.close()
