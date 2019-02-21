import numpy as np
import netCDF4
import sys

fname = "grid_Arctic_1.nc"
f2name = "popgrid_Arctic.nc"
ncid = netCDF4.Dataset(fname, "r")
nc = netCDF4.Dataset(f2name, "w", format='NETCDF3_CLASSIC')

ulon = ncid.variables['lon_psi'][1:,1:]
ulat = ncid.variables['lat_psi'][1:,1:]
#lon = ncid.variables['lon_rho'][:]
#lat = ncid.variables['lat_rho'][:]
pm = ncid.variables['pm'][:]
pn = ncid.variables['pn'][:]
kmt = ncid.variables['mask_rho'][1:-1,1:-1]
angle = ncid.variables['angle'][:]

ncid.close()

#lon = lon * np.pi/180.0
#lat = lat * np.pi/180.0
ulon = ulon * np.pi/180.0
ulat = ulat * np.pi/180.0
angleu = 0.25*(angle[1:-1,1:-1] + angle[2:,1:-1] + angle[1:-1,2:] + angle[2:,2:])

Mp, Lp  = angle.shape
nc.createDimension('x', Lp-2)
nc.createDimension('y', Mp-2)

nc.createVariable('ulon', 'f8', ('y', 'x'))
nc.variables['ulon'].units = 'radian'
nc.variables['ulon'][:] = ulon

nc.createVariable('ulat', 'f8', ('y', 'x'))
nc.variables['ulat'].units = 'radian'
nc.variables['ulat'][:] = ulat

nc.createVariable('angle', 'f8', ('y', 'x'))
nc.variables['angle'].units = 'radian'
nc.variables['angle'][:] = angleu

nc.createVariable('kmt', 'f8', ('y', 'x'))
nc.variables['kmt'].units = '1'
nc.variables['kmt'][:] = kmt

dx = 100.0 / pm
dy = 100.0 / pn

htn = 0.5*(dy[1:-1,1:-1] + dy[2:,1:-1])
hus = 0.5*(dx[1:-1,1:-1] + dx[1:-1,2:])
hte = 0.5*(dy[1:-1,1:-1] + dy[1:-1,2:])
huw = 0.5*(dx[1:-1,1:-1] + dx[2:,1:-1])

nc.createVariable('htn', 'f8', ('y', 'x'))
nc.variables['htn'].units = 'cm'
nc.variables['htn'][:] = htn

nc.createVariable('hus', 'f8', ('y', 'x'))
nc.variables['hus'].units = 'cm'
nc.variables['hus'][:] = hus

nc.createVariable('hte', 'f8', ('y', 'x'))
nc.variables['hte'].units = 'cm'
nc.variables['hte'][:] = hte

nc.createVariable('huw', 'f8', ('y', 'x'))
nc.variables['huw'].units = 'cm'
nc.variables['huw'][:] = huw

nc.close()
