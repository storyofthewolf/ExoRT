pro hybrid2height, nlon, nlat,nilev,PS,P0,g,R,hyai,hybi,hyam,hybm,T,lev_Z,ilev_Z

;AUTHOR: WOLF, E.T.
;7/14/2008
;---------------------------------------------------------------------
;PURPOSE:  Convert WACCM pressure coordinates to height coordinates.
;This programs assumes that hybrid level coordinates have already been
;transformed in pressure coordinates using hybrid2pressure.pro
;
;NOTES:  May need changes of passed parameters upon implementation
;
;---------------------------------------------------------------------
;ARGUMENTS:
;nlon                => number of elements,  longitude array
;nlat                => numbrt of elements, latitude array
;nilev               => number of elements, hybrid vertical levels at interfaces
;PS[lon,lat]         => surface pressure array at each grid point
;P0                  => reference pressure
;hyai                => hybrid A coefficient at layer interfaces
;g                   => gravity [m s-2]
;R                   => gas constant dry air [J kg-1 K-1] 
;hybi                => hybrid B coefficient at layer interfaces
;hyam                => hybrid A coefficient at midlayer
;hybm                => hybrid B coefficient at midlaye
;lev_Z[lon,lat,lev]  => returns mid layer height coordinate matrix in METERS,
;                       must be defined with proper dimensions before being passed 
;ilev_Z[lon,lat,nilev]  => returns interface height coordinate matrix in METERS,
;                       must be defined with proper dimensions before being passed 
;----------------------------------------------------------------------

;define some constants to be used
ilev_Z=fltarr(nlon,nlat,nilev)                             ;height array for interface levels       

delta_z=0.0                                               ;[m]            calculated via hypsometric equation

;get surface heights from initial conidtion topographical data set
;fname='/projects/btoon/wolfet/data/atm/cam/topo/topo_wa3_4x5_1950_spinup.cam2.i.1960-01-01-00000_boville.nc'
;fname = '/Users/wolfe/IDLWorkspace/CESMfiles/topo_wa3_4x5_1950_spinup.cam2.i.1960-01-01-00000_boville.nc'
;ncid=ncdf_open(fname, /nowrite)
;ncdf_varget,ncid,'PHIS', PHIS
;ncdf_close, ncid

PHIS=fltarr(nlon,nlat)
PHIS(*,*)=0.0                 
z_surf=PHIS/g

FOR x=0, nlon-1 DO BEGIN
  FOR y=0, nlat-1 DO BEGIN
    delta_Z=0.0
    ilev_Z(x,y,nilev-1)=z_surf(x,y)                         ;add geopotential height of surface geology from initiall conditions
    FOR z=nilev-1, 1,-1 DO BEGIN
      p1=hyai(z)*P0+hybi(z)*PS(x,y)                        ;calculate lower interface pressure, ie. higher pressure/lower altitude
      p2=hyai(z-1)*P0+hybi(z-1)*PS(x,y)                    ;calculate upper interface pressure  ie lower pressure/higher altitude
      delta_Z=R*T(x,y,z-1)/g*alog(p1/p2)                   ;calculate geopotential thickness of layer
      ilev_Z(x,y,z-1)=ilev_Z(x,y,z)+delta_Z                ;calculate height of each interface based on layer thicknesses
      p_prime=hyam(z-1)*P0+hybm(z-1)*PS(x,y)                   ;calculate pressure at midpoints of layers (lev_P)
      Z_SCALE=-R*T(x,y,z-1)/g*alog(p_prime/p1)             ;interpolate midlayer heights
      lev_Z(x,y,z-1)=ilev_Z(x,y,z)+Z_SCALE
    ENDFOR
  ENDFOR
ENDFOR


end
