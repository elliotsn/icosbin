#/bin/csh

# Script to bin Diviner modelled EFOVs from LRO orbits 14118 and 14119 onto a level 12 
# icosahedral grid. The region of interest is near the lunar north pole:
# 171 <= longitude < 176
# 83 <= latitude < 85

# Lat,lon bounds and spectral channel
set lonmin=171
set lonmax=176
set latmin=83
set latmax=85
set c=7

set ofile="c"$c"_"$lonmin"_"$lonmax"_"$latmin"_"$latmax".ofile"
set gfile="c"$c"_"$lonmin"_"$lonmax"_"$latmin"_"$latmax".gfile"
set rfile="c"$c"_"$lonmin"_"$lonmax"_"$latmin"_"$latmax".rfile"

# Triangle level of grid to create triangle numbers for.
set lev=12
set points=100

# Pass EFOVs through binning pipeline:
cat $rfile | ptrinum lat=tlat lon=tlon lev=$lev \
| ptrigather trinum1=trinum1 trinum2=trinum2 points=$points temis=orbit tinc=orbit tphase=orbit tlat=trilat tlon=trilon \
| ptrilink fields=tlat,tlon,weight weight=weight ofile=$ofile gfile=$gfile