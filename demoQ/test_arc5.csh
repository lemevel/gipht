#!/bin/csh

# test ARC function in grdmath Kurt Feigl
# 2014-08-10

\rm -f *.ps tmp.* *.grd
\rm *.cpt
gmtset -Ds # use SI units

# compile routine to make a mogi bump
# fails for unknown reasons
#gfortran m_gipht.f90 m_model.f90 t.f90 -o m_model

# run the routine to make a mogi bump
#./m_model  | awk '{print $1/1e3,$2/1e3,$3}' >! psp_mogi.xyp

# make a GMT GRD file
#gmtinfo psp_mogi.xyp -I1/1
#xyz2grd psp_mogi.xyp -R-25/25/-25/25 -Gpsp_mogi.grd -I0.1/0.1
#grdinfo psp_mogi.grd
# change units from radians to cycles
#grdmath psp_mogi.grd PI DIV 2.0 DIV = phase_cyc.grd

# make a GMT GRD file from scratch
grdmath -R-3/3/-4/4 -I0.1/0.1 X Y HYPOT 5.0 DIV 2.0 MUL PI MUL COS PI MUL  = cos.grd
grdmath -R-3/3/-4/4 -I0.1/0.1 X Y HYPOT 5.0 DIV 2.0 MUL PI MUL SIN PI MUL  = sin.grd
grdmath -R-3/3/-4/4 -I0.1/0.1 X PI MUL  = ramp.grd

# easy tests
grdmath -V sin.grd sin.grd ARC          = zero1.grd  
grdmath -V sin.grd 2.0 MUL sin.grd ARC  = nans1.grd
grdmath -V cos.grd PI      ARC          = test1.grd
grdmath -V 0 ramp.grd      ARC          = ramp1.grd

# make color table
makecpt -Crainbow -T-3.14159/3.14159/0.1 >! rad.cpt


foreach grdfile (*1.grd)

    echo "----- $grdfile -------"
    grdinfo $grdfile
   	# make an image of phase in cycles
   	# make a color table (like color_map.m)
	#grd2cpt $grdfile -V >! rad.cpt
	
	grdimage $grdfile -Crad.cpt -JX10 -P \
	-Bxf1a1+l"Easting (km)" \
	-Byf1a1+l"Northing (km)" \
	-BWSen+t"$grdfile" \
	-K >! $grdfile:r.ps
	# like colorbar
	psscale -Dx11/5/10/1 -Crad.cpt -E -Bf0.2a1":Phase:"/":Radians": -O -K -P >> $grdfile:r.ps

	grdtrack profile.xyd -G$grdfile | awk '{print $1,$4}' >! profile.dz
	
	#psxy profile.dz `gmtinfo profile.dz -I1/0.01` -JX10/5 -Y15 -Ss0.1 -G0 \
	psxy profile.dz -R-3/3/-3.5/3.5 -JX10/5 -Y15 -Ss0.1 -G0 \
	-Bxf1a1+l"Easting (km)" \
	-Byf0.1a1+l"radians" \
	-BWSen+t"$grdfile" \
	-O -P >> $grdfile:r.ps

end # loop over grid files

gs *.ps



