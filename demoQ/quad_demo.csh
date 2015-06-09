#!/bin/csh

# use GMT to exercise pha2qls
# Peter Sobol & Kurt Feigl
# 2014-08-10

\rm -f *.ps tmp.* *.grd
\rm *.cpt
gmtset MEASURE_UNIT            cm

# compile routine to make a mogi bump
# fails for unknown reasons
#gfortran m_gipht.f90 m_model.f90 t.f90 -o m_model

# run the routine to make a mogi bump
#./m_model  | awk '{print $1/1e3,$2/1e3,$3}' >! psp_mogi.xyp

# make a GMT GRD file
minmax psp_mogi.xyp -I1/1
xyz2grd psp_mogi.xyp -R-25/25/-25/25 -Gpsp_mogi.grd -I0.1/0.1
grdinfo psp_mogi.grd
# change units from radians to cycles
grdmath psp_mogi.grd PI DIV 2.0 DIV = phase_cyc.grd
grdinfo phase_cyc.grd

# make an image of phase in cycles
# make a color table (like color_map.m)
#grd2cpt phase_cyc.grd -V >! cyc.cpt
makecpt -Crainbow -T-0.5/0.5/0.05 >! cyc.cpt
makecpt -Crainbow -T-128/127/4 >! dn.cpt
# show 1 sig fig
# gmtset D_FORMAT                = %.2f
# like image
grdimage phase_cyc.grd -Xc -Yc -Ccyc.cpt -JX5 -R-24/24/-24/24 -P -JX12/12 \
-Bf1a5:"Easting (km)":/f1a5:"Northing (km)":WSne -K >! phase_cyc.ps
# like colorbar
psscale -D12/6/12/1 -Ccyc.cpt -Bf0.1a0.1":Phase:"/":Cycles": -O >> phase_cyc.ps 
#gv phase_cyc.ps &

grdimage phase_cyc.grd -X3 -Y21 -Ccyc.cpt -JX5 -R-24/24/-24/24 -P -JX06/06 \
-Bf1a10:"Original":/f1a10:::.:WSne -K -0 >! summary.ps

psscale -D06/3/06/0.5 -Ccyc.cpt -Bf0.2a0.2":Phase:"/":Cycles": -O -K >> summary.ps 


# make phase file in DN
grdmath psp_mogi.grd PI DIV 2.0 DIV 256.0 MUL = psp_test.grd
grdinfo psp_test.grd 
# grd2pha is home-made script 
#/usr1/feigl/diapason4/dtools/grd2pha.gmt
#grd2pha.gmt psp_test.grd  # writes psp_test.pha
  # output so that first line is at the NORTHERN EDGE
  #grd2xyz tmp.grd -ZcTL >! output.pha
grd2xyz psp_test.grd -ZcTL >! psp_test.pha

# replace with grd2xyz with binary option -Z
#        -Z sets exact specification of resulting 1-column output z-table
#           If data is in row format, state if first row is at T(op) or B(ottom)
#             Then, append L or R to indicate starting point in row
#           If data is in column format, state if first columns is L(left) or R(ight)
#             Then, append T or B to indicate starting point in column
#           Append x if gridline-registered, periodic data in x without repeating column at xmax
#           Append y if gridline-registered, periodic data in y without repeating row at ymax
#          Specify one of the following data types (all binary except a):
#             a  Ascii
#             c  signed 1-byte character **** .pha files are like this 1 DN = 2^-8 values range from -2^7 to +2^7 - 1 [-128, +127]
#             u  unsigned 1-byte character
#             h  signed short 2-byte integer  **** .i2 files are like this  1 DN = 2^-16 values range from -2^15 to +2^15 - 1 [-32768, +32768]
#             H  unsigned short 2-byte integer
#             i  signed 4-byte integer
#             I  unsigned 4-byte integer
#             l  long (4- or 8-byte) integer
#             f  4-byte floating point single precision
#             d  8-byte floating point double precision
 #          [Default format is scanline orientation in ascii representation: -ZTLa]
histogram.gmt psp_test.pha
#exit


# PERFORM QUAD-TREE COMPRESSION
#~feigl/GIPhT1.8/src/pha2qls3.a64 psp_test.pha 501 501 qsp_test.pha grx_test.i2 gry_test.i2 quad.qls 12 4 64
#~feigl/GIPhT1.8/src/pha2qls3.a64 psp_test.pha 501 501 qsp_test.pha grx_test.i2 gry_test.i2 quad.qls 12 12 64
echo "###Calculating QLS###"
../src/pha2qls.maci64 psp_test.pha 501 501 -P qsp_test.pha -G -L 16 -N 4 -M 32

echo "###Plot Quad Tree Phase###"
# plot quad-tree phase
xyz2grd qsp_test.pha -ZcTL -R-25/25/-25/25 -Gqsp_test.grd -I0.1/0.1
grdmath qsp_test.grd 256 DIV = qsp_test_cyc.grd
gmtset D_FORMAT                = %.2f
grdimage qsp_test_cyc.grd  -Ccyc.cpt  -R-24/24/-24/24 -JX12/12 -U/2/-4/$0 \
-Bf1a5:"Easting (km)":/f1a5:"After Quadtree"::.:wSne -K >! qsp_test.ps
psscale -D12/6/12/1 -Ccyc.cpt -Bf0.1a0.1":Phase:"/"::" -O >> qsp_test.ps 
#gv qsp_test.ps  &

grdimage qsp_test_cyc.grd -X7 -Y0  -Ccyc.cpt  -R-24/24/-24/24 -JX6/6            \
-Bf1a10:"After Quadtree (flat)":/f1a5:::.:wSne -K -O -P >> summary.ps
psscale -D6/3/6/.5 -Ccyc.cpt -Bf0.2a0.2":Phase:"/"::" -O -K >> summary.ps 


echo "###Plot East Gradient###"
# plot  East component of gradient 
xyz2grd grx_test.i2 -ZhTL -R-25/25/-25/25 -Ggrx_test.grd -I0.1/0.1
grdmath grx_test.grd PI MUL 2.0 MUL 256 DIV = grx_test_rad.grd
grdmath grx_test.grd 256 DIV 256 DIV        = grx_test_cyc.grd
gmtset D_FORMAT                = %.3f
makecpt -T-0.01/0.01/0.001 -D > xgrad.cpt
#grd2cpt grx_test_cyc.grd -V -L-0.1/0.1 >! xgrad.cpt


gmtset D_FORMAT                = %.2f
grdimage grx_test_cyc.grd -Cxgrad.cpt -P -Xc -Yc -R-24/24/-24/24 -JX12/12 \
-Bf1a5:"Easting (km)":/f1a5:"Northing (km)":wSne -K  >! xgrad.ps
gmtset D_FORMAT                = %.3f
psscale -D12/6/12/1 -Cxgrad.cpt -Bf0.001a0.001":Gradient:"/"::" -O >> xgrad.ps 
#gv xgrad..ps &

grdimage grx_test_cyc.grd -Cxgrad.cpt -P -X-7 -Y-9 -R-24/24/-24/24 -JX06/06 \
-Bf1a10:"X-Gradient":/f1a5::wSne -K -O >> summary.ps
psscale -D06/3/06/.5 -Cxgrad.cpt -Bf0.001a0.001":Gradient:"/"::" -K -O >> summary.ps 


echo "###Plot North Gradient###"
# plot  North component of gradient 
xyz2grd gry_test.i2 -ZhTL -R-25/25/-25/25 -Ggry_test.grd -I0.1/0.1
grdmath gry_test.grd PI MUL 2.0 MUL 256 DIV = gry_test_rad.grd
grdmath gry_test.grd 256 DIV 256 DIV        = gry_test_cyc.grd
gmtset D_FORMAT                = %.4f
makecpt -T-0.01/0.01/0.0005 -D > ygrad.cpt
#grd2cpt gry_test_cyc.grd -V -L-0.1/0.1 >! ygrad.cpt


gmtset D_FORMAT                = %.2f
grdimage gry_test_cyc.grd -Cygrad.cpt -P -Xc -Yc \
-R-24/24/-24/24 -JX12/12 \
-Bf1a5:"Easting (km)":/f1a5:"Northing (km)":wSne -K  >! ygrad.ps
gmtset D_FORMAT                = %.2f
psscale -D12/6/12/1 -Cygrad.cpt -Bf0.001a0.001":Gradient:"/"::" -O >> ygrad.ps 
#gv ygrad.ps &
grdimage gry_test_cyc.grd -Cygrad.cpt -P -X7  -Y0   -R-24/24/-24/24 -JX06/06 \
-Bf1a10:"Y-Gradient":/f1a10::wSne -K -O >> summary.ps
psscale -D06/3/06/.5 -Cxgrad.cpt -Bf0.002a0.002":Gradient:"/"::" -K -O >> summary.ps 


echo "###Reconstruct form QLS data###"
# plot reconstruction
../src/qls2pha.maci64 psp_test.qls -o qha_test.pha

echo "###Plot Reconstruction###"

xyz2grd qha_test.pha -ZcTL -R-25/25/-25/25 -Gqha_test.grd -I0.1/0.1
grdmath qha_test.grd 256 DIV = qha_test_cyc.grd
gmtset D_FORMAT                = %.2f
grdimage qha_test_cyc.grd -X3 -Y5 -Ccyc.cpt  -R-24/24/-24/24 -P -JX12/12 -U/2/-2/$0 \
-Bf1a5:"Easting (km)":/f1a5:"Northing (km)"::."Reconstruction from QLS File":wSne -K >! qha_test.ps
psscale -D12/6/12/1 -Ccyc.cpt -Bf0.1a0.1":Phase:"/"::" -O >> qha_test.ps 
#gv qha_test.ps &

grdimage qha_test_cyc.grd -X-7 -Y-9 -Ccyc.cpt  -R-24/24/-24/24 -P -JX06/06 \
-Bf1a10:"Reconstruction from QLS":/f1a10:::.:wSne -K  -O  >> summary.ps
psscale -D6/3/6/.5 -Ccyc.cpt -Bf0.1a0.1":Phase:"/"::" -K -O  >> summary.ps



# plot reconstruction error
../src/diffpha.maci64 psp_test.pha qha_test.pha  501 501 dif_test.pha -0
echo "###Plot reconstruction error###"
xyz2grd dif_test.pha -ZcTL -R-25/25/-25/25 -Gdif_test.grd -I0.1/0.1


#grdmath qsp_test.grd qha_test.grd  SUB  = dif_test.grd
#grdmath dif_test.grd 128 ADD 256 MOD 128 SUB = dif_test.grd
grdmath dif_test.grd 256 DIV = dif_test_cyc.grd
grdimage dif_test_cyc.grd -X3 -Y5 -Ccyc.cpt  -R-24/24/-24/24 -P -JX12/12 -U/2/-2/$0 \
-Bf1a5:"Original-Reconstruction":/f1a5:::.:wSne -K >! dif_test.ps
psscale -D12/6/12/1 -Ccyc.cpt -Bf0.1a0.1":Phase:"/"::" -O >> dif_test.ps 
#gv dif_test.ps &

grdimage dif_test_cyc.grd -X7 -Y0 -Ccyc.cpt  -R-24/24/-24/24 -P -JX06/06 \
-Bf1a10:"Original-Reconstruction":/f1a10:::.:wSne -K -O  >> summary.ps
psscale -D06/3/06/.5 -Ccyc.cpt -Bf0.2a0.2":Phase:"/"::" -O >> summary.ps 




ls -l *.ps
echo gs *.ps

#gv summary.ps -media=LETTER &
