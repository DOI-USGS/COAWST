# writing ROMS input file

TINI="template_ini.in"

#***********************************************************************************

CODE="100"
TEMPLATE="template_ke.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="101"
TEMPLATE="template_ke.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="102"
TEMPLATE="template_ke.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="103"
TEMPLATE="template_ke.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="104"
TEMPLATE="template_ke.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="105"
TEMPLATE="template_ke.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="106"
TEMPLATE="template_kw.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="107"
TEMPLATE="template_kw.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="108"
TEMPLATE="template_kw.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="109"
TEMPLATE="template_kw.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="110"
TEMPLATE="template_kw.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="111"
TEMPLATE="template_kw.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="112"
TEMPLATE="template_gen.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="113"
TEMPLATE="template_gen.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="114"
TEMPLATE="template_gen.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="115"
TEMPLATE="template_gen.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 1400.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="116"
TEMPLATE="template_gen.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 14000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#***********************************************************************************

CODE="117"
TEMPLATE="template_gen.in"

# ----------------------------------------------------------------------------------
LEV="3p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 16              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 10.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="1p0"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 30              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p3"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 80              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
LEV="0p1"
NOMEFILE="ideal_float_test_""$CODE""_""$LEV"".in"
OUTFILE="ideal_float_test_""$CODE"

cp $TINI $NOMEFILE

echo '           N == 300              ! Number of vertical levels' >> $NOMEFILE

cat $TEMPLATE >> $NOMEFILE
 
echo '  CHARNOK_ALPHA == 56000.0d0         ! Charnok surface roughness' >> $NOMEFILE
echo '      TCLINE == 5.0d0                     ! m' >> $NOMEFILE

rst="       RSTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_rst.nc"
his="       HISNAME == Idealized/Output/""$OUTFILE""_""$LEV""_his.nc"
avg="       AVGNAME == Idealized/Output/""$OUTFILE""_""$LEV""_avg.nc"
flt="       FLTNAME == Idealized/Output/""$OUTFILE""_""$LEV""_flt.nc"

echo $rst >> $NOMEFILE
echo $his >> $NOMEFILE
echo $avg >> $NOMEFILE
echo $flt >> $NOMEFILE
#------------------------------------------------------------------------------------
#************************************************************************************
#************************************************************************************
#************************************************************************************
