#!/bin/bash
#
# svn $Id$
#######################################################################
## Copyright (c) 2002-2026 The ROMS Group                             #
##   Licensed under a MIT/X style license                             #
##   See License_ROMS.md                                              #
################################################## Hernan G. Arango ###
##                                                                    #
## ROMS Split RBL4D-Var Data Assimilation running BASH script:        #
##                                                                    #
## The RBL4D-Var is split in several phases and executables:          #
##                                                                    #
##   1) background                          ROMS_EXE_A                #
##   2) increment                           ROMS_EXE_B                #
##   3) analysis                            ROMS_EXE_A                #
##                                                                    #
## The 4D-Var algorithm is split into Executable A and B:             #
##                                                                    #
## ROMS_EXE_A: Computes ROMS NLM trajectory in the "background" phase #
##              to linearize the TLM and ADM kernels used in the      #
## 4D-Var minimization. It interpolates the NLM solution to the       #
## observation locations in space and time, H(x). The NLM  is run at  #
## a higher resolution and single or double precision. It could be a  #
## part of a coupling system and include nested grids.                #
##                                                                    #
## ROMS_EXE_B: It is used in the "increment" phase. The 4D-Var        #
##             increment is obtained by minimizing the cost function  #
## over Ninner loops. It is possible to use a coarser grid resolution #
## in the iterations of the inner loops. If so, the NLM trajectory    #
## was decimated by a factor of two in the "background" phase. Then,  #
## after the inner loops, the coarser grid increment must be          #
## interpolated to the finer grid before the "analysis" phase. The    #
## "increment" phase is run at a lower resolution to accelerate the   #
## computations. It could be run in single precision.                 #
##                                                                    #
## Notice that this script has a customizable function, My4DVarScript,#
## at the top that creates "rbl4dvar_*.in" from a template, which is  #
## overwritten in each phase. The "rbl4dvar_*.in" passes the phase    #
## and outer loop value to ROMS to compute.                           #
##                                                                    #
## Then, there is the customizable section for the computer batch     #
## directives and the tunable parameters.                             #
##                                                                    #
## RBL4D-Var phases workflow:                                         #
##                                                                    #
##   CALL prior_error                                                 #
##                                                                    #
##   CALL background (outer=0, RunInterval)                           #
##                                                                    #
##   OuterLoop : DO outer=1:Nouter                                    #
##     CALL increment (outer, RunInterval)               Inner-loops  #
##     CALL analysis (outer, RunInterval)                             #
##   END DO OuterLoop                                                 #
##                                                                    #
##   CALL posterior_error (RunInterval)                  if requested #
##                                                                    #
## Check Test Case at:                                                #
##                                                                    #
## https://github.com/myroms/roms_test/tree/main/USEC/RBL4DVAR_mixres #
##                                                                    #
#######################################################################

##  SLURM configuration for amarel:
##
##  Use 'sbatch submit_split_rbl4dvar.sh'  to queue the job
##  Use 'squeue -p p_omg_1'                to check our group jobs (including JOBID)
##  Use 'sacct -j JOBID -l'                to check job accounting data
##  Use 'scontrol show job JOBID'          to check job configuration
##  Use 'sinfo JOBID'                      to check job information
##  Use 'scancel JOBID'                    to cancel a job

#SBATCH --exclusive                      # don't run on nodes with other jobs running
#SBATCH --partition=p_omg_1              # Partition (job queue), NodeList: 108,116-120
#SBATCH --requeue                        # Return job to the queue if preempted
#SBATCH --job-name=ROMS_mixres_rb4dvar   # Assign an short name to your job
#SBATCH --nodes=1                        # Number of nodes you require (each has 32 PETs)
#SBATCH --ntasks=12                      # Total number of tasks you'll launch
#SBATCH --ntasks-per-node=32             # Number of tasks you'll launch on each node
#SBATCH --cpus-per-task=1                # Cores per task (>1 if multithread tasks)
#SBATCH --mem=177000                     # Real memory (RAM) required (MB)
#SBATCH --time=00-01:00:00               # Total run time limit (DD-HH:MM:SS)
#SBATCH --output=log.%N.%j               # STDOUT output file
#SBATCH --error=err.%N.%j                # STDERR output file (optional)
#SBATCH --export=ALL                     # Export you current env to the job env

##---------------------------------------------------------------------
## Control switches: What do you want to do?
##---------------------------------------------------------------------

#      DRYRUN=1                # Print configuration but do not execute
       DRYRUN=0                # Run 4D-Var cycle

        BATCH=0                # No batch system submission
#       BATCH=1                # Use batch system SLURM to submit:
                               #   "sbatch submit_mixres_rbl4dvar.sh"

    ROMS_ROOT=${HOME}/ocean/repository/git/roms    # Set ROMS location:
                                                   # uses ROMS/Bin
                                                   # Perl scripts

##---------------------------------------------------------------------
## USER TUNABLE PARAMETERS. If you follow recommendations, this is
## the only section that you need to customize..
##---------------------------------------------------------------------

      ROMS_APP="USEC"                # ROMS Application CPP

      roms_app=`echo ${ROMS_APP} | tr '[:upper:]' '[:lower:]'` # lowercase

    roms_app_F=${roms_app}"3km"        # outer loop grid
    roms_app_C=${roms_app}"6km"        # inner loop grid

          ResF=3                       # outer loop resolution
          ResC=6                       # innet loop resolution

       HereDir=${PWD}                  # current directory

       DataDir="../../Data"            # data directory 

        ObsDir=${DataDir}/OBS          # observations directory

 if [[ "$OSTYPE" == "darwin"* ]]; then
      DATE_EXE=gdate                   # macOS system GNU date
 else
      DATE_EXE=date                    # Linux system date
 fi

 if [ ${BATCH} -eq 1 ]; then
          SRUN="srun --mpi=pmi2"       # SLURM workload manager
 else    
        MPIrun="mpirun -np"            # Basic MPI workload manager
 fi

    ROMS_EXE_A="romsM_nl"              # ROMS executable A
    ROMS_EXE_B="romsM_da"              # ROMS executable B

    START_DATE="2019-08-27"            # 4D-Var starting date

  FRST_INI_DAY="${START_DATE}"         # first cycle initialization date
# FRST_INI_DAY="2019-08-27"            # restart initialization date

# LAST_INI_DAY="${START_DATE}"         # last cycle initialization date
  LAST_INI_DAY="2019-08-30"            # last cycle initialization date

  ROMS_TIMEREF="2006-01-01"            # ROMS time reference date

      INTERVAL=3                       # 4D-Var interval window (days)

        nPETsX=3                       # number PETs in the X-direction
        nPETsY=4                       # number PETs in the Y-direction

      MyNouter=1                       # number of 4D-Var outer loops
#     MyNouter=2                       # number of 4D-Var outer loops

      MyNinner=16                      # number of 4D-Var inner loops: RPCG
#     MyNinner=8                       # number of 4D-Var inner loops: RPCG

     MyTimeIAU=0.0d0                   # Incremental Analysis Update window (days)
#    MyTimeIAU=0.03125d0               # Incremental Analysis Update window (days)

     MyNHIS_nl=30                      # NLM trajectory is saved hourly
  MyNDEFHIS_nl=0                       # No multi-file NLM trajectory

     MyNXTR_nl=30                      # NLM decimation is saved hourly
  MyNDEFXTR_nl=0                       # No multi-file NLM decimation

     MyNQCK_nl=30                      # NLM quicksave trajectory is saved hourly
  MyNDEFQCK_nl=0                       # No multi-file NLM trajectory
 
     MyNTLM_nl=30                      # TLM trajectory is saved hourly
     MyNADJ_nl=2160                    # strong contraint, ADM saved at end
     MyNSFF_nl=30                      # SFF asjustment is saved hourly
     MyNOBC_nl=30                      # OBC asjustment is saved hourly

     MyINP_LIB=2                       # reading library: [1] standard [2] PIO
     MyOUT_LIB=2                       # writing library: [1] standard [2] PIO
  MyPIO_METHOD=3                       # [1] NetCDF3, ...
 MyPIO_IOTASKS=2                       # number of I/O processes
  MyPIO_STRIDE=5                       # stride in MPI-rank between I/O tasks
    MyPIO_BASE=0                       # offset for the first I/O task
   MyPIO_REARR=1                       # rearranger method: [1] box [2] subset
MyPIO_REARRCOM=1                       # rearranger communications: [0] p2p [1] coll
MyPIO_REARRDIR=0                       # rearranger direction: [0] I2C/C2I, ... [3]

       restart=0                       # restart 4D-Var cycle (0:no, 1:yes)
#      restart=1                       # restart 4D-Var cycle (0:no, 1:yes)

#   ROMS_NLpre="roms_nl_${roms_app}_era5"          # ROMS NLM stdinp prefix, ERA5
    ROMS_NLpre="roms_nl_${roms_app}_nam"           # ROMS NLM stdinp prefix, NAM
    ROMS_NLtmp="${ROMS_NLpre}.tmpl"                # ROMS NLM stdinp template

#   ROMS_DApre="roms_da_${roms_app}_era5"          # ROMS ADM/TLM stdinp prefix, ERA5
    ROMS_DApre="roms_da_${roms_app}_nam"           # ROMS ADM/TLM stdinp prefix, NAM
    ROMS_DAtmp="${ROMS_DApre}.tmpl"                # ROMS ADM/TLM stdinp template

 if [[ ${restart} -eq 0 ]]; then
      ROMSiniF="usec3km_roms_ini_20190827.nc4"     # ROMS outer loop IC
   ROMSgeniniC="usec6km_roms_ini.nc4"              # ROMS Generic inner loop IC
    ROMSiniDir=${DataDir}/INI                      # ROMS IC directory
 else
       ROMSini="usec3km_roms_dai_20190827.nc"      # restart ROMS IC
    ROMSiniDir=../2019.08.27                       # restart ROMS IC directory
 fi

   Inp4DVAR_nl="rbl4dvar_nl.in"                    # ROMS NL 4D-Var input script
   Inp4DVAR_da="rbl4dvar_da.in"                    # ROMS DA 4D-Var input script

 if [[ ${BATCH} -eq 1 ]]; then
         etime=00:00:00                            # elapsed time
         ptime=${etime}                            # previous time
 else
         stime=`${DATE_EXE} -u +"%s"`              # start time (sec) since epoch
         ptime=$stime                              # previous time
 fi

## Inner loops grid parameters:


 if [[ $ResF == $ResC ]]; then
     MyNHIS_da=${MyNHIS_nl}                        # same resolution for outer
  MyNDEFHIS_da=${MyNDEFHIS_nl}                     # and inner loops

     MyNXTR_da=${MyNXTR_nl}
  MyNDEFXTR_da=${MyNDEFXTR_nl}

     MyNQCK_da=${MyNQCK_nl}
  MyNDEFQCK_da=${MyNDEFQCK_nl}

     MyNTLM_da=${MyNTLM_nl}
     MyNADJ_da=${MyNADJ_nl}
     MyNSFF_da=${MyNSFF_nl}
     MyNOBC_da=${MyNOBC_nl}
 else
     MyNHIS_da=$(( ${MyNHIS_nl} / 2 ))            # finer   loops grid
  MyNDEFHIS_da=$(( ${MyNDEFHIS_nl} / 2 ))         # coarser inner loops
                                                  # divide by 2
     MyNXTR_da=$(( ${MyNXTR_nl} / 2 ))
  MyNDEFXTR_da=$(( ${MyNDEFXTR_nl} / 2 ))

     MyNQCK_da=$(( ${MyNQCK_nl} / 2 ))
  MyNDEFQCK_da=$(( ${MyNDEFQCK_nl} / 2 ))

     MyNTLM_da=$(( ${MyNTLM_nl} / 2 ))
     MyNADJ_da=$(( ${MyNADJ_nl} / 2 ))
     MyNSFF_da=$(( ${MyNSFF_nl} / 2 ))
     MyNOBC_da=$(( ${MyNOBC_nl} / 2 ))
 fi

##---------------------------------------------------------------------
## END OF USER TUNABLE PARAMETERS.
##---------------------------------------------------------------------

########################################################################
## RBL4D-Var data assimilation input script FUNCTION.  It generates    #
## 'rbl4dvar_nl.in' or 'rbl4dvar_nl.in' from the 's4dvar.in' template. #
########################################################################

## Start of My4DVarScript() FUNCTION definition

My4DVarScript() {

     DataDir=${1}              # Data directory
  SUBSTITUTE=${2}              # ROMS Perl substitution function
   OuterLoop=${3}              # current outer loop counter
  Phase4DVAR=${4}              # current 4D-Var computation phase
     TimeIAU=${5}              # Incremental Analysis Update window
     OBSname=${6}              # 4D-Var observations NetCDF file
     Fprefix=${7}              # ROMS output files prefix (use roms_app)
     Fsuffix=${8}              # ROMS output files suffix
    Inp4DVAR=${9}              # 4D-Var standard input
         res=${10}             # Grid resolution

  echo
  echo "   Creating 4D-Var Input Script from Template: ${Inp4DVAR}" \
       "  Outer = ${OuterLoop}  Phase = ${Phase4DVAR}"
  echo "     (Resolution = ${res} km, Fprefix = ${Fprefix}, Fsuffix = ${Fsuffix})"
  echo

## Set model, initial conditions, boundary conditions and surface
## forcing error covariance standard deviations files. For model
## error, use the same as initial condition since we are not
## running in weak-constraint mode.

if [[ $res -eq 6 ]]; then
 STDnameM=${DataDir}/STD/usec6km_roms_std_i_${Fsuffix}.nc4
 STDnameI=${DataDir}/STD/usec6km_roms_std_i_${Fsuffix}.nc4
 STDnameB=${DataDir}/STD/usec6km_roms_std_b_${Fsuffix}.nc4
 STDnameF=${DataDir}/STD/usec6km_roms_std_f_${Fsuffix}.nc4
else
 STDnameM=${DataDir}/STD/usec3km_roms_std_i_${Fsuffix}.nc4
 STDnameI=${DataDir}/STD/usec3km_roms_std_i_${Fsuffix}.nc4
 STDnameB=${DataDir}/STD/usec3km_roms_std_b_${Fsuffix}.nc4
 STDnameF=${DataDir}/STD/usec3km_roms_std_f_${Fsuffix}.nc4
fi

## Set output file for standard deviation computed/modeled from background
## (prior) state.

if [[ $res -eq 6 ]]; then
 STDnameC=usec6km_roms_std_computed.nc4
else
 STDnameC=usec3km_roms_std_computed.nc4
fi

## Set model, initial conditions, boundary conditions and surface
## forcing error covariance normalization factors files. For model
## error, use the same as initial condition since we are not
## running in weak-constraint mode.

if [[ $res -eq 6 ]]; then
 NRMnameM=${DataDir}/NRM/usec6km_roms_nrm_i.nc4
 NRMnameI=${DataDir}/NRM/usec6km_roms_nrm_i.nc4
 NRMnameB=${DataDir}/NRM/usec6km_roms_nrm_b.nc4
 NRMnameF=${DataDir}/NRM/usec6km_roms_nrm_f.nc4
else
 NRMnameM=${DataDir}/NRM/usec3km_roms_nrm_i.nc4
 NRMnameI=${DataDir}/NRM/usec3km_roms_nrm_i.nc4
 NRMnameB=${DataDir}/NRM/usec3km_roms_nrm_b.nc4
 NRMnameF=${DataDir}/NRM/usec3km_roms_nrm_f.nc4
fi

## Modify 4D-Var template input script and specify above files.

 if [ -f $Inp4DVAR ]; then
   /bin/rm ${Inp4DVAR}
 fi

 cp ../s4dvar.in ${Inp4DVAR}

 $SUBSTITUTE $Inp4DVAR MyOuterLoop   ${OuterLoop}
 $SUBSTITUTE $Inp4DVAR MyPhase4DVAR  ${Phase4DVAR}
 $SUBSTITUTE $Inp4DVAR MyTimeIAU     ${TimeIAU}
 $SUBSTITUTE $Inp4DVAR roms_std_i.nc ${STDnameI}
 $SUBSTITUTE $Inp4DVAR roms_std_m.nc ${STDnameM}
 $SUBSTITUTE $Inp4DVAR roms_std_b.nc ${STDnameB}
 $SUBSTITUTE $Inp4DVAR roms_std_f.nc ${STDnameF}
 $SUBSTITUTE $Inp4DVAR roms_std_c.nc ${STDnameC}
 $SUBSTITUTE $Inp4DVAR roms_nrm_i.nc ${NRMnameI}
 $SUBSTITUTE $Inp4DVAR roms_nrm_m.nc ${NRMnameM}
 $SUBSTITUTE $Inp4DVAR roms_nrm_b.nc ${NRMnameB}
 $SUBSTITUTE $Inp4DVAR roms_nrm_f.nc ${NRMnameF}
 $SUBSTITUTE $Inp4DVAR roms_obs.nc   ${OBSname}
 $SUBSTITUTE $Inp4DVAR roms_hss.nc   ${Fprefix}_roms_hss_${Fsuffix}.nc
 $SUBSTITUTE $Inp4DVAR roms_lcz.nc   ${Fprefix}_roms_lcz_${Fsuffix}.nc
 $SUBSTITUTE $Inp4DVAR roms_lze.nc   ${Fprefix}_roms_lze_${Fsuffix}.nc

 if [[ $res -eq 6 ]]; then           # same as outer loops grid
   $SUBSTITUTE $Inp4DVAR roms_mod.nc usec3km_roms_mod_${Fsuffix}.nc
   $SUBSTITUTE $Inp4DVAR roms_inc.nc ${Fprefix}_roms_inc_${Fsuffix}.nc
 else
   $SUBSTITUTE $Inp4DVAR roms_mod.nc ${Fprefix}_roms_mod_${Fsuffix}.nc
   $SUBSTITUTE $Inp4DVAR roms_inc.nc usec6km_roms_itl_${Fsuffix}.nc
 fi

 $SUBSTITUTE $Inp4DVAR roms_err.nc   ${Fprefix}_roms_err_${Fsuffix}.nc
}

## End of My4DVarScript() FUNCTION definition

#######################################################################
## Main body of script starts here. It is very unlikely that the USER
## needs to modify it.
#######################################################################

   SUBSTITUTE=${ROMS_ROOT}/ROMS/Bin/substitute   # Perl substitution
   separator1=`perl -e "print ':' x 100;"`       # title sparator
   separator2=`perl -e "print '-' x 100;"`       # run sparator

        nPETs=$(( $nPETsX * $nPETsY ))

echo
echo "${separator1}"
echo " ROMS Split RBL4D-Var Data Assimilation: ${ROMS_APP}"
echo "${separator1}"
echo
echo "                    ROMS Root: ${ROMS_ROOT}"
echo "            ROMS Executable A: ${ROMS_EXE_A}  (background, analysis)"
echo "            ROMS Executable B: ${ROMS_EXE_B}  (increment, post_error)"
echo

##---------------------------------------------------------------------
## Compute date number for reference date, and first and last
## initialization dates. It uses ROMS/Bin/dates Perl script.
##---------------------------------------------------------------------

         S_DN=`${ROMS_ROOT}/ROMS/Bin/dates datenum ${START_DATE}`
         F_DN=`${ROMS_ROOT}/ROMS/Bin/dates datenum ${FRST_INI_DAY}`
         L_DN=`${ROMS_ROOT}/ROMS/Bin/dates datenum ${LAST_INI_DAY}`
       REF_DN=`${ROMS_ROOT}/ROMS/Bin/dates datenum ${ROMS_TIMEREF}`

echo "      RBL4D-Var Starting Date: ${START_DATE}  datenum = ${S_DN}"
echo "   First RBL4D-Var Cycle Date: ${FRST_INI_DAY}  datenum = ${F_DN}"
echo "    Last RBL4D-Var Cycle Date: ${LAST_INI_DAY}  datenum = ${L_DN}"
echo "          ROMS Reference Date: ${ROMS_TIMEREF}  datenum = ${REF_DN}"
echo "       RBL4D-Var Cycle Window: ${INTERVAL} days"
echo "      Number of parallel PETs: ${nPETs}  (${nPETsX}x${nPETsY})"
echo "   Current Starting Directory: ${HereDir}"
echo "         ROMS Application CPP: ${ROMS_APP}"
echo "      Descriptor in filenames: ${roms_app_F}  (outer loops grid)"
echo "      Descriptor in filenames: ${roms_app_C}  (inner loops grid)"
echo

##=====================================================================
## Loop over all 4D-Var Cycles: FRST_INI_DAY to LAST_INI_DAY
##=====================================================================

if [ ${restart} -eq 0 ]; then
  Cycle=0                               # 4D-Var cycle counter
else
  EDAYS=`${ROMS_ROOT}/ROMS/Bin/dates daysdiff ${START_DATE} ${FRST_INI_DAY}`

  let "Cycle=${EDAYS} / ${INTERVAL}"    # restart cycle counter
fi

SDAY=${F_DN}                            # initialize starting day

while [ $SDAY -le $L_DN ]; do

## Set configuration parameters.

          Cycle=$(( ${Cycle} + 1 ))          # advance cycle by one
         DSTART=$(( ${SDAY} - ${REF_DN} ))   # ROMS DSTART parameter
           EDAY=$(( ${SDAY} + ${INTERVAL} )) # end day for 4D-Var cycle

          RDATE=`${ROMS_ROOT}/ROMS/Bin/dates numdate ${REF_DN}`
          SDATE=`${ROMS_ROOT}/ROMS/Bin/dates numdate ${SDAY}`
          EDATE=`${ROMS_ROOT}/ROMS/Bin/dates numdate ${EDAY}`
            DOY=`${ROMS_ROOT}/ROMS/Bin/dates yday ${SDATE}`
           yday=`printf %03d $DOY`

  ReferenceTime=`${DATE_EXE} -d "${RDATE}" '+%Y %m %d %H %M %S'`
      StartTime=`${DATE_EXE} -d "${SDATE}" '+%Y %m %d %H %M %S'`
       StopTime=`${DATE_EXE} -d "${EDATE}" '+%Y %m %d %H %M %S'`
       FprefixF="${roms_app_F}"              # outer loops grid
       FprefixC="${roms_app_C}"              # inner loopd grid
        Fsuffix=`${DATE_EXE} -d "${SDATE}" '+%Y%m%d'`
         RunDir=`${DATE_EXE} -d "${SDATE}" '+%Y.%m.%d'`
     ROMS_INI_F="${ROMSiniF}"
     ROMS_INI_C="${FprefixC}_roms_ini_${Fsuffix}.nc4"

  echo "${separator2}"
  echo
  echo "       RBL4D-Var Cycle Date: ${SDATE}  DayOfYear = ${yday}" \
                                     " Cycle = ${Cycle}"
  echo "         Data sub-directory: ${DataDir}"
  echo "          Run sub-directory: ${RunDir}"
  echo "      Number of outer loops: ${MyNouter}"
  echo "      Number of inner loops: ${MyNinner}"
  echo "    NLM trajectory  writing: ${MyNHIS_nl}, ${MyNHIS_da} timesteps"
  echo "    NLM quicksave   writing: ${MyNQCK_nl}, ${MyNQCK_da} timesteps"
  echo "    NLM decimation  writing: ${MyNXTR_nl}, ${MyNXTR_da} timesteps"
  echo "    TLM trajectory  writing: ${MyNTLM_nl}, ${MyNTLM_da} timesteps"
  echo "    ADM trajectory  writing: ${MyNADJ_nl}, ${MyNADJ_da} timesteps"
  echo "    SFF adjustment  writing: ${MyNSFF_nl}, ${MyNSFF_da} timesteps"
  echo "    OBC adjustment  writing: ${MyNOBC_nl}, ${MyNOBC_da} timesteps"
  echo "  NLM multi-file trajectory: ${MyNDEFHIS_nl}, ${MyNDEFHIS_da} timesteps"
  echo "                ROMS DSTART: ${DSTART}.0d0"
  echo "     Outer Loops Resolution: ${ResF} km"
  echo "     Inner Loops Resolution: ${ResC} km"
  echo "           I/O Files Prefix: ${FprefixF}  (outer loops grid)"
  echo "           I/O Files Prefix: ${FprefixC}  (inner loops grid)"
  echo "           I/O Files Suffix: ${Fsuffix}"
  echo "              ReferenceTime: ${ReferenceTime}"
  echo "        RBL4D-Var StartTime: ${StartTime}"
  echo "        RBL4D-Var  StopTime: ${StopTime}"
  echo "   ROMS outer loops Grid IC: ${ROMSiniDir}/${ROMS_INI_F}"
  echo "   ROMS inner loops Grid IC: ${ROMSiniDir}/${ROMS_INI_C}"

##---------------------------------------------------------------------
## Create run sub-directory based on cycle date (YYYY.MM.DD) and create
## ROMS standard input script from template.
##---------------------------------------------------------------------

  ROMS_NLinp=`echo ${ROMS_NLpre}_${Fsuffix}'.in'`
  ROMS_DAinp=`echo ${ROMS_DApre}_${Fsuffix}'.in'`

  echo "   NL Standard Input Script: ${ROMS_NLinp}   (outer loops grid)"
  echo "   DA Standard Input Script: ${ROMS_DAinp}   (inner loops grid)"
  echo "  NL RBL4D-Var Input Script: ${Inp4DVAR_nl}  (outer loops grid)"
  echo "  DA RBL4D-Var Input Script: ${Inp4DVAR_da}  (inner loops grid)"
  echo

  if [[ ${DRYRUN} -eq 1 ]]; then              # if dry-run, remove run
    if [ -d ${RunDir} ]; then                 # sub-directory if exist
      /bin/rm -rf ${RunDir}
    fi
  fi

  if [ ! -d ./${RunDir} ]; then
    mkdir ${RunDir}
    echo "Cycle ${Cycle}, Creating run sub-directory: ${RunDir}"
    echo
  fi

  echo "Changing to directory: ${HereDir}/${RunDir}"
  echo

  cd ${RunDir}

  echo "   Creating NL ROMS Standard Input Script: ${ROMS_NLinp}"

  if [ -f ${ROMS_NLinp} ]; then
    /bin/rm ${ROMS_NLinp}
  fi
  cp -f ../${ROMS_NLtmp} ${ROMS_NLinp}

  $SUBSTITUTE ${ROMS_NLinp} MyNtileI       ${nPETsX}
  $SUBSTITUTE ${ROMS_NLinp} MyNtileJ       ${nPETsY}
  $SUBSTITUTE ${ROMS_NLinp} MyNouter       ${MyNouter}
  $SUBSTITUTE ${ROMS_NLinp} MyNinner       ${MyNinner}
  $SUBSTITUTE ${ROMS_NLinp} MyNHIS         ${MyNHIS_nl}
  $SUBSTITUTE ${ROMS_NLinp} MyNDEFHIS      ${MyNDEFHIS_nl}
  $SUBSTITUTE ${ROMS_NLinp} MyNXTR         ${MyNXTR_nl}
  $SUBSTITUTE ${ROMS_NLinp} MyNDEFXTR      ${MyNDEFXTR_nl}
  $SUBSTITUTE ${ROMS_NLinp} MyNQCK         ${MyNQCK_nl}
  $SUBSTITUTE ${ROMS_NLinp} MyNDEFQCK      ${MyNDEFQCK_nl}
  $SUBSTITUTE ${ROMS_NLinp} MyNTLM         ${MyNTLM_nl}
  $SUBSTITUTE ${ROMS_NLinp} MyNADJ         ${MyNADJ_nl}
  $SUBSTITUTE ${ROMS_NLinp} MyNSFF         ${MyNSFF_nl}
  $SUBSTITUTE ${ROMS_NLinp} MyNOBC         ${MyNOBC_nl}
  $SUBSTITUTE ${ROMS_NLinp} MyDSTART       "${DSTART}.0d0"
  $SUBSTITUTE ${ROMS_NLinp} MyINP_LIB      ${MyINP_LIB}
  $SUBSTITUTE ${ROMS_NLinp} MyOUT_LIB      ${MyOUT_LIB}
  $SUBSTITUTE ${ROMS_NLinp} MyPIO_METHOD   ${MyPIO_METHOD}
  $SUBSTITUTE ${ROMS_NLinp} MyPIO_IOTASKS  ${MyPIO_IOTASKS}
  $SUBSTITUTE ${ROMS_NLinp} MyPIO_STRIDE   ${MyPIO_STRIDE}
  $SUBSTITUTE ${ROMS_NLinp} MyPIO_BASE     ${MyPIO_BASE}
  $SUBSTITUTE ${ROMS_NLinp} MyPIO_REARRCOM ${MyPIO_REARRCOM}
  $SUBSTITUTE ${ROMS_NLinp} MyPIO_REARRDIR ${MyPIO_REARRDIR}
  $SUBSTITUTE ${ROMS_NLinp} MyPIO_REARR    ${MyPIO_REARR}
  $SUBSTITUTE ${ROMS_NLinp} MyFprefix      "${FprefixF}"
  $SUBSTITUTE ${ROMS_NLinp} MyFsuffix      "${Fsuffix}"
  $SUBSTITUTE ${ROMS_NLinp} MyININAME      "${ROMS_INI_F}"
  $SUBSTITUTE ${ROMS_NLinp} MyAPARNAM      "${Inp4DVAR_nl}"

  echo "   Creating DA ROMS Standard Input Script: ${ROMS_DAinp}"

  if [ -f ${ROMS_DAinp} ]; then
    /bin/rm ${ROMS_DAinp}
  fi
  cp -f ../${ROMS_DAtmp} ${ROMS_DAinp}

  $SUBSTITUTE ${ROMS_DAinp} MyNtileI       ${nPETsX}
  $SUBSTITUTE ${ROMS_DAinp} MyNtileJ       ${nPETsY}
  $SUBSTITUTE ${ROMS_DAinp} MyNouter       ${MyNouter}
  $SUBSTITUTE ${ROMS_DAinp} MyNinner       ${MyNinner}
  $SUBSTITUTE ${ROMS_DAinp} MyNHIS         ${MyNHIS_da}
  $SUBSTITUTE ${ROMS_DAinp} MyNDEFHIS      ${MyNDEFHIS_da}
  $SUBSTITUTE ${ROMS_DAinp} MyNXTR         ${MyNXTR_da}
  $SUBSTITUTE ${ROMS_DAinp} MyNDEFXTR      ${MyNDEFXTR_da}
  $SUBSTITUTE ${ROMS_DAinp} MyNQCK         ${MyNQCK_da}
  $SUBSTITUTE ${ROMS_DAinp} MyNDEFQCK      ${MyNDEFQCK_da}
  $SUBSTITUTE ${ROMS_DAinp} MyNTLM         ${MyNTLM_da}
  $SUBSTITUTE ${ROMS_DAinp} MyNADJ         ${MyNADJ_da}
  $SUBSTITUTE ${ROMS_DAinp} MyNSFF         ${MyNSFF_da}
  $SUBSTITUTE ${ROMS_DAinp} MyNOBC         ${MyNOBC_da}
  $SUBSTITUTE ${ROMS_DAinp} MyDSTART       "${DSTART}.0d0"
  $SUBSTITUTE ${ROMS_DAinp} MyINP_LIB      ${MyINP_LIB}
  $SUBSTITUTE ${ROMS_DAinp} MyOUT_LIB      ${MyOUT_LIB}
  $SUBSTITUTE ${ROMS_DAinp} MyPIO_METHOD   ${MyPIO_METHOD}
  $SUBSTITUTE ${ROMS_DAinp} MyPIO_IOTASKS  ${MyPIO_IOTASKS}
  $SUBSTITUTE ${ROMS_DAinp} MyPIO_STRIDE   ${MyPIO_STRIDE}
  $SUBSTITUTE ${ROMS_DAinp} MyPIO_BASE     ${MyPIO_BASE}
  $SUBSTITUTE ${ROMS_DAinp} MyPIO_REARRCOM ${MyPIO_REARRCOM}
  $SUBSTITUTE ${ROMS_DAinp} MyPIO_REARRDIR ${MyPIO_REARRDIR}
  $SUBSTITUTE ${ROMS_DAinp} MyPIO_REARR    ${MyPIO_REARR}
  $SUBSTITUTE ${ROMS_DAinp} MyFprefix      "${FprefixC}"
  $SUBSTITUTE ${ROMS_DAinp} MyFsuffix      "${Fsuffix}"
  $SUBSTITUTE ${ROMS_DAinp} MyININAME      "${ROMS_INI_C}"
  $SUBSTITUTE ${ROMS_DAinp} MyAPARNAM      "${Inp4DVAR_da}"

##---------------------------------------------------------------------
## Run RBL4D-Var for the current time window
##---------------------------------------------------------------------

  OuterLoop=0                        # initialize outer loop counter
  Phase4DVAR="background"            # initialize 4D-Var phase

## Set observations NetCDF filename.

  OBSnameF="${FprefixF}_roms_obs_${Fsuffix}.nc4"
  OBSnameC="${FprefixC}_roms_obs_${Fsuffix}.nc4"

## Copy outer loops nonlinear model initial conditions file.

  if [[ ${ResF} -ne ${ResC} ]]; then
    echo "   Copying NLM IC file ${ROMSiniDir}/${ROMS_INI_F}  as  ${ROMS_INI_F}"

    cp ${ROMSiniDir}/${ROMS_INI_F} ${ROMS_INI_F}
    chmod u+w ${ROMS_INI_F}                  # change protection
  fi

## Copy inner loops nonlinear model initial conditions file.

  echo "   Copying NLM IC file ${DataDir}/INI/${ROMSgeniniC}  as  ${ROMS_INI_C}"

  cp ${DataDir}/INI/${ROMSgeniniC} ${ROMS_INI_C}
  chmod u+w ${ROMS_INI_C}                    # change protection

## Get a clean copy of the observation file.  This is really important
## since this file will be modified.

  echo "   Copying OBS    file ${ObsDir}/${OBSnameF}  as  ${OBSnameF}"

  cp -p ${ObsDir}/${OBSnameF} .
  chmod u+w ${OBSnameF}                      # change protection

  if [[ ${ResF} -ne ${ResC} ]]; then
    echo "   Copying OBS    file ${ObsDir}/${OBSnameC}  as  ${OBSnameC}"

    cp -p ${ObsDir}/${OBSnameC} .
    chmod u+w ${OBSnameC}                    # change protection
  fi

## Set ROMS executable file links.

  if [ "$ROMS_EXE_A" = "$ROMS_EXE_B" ]; then
    ln -sf "../${ROMS_EXE_A}" .
  else
    ln -sf "../${ROMS_EXE_A}" .
    ln -sf "../${ROMS_EXE_B}" .
  fi

  if [ ${BATCH} -eq 1 ]; then
    EXECUTE_A="${SRUN} ${ROMS_EXE_A} ${ROMS_NLinp}"
    EXECUTE_B="${SRUN} ${ROMS_EXE_B} ${ROMS_DAinp}"
  else
    EXECUTE_A="${MPIrun} ${nPETs} ${ROMS_EXE_A} ${ROMS_NLinp}"
    EXECUTE_B="${MPIrun} ${nPETs} ${ROMS_EXE_B} ${ROMS_DAinp}"
  fi


## Run 4D-Var 'background' phase ......................................

  echo
  echo "Running 4D-Var System:  Cycle = ${Cycle}" \
                             "  Outer = ${OuterLoop}" \
                             "  Phase = ${Phase4DVAR}"

## Create ROMS 4D-Var input script 'rbl4dvar_nl.in' from template.

  My4DVarScript ${DataDir} ${SUBSTITUTE} ${OuterLoop} ${Phase4DVAR} \
                ${MyTimeIAU} ${OBSnameF} ${FprefixF} ${Fsuffix} \
                ${Inp4DVAR_nl} ${ResF}

  echo "   ${EXECUTE_A}"

  if [ ${DRYRUN} -eq 0 ]; then

    nl_log="log_outer${OuterLoop}.nl"

    if [ ${BATCH} -eq 1 ]; then
      ${SRUN} ${ROMS_EXE_A} ${ROMS_NLinp} > ${nl_log}
    else
      ${MPIrun} ${nPETs} ${ROMS_EXE_A} ${ROMS_NLinp} > ${nl_log}
    fi

    if [ $? -ne 0 ] ; then
      echo 
      echo "Error while running 4D-Var System:  Cycle = ${Cycle}" \
                                             "  Outer = ${OuterLoop}" \
                                             "  Phase = ${Phase4DVAR}"
      echo "Check ${nl_log} for details ..."
      exit 1
    fi
  fi

## Rename extracted coarse-resolution NLM trajectory NetCDF file used
## to linearize coarser TLM and ADM kernels.

  if [[ $ResF -ne $ResC ]]; then
    FWD_extract="${FprefixC}_roms_fwd_${Fsuffix}.nc"
    FWD_traject="${FprefixC}_roms_fwd_${Fsuffix}_outer${OuterLoop}.nc"

    echo
    echo "   Renaming NLM trajectory ${FWD_extract}  to  ${FWD_traject}"

    if [ ${DRYRUN} -eq 0 ]; then
      mv  ${FWD_extract} ${FWD_traject} 
    fi
  fi

## Start 4D-Var outer loops :::::::::::::::::::::::::::::::::::::::::::

  while [ $OuterLoop -lt $MyNouter ]; do
  
    OuterLoop=$(( $OuterLoop + 1 ))

## Run 4D-Var 'increment' phase .......................................

    Phase4DVAR="increment"

    echo
    echo "Running 4D-Var System:  Cycle = ${Cycle}" \
                               "  Outer = ${OuterLoop}" \
                               "  Phase = ${Phase4DVAR}"

## Create ROMS 4D-Var input script 'rbl4dvar_da.in' from template.

    My4DVarScript ${DataDir} ${SUBSTITUTE} ${OuterLoop} ${Phase4DVAR} \
                  ${MyTimeIAU} ${OBSnameC} ${FprefixC} ${Fsuffix} \
                  ${Inp4DVAR_da} ${ResC}

    echo "   ${EXECUTE_B}"

    if [ ${DRYRUN} -eq 0 ]; then

      da_log="log_outer${OuterLoop}.da"

      if [ ${BATCH} -eq 1 ]; then
        ${SRUN} ${ROMS_EXE_B} ${ROMS_DAinp} > ${da_log}
      else
        ${MPIrun} ${nPETs} ${ROMS_EXE_B} ${ROMS_DAinp} > ${da_log}
      fi

      if [ $? -ne 0 ] ; then
        echo 
        echo "Error while running 4D-Var System:  Cycle = ${Cycle}" \
                                               "  Outer = ${OuterLoop}" \
			                       "  Phase = ${Phase4DVAR}"
        echo "Check ${da_log} for details ..."
        exit 1
      fi
    fi

## Run 4D-Var 'analysis' phase ........................................

    Phase4DVAR="analysis"

    echo
    echo "Running 4D-Var System:  Cycle = ${Cycle}" \
                               "  Outer = ${OuterLoop}" \
			       "  Phase = ${Phase4DVAR}"

    My4DVarScript ${DataDir} ${SUBSTITUTE} ${OuterLoop} ${Phase4DVAR} \
                  ${MyTimeIAU} ${OBSnameF} ${FprefixF} ${Fsuffix} \
                  ${Inp4DVAR_nl} ${ResF}

    echo "   ${EXECUTE_A}"

    if [ ${DRYRUN} -eq 0 ]; then

      nl_log="log_outer${OuterLoop}.nl"

      if [ ${BATCH} -eq 1 ]; then
        ${SRUN} ${ROMS_EXE_A} ${ROMS_NLinp} > ${nl_log}
      else
        ${MPIrun} ${nPETs} ${ROMS_EXE_A} ${ROMS_NLinp} > ${nl_log}
      fi

      if [ $? -ne 0 ] ; then
        echo 
        echo "Error while running 4D-Var System:  Cycle = ${Cycle}" \
                                               "  Outer = ${OuterLoop}" \
			                       "  Phase = ${Phase4DVAR}"
        echo "Check ${nl_log} for details ..."
        exit 1
      fi
    fi

## Rename extracted coarse-resolution NLM trajectory NetCDF file used
## to linearize coarser TLM and ADM kernels.

    if [[ $ResF -ne $ResC ]]; then
      FWD_extract="${FprefixC}_roms_fwd_${Fsuffix}.nc"
      FWD_traject="${FprefixC}_roms_fwd_${Fsuffix}_outer${OuterLoop}.nc"

      echo
      echo "   Renaming NLM trajectory ${FWD_extract}  to  ${FWD_traject}"

      if [ ${DRYRUN} -eq 0 ]; then
        mv  ${FWD_extract} ${FWD_traject} 
      fi
    fi

## End of outer loops :::::::::::::::::::::::::::::::::::::::::::::::::

  done

  echo
  echo "Finished 4D-Var outer loops iterations"

##---------------------------------------------------------------------
## Advance to the next RBL4D-Var cycle, if any.
##---------------------------------------------------------------------

  echo

  if [ ${BATCH} -eq 1 ]; then
    etime=`sacct -n -X -j $SLURM_JOBID --format=Elapsed | sed 's/-/ days /'`
    ptime_sec=$(date -u -d "$ptime" +"%s")
    etime_sec=$(date -u -d "$etime" +"%s")
    time_diff=`date -u -d "0 ${etime_sec} sec - ${ptime_sec} sec" +"%H:%M:%S"`
    ptime=$etime
  else
    now=`${DATE_EXE} -u +"%s"`
    time_diff=`${DATE_EXE} -u -d "0 ${now} sec - ${ptime} sec" +"%H:%M:%S"`
    ptime=$now
  fi

  echo "Finished RBL4D-Var Cycle ${Cycle},  Elapsed time = $time_diff"

  echo
  echo "Changing to directory: ${HereDir}"

  cd ../                                          # go back to start directory

  SDAY=$(( $SDAY + $INTERVAL ))

## Notice that the IC file for the coarser inner loops is generic and
## updated in ROMS "increment phase" with the first record of the
## decimated nonlinear trajectory. It is only used to process the
## coarse grid increment in the inner loops.

  ROMSiniF="${FprefixF}_roms_dai_${Fsuffix}.nc"   # new ROMS finer IC (DAI file)
  ROMSiniDir="../${RunDir}"                       # new ROMS IC directory

  echo

## End of RBL4D-Var cycle.

done

##---------------------------------------------------------------------
## Done with computations.
##---------------------------------------------------------------------

if [ ${BATCH} -eq 1 ]; then
  total_time=`sacct -n -X -j $SLURM_JOBID --format=Elapsed`
else
  days=$(( (${ptime} - ${stime}) / 86400 ))
  hms=`${DATE_EXE} -u -d "0 ${ptime} sec - ${stime} sec" +"%H:%M:%S"`

  if (( ${days} == 0 )); then
    total_time=$hms
  else
    total_time="${days}-${hms}"
  fi
fi

echo "Finished computations, Total time = $total_time"

exit 0
