#!/bin/bash
#
# git $Id$
# svn $Id: submit_split_i4dvar.sh 1151 2023-02-09 03:08:53Z arango $
#######################################################################
## Copyright (c) 2002-2023 The ROMS/TOMS Group                        #
##   Licensed under a MIT/X style license                             #
##   See License_ROMS.txt                                             #
################################################## Hernan G. Arango ###
##                                                                    #
## ROMS Split I4D-Var Data Assimilation running BASH script:          #
##                                                                    #
## The I4D-Var is split in several phases and executables:            #
##                                                                    #
##   1) background                          ROMS_EXE_A                #
##   2) increment                           ROMS_EXE_B                #
##   3) analysis                            ROMS_EXE_A                #
##   4) posterior_analysis                  ROMS_EXE_A                #
##                                                                    #
## The 4D-Var algorithm is split into Executable A and B:             #
##                                                                    #
## ROMS_EXE_A: Computes ROMS NLM "background" trajectory used to      #
##             linearize the TLM and ADM kernels used in the 4D-Var   #
##             minimization. It also computes the "posterior_analysis"#
##             solution. It interpolates the NLM solution to the      #
##             observation locations in space and time. The NLM could #
##             be a part of a coupling system and or include nested   #
##             grids.                                                 #
##                                                                    #
##             It is also used in the "analysis" phase to generate    #
##             the NLM initial conditions by adding the 4D-Var        #
##             increments to the background.                          #
##                                                                    #
## ROMS_EXE_B: It is used in the "increment" phase. The 4D-Var        #
##             increment is obtained by minimizing the cost function  #
##             over Ninner loops. It is possible to use a coarser     #
##             grid resolution in the iterations of the inner loops.  #
##             If so, the "background" trajectory needs to be         #
##             interpolated into the coarser grid. Then, after the    #
##             inner loops, the coarser grid increment needs to be    #
##             interpolated to the finer grid before the "analysis"   #
##             phase. The "increment" phase could be run at a lower   #
##             precision.                                             #
##                                                                    #
## Notice that this script has a customizable function, My4DVarScript,#
## at the top that creates "i4dvar.in" from a template, which is      #
## overwritten in each phase. The "i4dvar.in" passes to ROMS the      #
## phase and outer loop value to compute.                             #
##                                                                    #
## Then, there is the customizable section for the computer batch     #
## directives and the tunable parameters.                             #
##                                                                    #
## I4D-Var phases workflow:                                           #
##                                                                    #
##   CALL prior_error                                                 #
##                                                                    #
##   OuterLoop : DO outer=1:Nouter                                    #
##     CALL background (outer, RunInterval)                           #
##     CALL increment (outer, RunInterval)                Inner-loops #
##     CALL analysis (outer, RunInterval)                             #
##   END DO OuterLoop                                                 #
##                                                                    #
##   CALL posterior_analysis (RunInterval)                            #
##                                                                    #
#######################################################################

##  SLURM configuration for amarel:
##
##  Use 'sbatch submit_split_i4dvar.sh'    to queue the job
##  Use 'squeue -p p_omg_1'                to check our group jobs (including JOBID)
##  Use 'sacct -j JOBID -l'                to check job accounting data
##  Use 'scontrol show job JOBID'          to check job configuration
##  Use 'sinfo JOBID'                      to check job information
##  Use 'scancel JOBID'                    to cancel a job

#SBATCH --exclusive                      # don't run on nodes with other jobs running
#SBATCH --partition=p_omg_1              # Partition (job queue), NodeList: 108,116-120
#SBATCH --requeue                        # Return job to the queue if preempted
#SBATCH --job-name=WC13_split_i4dvar     # Assign an short name to your job
#SBATCH --nodes=1                        # Number of nodes you require (each has 32 PETs)
#SBATCH --ntasks=12                      # Total number of tasks you'll launch
#SBATCH --ntasks-per-node=32             # Number of tasks you'll launch on each node
#SBATCH --cpus-per-task=1                # Cores per task (>1 if multithread tasks)
#SBATCH --mem=177000                     # Real memory (RAM) required (MB)
#SBATCH --time=00-01:00:00               # Total run time limit (DD-HH:MM:SS)
#SBATCH --output=log.%N.%j               # STDOUT output file
#SBATCH --error=err.%N.%j                # STDERR output file (optional)
#SBATCH --export=ALL                     # Export you current env to the job env

#######################################################################
## I4D-Var data assimilation input script function.  It generates     #
## 'i4dvar.in' from the 's4dvar.in' template.                         #
#######################################################################

## Start of My4DVarScript() function definition

My4DVarScript() {

     DataDir=$1                # Data directory
  SUBSTITUTE=$2                # ROMS Perl subtitution function
   OuterLoop=$3                # current outer loop counter
  Phase4DVAR=$4                # current 4D-Var computation phase
     OBSname=$5                # 4D-Var observations NetCDF file
     Fprefix=$6                # ROMS output files prefix (use roms_app)
     Fsuffix=$7                # ROMS output files suffix
    Inp4DVAR=$8                # 4D-Var standard input

  echo
  if [ "${Phase4DVAR}" = "post_analysis" ]; then
    echo "   Creating 4D-Var Input Script from Template: ${Inp4DVAR}" \
          "  Phase = ${Phase4DVAR}"
  else
    echo "   Creating 4D-Var Input Script from Template: ${Inp4DVAR}" \
          "  Outer = ${OuterLoop}  Phase = ${Phase4DVAR}"
  fi
  echo

## Set model, initial conditions, boundary conditions and surface
## forcing error covariance standard deviations files.

 STDnameM=${DataDir}/wc13_std_m.nc
 STDnameI=${DataDir}/wc13_std_i.nc
 STDnameB=${DataDir}/wc13_std_b.nc
 STDnameF=${DataDir}/wc13_std_f.nc

## Set model, initial conditions, boundary conditions and surface
## forcing error covariance normalization factors files.

 NRMnameM=${DataDir}/wc13_nrm_m.nc
 NRMnameI=${DataDir}/wc13_nrm_i.nc
 NRMnameB=${DataDir}/wc13_nrm_b.nc
 NRMnameF=${DataDir}/wc13_nrm_f.nc

## Modify 4D-Var template input script and specify above files.

 if [ -f $Inp4DVAR ]; then
   /bin/rm ${Inp4DVAR}
 fi

 cp ../s4dvar.in ${Inp4DVAR}

 $SUBSTITUTE $Inp4DVAR MyOuterLoop   ${OuterLoop}
 $SUBSTITUTE $Inp4DVAR MyPhase4DVAR  ${Phase4DVAR}
 $SUBSTITUTE $Inp4DVAR roms_std_m.nc ${STDnameM}
 $SUBSTITUTE $Inp4DVAR roms_std_i.nc ${STDnameI}
 $SUBSTITUTE $Inp4DVAR roms_std_b.nc ${STDnameB}
 $SUBSTITUTE $Inp4DVAR roms_std_f.nc ${STDnameF}
 $SUBSTITUTE $Inp4DVAR roms_nrm_m.nc ${NRMnameM}
 $SUBSTITUTE $Inp4DVAR roms_nrm_i.nc ${NRMnameI}
 $SUBSTITUTE $Inp4DVAR roms_nrm_b.nc ${NRMnameB}
 $SUBSTITUTE $Inp4DVAR roms_nrm_f.nc ${NRMnameF}
 $SUBSTITUTE $Inp4DVAR roms_obs.nc   ${OBSname}
 $SUBSTITUTE $Inp4DVAR roms_hss.nc   ${Fprefix}_hss_${Fsuffix}.nc
 $SUBSTITUTE $Inp4DVAR roms_lcz.nc   ${Fprefix}_lcz_${Fsuffix}.nc
 $SUBSTITUTE $Inp4DVAR roms_lze.nc   ${Fprefix}_lze_${Fsuffix}.nc
 $SUBSTITUTE $Inp4DVAR roms_mod.nc   ${Fprefix}_mod_${Fsuffix}.nc
 $SUBSTITUTE $Inp4DVAR roms_err.nc   ${Fprefix}_err_${Fsuffix}.nc
}

## End of My4DVarScript() function definition

##---------------------------------------------------------------------
## Control switches: What do you want to do?
##---------------------------------------------------------------------

#      DRYRUN=1                # Print configuration but do not execute
       DRYRUN=0                # Run 4D-Var cycle

        BATCH=0                # No batch system submission
#       BATCH=1                # Use batch system SLURM to submit

##---------------------------------------------------------------------
## User tunable parameters. If you follow recommendations, this is
## the only section that you need to customize..
##---------------------------------------------------------------------

     ROMS_APP="WC13"                   # ROMS Application CPP

     roms_app=`echo ${ROMS_APP} | tr '[:upper:]' '[:lower:]'` # lowercase

    ROMS_ROOT=${HOME}/ocean/repository/4dvar

      HereDir=${PWD}                   # current directory

      DataDir="../../Data"             # data directory

       ObsDir=${DataDir}               # observations directory

 if [[ "$OSTYPE" == "darwin"* ]]; then
     DATE_EXE=gdate                    # macOS system GNU date
 else
     DATE_EXE=date                     # Linux system date
 fi

 if [ ${BATCH} -eq 1 ]; then
         SRUN="srun --mpi=pmi2"        # SLURM workload manager
 else
       MPIrun="mpirunI -np"            # Basic MPI workload manager
 fi

   ROMS_EXE_A="romsM"                  # ROMS executable A
   ROMS_EXE_B="romsM"                  # ROMS executable B

   START_DATE="2004-01-03"             # 4D-Var starting date

 FRST_INI_DAY="${START_DATE}"          # first cycle initialization date
#FRST_INI_DAY="2004-01-07"             # restart initialization date

 LAST_INI_DAY="${START_DATE}"          # last cycle initialzation date
#LAST_INI_DAY="2004-01-11"             # last cycle initialzation date

 ROMS_TIMEREF="1968-05-23"             # ROMS time reference date

     INTERVAL=4                        # 4D-Var interval window (days)

       nPETsX=3                        # number PETs in the X-direction
       nPETsY=4                        # number PETs in the Y-direction

     MyNouter=1                        # number of 4D-Var outer loops
#    MyNouter=2                        # number of 4D-Var outer loops
     MyNinner=25                       # number of 4D-Var inner loops

       MyNHIS=4                        # NLM trajectory is saved every 2 hours
    MyNDEFHIS=0                        # No multi-file NLM trajectory
#      MyNHIS=48                       # NLM trajectory is saved every 24 hours

       MyNQCK=4                        # NLM quicksave trajectory is saved every 2 hours

      restart=0                        # restart 4D-Var windows (0:no, 1:yes)
#     restart=1                        # restart 4D-Var windows (0:no, 1:yes)

   ROMS_NLpre="roms_nl_${roms_app}"    # ROMS NLM stdinp prefix
   ROMS_NLtmp="${ROMS_NLpre}.tmp"      # ROMS NLM stdinp template

   ROMS_DApre="roms_da_${roms_app}"    # ROMS ADM/TLM stdinp prefix
   ROMS_DAtmp="${ROMS_DApre}.tmp"      # ROMS ADM/TLM stdinp template

 if [ ${restart} -eq 0 ]; then
      ROMSini="wc13_roms_ini_20040103.nc"    # ROMS IC
   ROMSiniDir=${DataDir}                     # ROMS IC directory
 else
      ROMSini="wc13_roms_dai_20040103.nc"    # restart ROMS IC
   ROMSiniDir=../2004.01.03                  # restart ROMS IC directory
 fi

     Inp4DVAR="i4dvar.in"                    # ROMS 4D-Var input script

 if [ ${BATCH} -eq 1 ]; then
        etime=00:00:00                       # elapsed time
        ptime=${etime}                       # previous time
 else
        stime=`${DATE_EXE} -u +"%s"`         # start time (sec) since epoch
        ptime=$stime                         # previous time
 fi

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
echo " ROMS Split I4D-Var Data Assimilation: ${ROMS_APP}"
echo "${separator1}"
echo
echo "                    ROMS Root: ${ROMS_ROOT}"
echo "            ROMS Executable A: ${ROMS_EXE_A}  (background," \
                                     "analysis, post_analysis)"
echo "            ROMS Executable B: ${ROMS_EXE_B}  (increment)"
echo

##---------------------------------------------------------------------
## Compute date number for reference date, and first and last
## initialization dates.
##---------------------------------------------------------------------

         S_DN=`${ROMS_ROOT}/ROMS/Bin/dates datenum ${START_DATE}`
         F_DN=`${ROMS_ROOT}/ROMS/Bin/dates datenum ${FRST_INI_DAY}`
         L_DN=`${ROMS_ROOT}/ROMS/Bin/dates datenum ${LAST_INI_DAY}`
       REF_DN=`${ROMS_ROOT}/ROMS/Bin/dates datenum ${ROMS_TIMEREF}`

echo "         4D-Var Starting Date: ${START_DATE}  datenum = ${S_DN}"
echo "      First 4D-Var Cycle Date: ${FRST_INI_DAY}  datenum = ${F_DN}"
echo "      Last  4D-Var Cycle Date: ${LAST_INI_DAY}  datenum = ${L_DN}"
echo "          ROMS Reference Date: ${ROMS_TIMEREF}  datenum = ${REF_DN}"
echo "          4D-Var Cycle Window: ${INTERVAL} days"
echo "      Number of parallel PETs: ${nPETs}  (${nPETsX}x${nPETsY})"
echo "    Current Stating Directory: ${HereDir}"
echo "         ROMS Application CPP: ${ROMS_APP}"
echo "      Descriptor in filenames: ${roms_app}"
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

## Set coupling parameters.

          Cycle=$(( $Cycle + 1 ))            # advance cycle by one
         DSTART=$(( $SDAY - $REF_DN ))       # ROMS DSTART parameter
           EDAY=$(( $SDAY + $INTERVAL ))     # end day for 4D-Var cycle

          RDATE=`${ROMS_ROOT}/ROMS/Bin/dates numdate ${REF_DN}`
          SDATE=`${ROMS_ROOT}/ROMS/Bin/dates numdate ${SDAY}`
          EDATE=`${ROMS_ROOT}/ROMS/Bin/dates numdate ${EDAY}`
            DOY=`${ROMS_ROOT}/ROMS/Bin/dates yday ${SDATE}`
           yday=`printf %03d $DOY`

  ReferenceTime=`${DATE_EXE} -d "${RDATE}" '+%Y %m %d %H %M %S'`
      StartTime=`${DATE_EXE} -d "${SDATE}" '+%Y %m %d %H %M %S'`
       StopTime=`${DATE_EXE} -d "${EDATE}" '+%Y %m %d %H %M %S'`
    RestartTime="${StartTime}"
        Fprefix="${roms_app}"
        Fsuffix=`${DATE_EXE} -d "${SDATE}" '+%Y%m%d'`
         RunDir=`${DATE_EXE} -d "${SDATE}" '+%Y.%m.%d'`
       ROMS_INI="${ROMSini}"

  echo "${separator2}"
  echo
  echo "            4D-Var Cycle Date: ${SDATE}  DayOfYear = ${yday}" \
                                       " Cycle = ${Cycle}"
  echo "            Run sub-directory: ${RunDir}"
  echo "        Number of outer loops: ${MyNouter}"
  echo "        Number of inner loops: ${MyNinner}"
  echo "       NLM trajectory writing: ${MyNHIS} NHIS timesteps"
  echo "       NLM quicksave  writing: ${MyNQCK} NQCK timesteps"
  echo "    NLM multi-file trajectory: ${MyNDEFHIS} NDEFHIS timesteps"
  echo "                  ROMS DSTART: ${DSTART}.0d0"
  echo "             I/O Files Prefix: ${Fprefix}"
  echo "             I/O Files Suffix: ${Fsuffix}"
  echo "                ReferenceTime: ${ReferenceTime}"
  echo "             4D-Var StartTime: ${StartTime}"
  echo "              4D-Var StopTime: ${StopTime}"
  echo "      ROMS Initial Conditions: ${ROMSiniDir}/${ROMS_INI}"

##---------------------------------------------------------------------
## Create run sub-directory based on cycle date (YYYY.MM.DD) and create
## ROMS standard input script from template.
##---------------------------------------------------------------------

  ROMS_NLinp=`echo ${ROMS_NLpre}_${Fsuffix}'.in'`
  ROMS_DAinp=`echo ${ROMS_DApre}_${Fsuffix}'.in'`

  echo "NL ROMS Standard Input Script: ${ROMS_NLinp}"
  echo "DA ROMS Standard Input Script: ${ROMS_DAinp}"
  echo "          4D-Var Input Script: ${Inp4DVAR}"
  echo

  if [ ${DRYRUN} -eq 1 ]; then                # if dry-run, remove run
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

  echo "   Creating NL ROMS Standart Input Script: ${ROMS_NLinp}"

  if [ -f ${ROMS_NLinp} ]; then
    /bin/rm ${ROMS_NLinp}
  fi
  cp -f ../${ROMS_NLtmp} ${ROMS_NLinp}

  $SUBSTITUTE ${ROMS_NLinp} MyNtileI  ${nPETsX}
  $SUBSTITUTE ${ROMS_NLinp} MyNtileJ  ${nPETsY}
  $SUBSTITUTE ${ROMS_NLinp} MyNouter  ${MyNouter}
  $SUBSTITUTE ${ROMS_NLinp} MyNinner  ${MyNinner}
  $SUBSTITUTE ${ROMS_NLinp} MyNHIS    ${MyNHIS}
  $SUBSTITUTE ${ROMS_NLinp} MyNDEFHIS ${MyNDEFHIS}
  $SUBSTITUTE ${ROMS_NLinp} MyNQCK    ${MyNQCK}
  $SUBSTITUTE ${ROMS_NLinp} MyDSTART  "${DSTART}.0d0"
  $SUBSTITUTE ${ROMS_NLinp} MyFprefix "${Fprefix}"
  $SUBSTITUTE ${ROMS_NLinp} MyFsuffix "${Fsuffix}"
  $SUBSTITUTE ${ROMS_NLinp} MyININAME "${ROMS_INI}"
  $SUBSTITUTE ${ROMS_NLinp} MyAPARNAM "${Inp4DVAR}"

  echo "   Creating DA ROMS Standart Input Script: ${ROMS_DAinp}"

  if [ -f ${ROMS_DAinp} ]; then
    /bin/rm ${ROMS_DAinp}
  fi
  cp -f ../${ROMS_DAtmp} ${ROMS_DAinp}

  $SUBSTITUTE ${ROMS_DAinp} MyNtileI  ${nPETsX}
  $SUBSTITUTE ${ROMS_DAinp} MyNtileJ  ${nPETsY}
  $SUBSTITUTE ${ROMS_DAinp} MyNouter  ${MyNouter}
  $SUBSTITUTE ${ROMS_DAinp} MyNinner  ${MyNinner}
  $SUBSTITUTE ${ROMS_DAinp} MyNHIS    ${MyNHIS}
  $SUBSTITUTE ${ROMS_DAinp} MyNDEFHIS ${MyNDEFHIS}
  $SUBSTITUTE ${ROMS_DAinp} MyNQCK    ${MyNQCK}
  $SUBSTITUTE ${ROMS_DAinp} MyDSTART  "${DSTART}.0d0"
  $SUBSTITUTE ${ROMS_DAinp} MyFprefix "${Fprefix}"
  $SUBSTITUTE ${ROMS_DAinp} MyFsuffix "${Fsuffix}"
  $SUBSTITUTE ${ROMS_DAinp} MyININAME "${ROMS_INI}"
  $SUBSTITUTE ${ROMS_DAinp} MyAPARNAM "${Inp4DVAR}"

##---------------------------------------------------------------------
## Run 4D-Var for the current time window
##---------------------------------------------------------------------

  OuterLoop=0                        # initialize outer loop counter

## Set observations NetCDF filename.

  OBSname="${Fprefix}_obs_${Fsuffix}.nc"

## Copy nonlinear model initial conditions file.

  echo "   Copying NLM IC file ${ROMSiniDir}/${ROMS_INI}  as  ${ROMS_INI}"

  cp ${ROMSiniDir}/${ROMS_INI} ${ROMS_INI}
  chmod u+w ${ROMS_INI}                      # change protection

## Get a clean copy of the observation file.  This is really important
## since this file will be modified.

  echo "   Copying OBS file ${ObsDir}/${OBSname}  as  ${OBSname}"

  cp -p ${ObsDir}/${OBSname} .
  chmod u+w ${OBSname}                       # change protection

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

## Start 4D-Var outer loops :::::::::::::::::::::::::::::::::::::::::::

  while [ $OuterLoop -lt $MyNouter ]; do

    OuterLoop=$(( $OuterLoop + 1 ))

## Run 4D-Var 'background' phase ......................................

    Phase4DVAR="background"            # background 4D-Var phase

    echo
    echo "Running 4D-Var System:  Cycle = ${Cycle}" \
                               "  Outer = ${OuterLoop}" \
                               "  Phase = ${Phase4DVAR}"

## Create ROMS 4D-Var input script 'i4dvar.in' from template.

    My4DVarScript ${DataDir} ${SUBSTITUTE} ${OuterLoop} ${Phase4DVAR} \
                  ${OBSname} ${Fprefix} ${Fsuffix} ${Inp4DVAR}

    echo "   ${EXECUTE_A}"

    if [ ${DRYRUN} -eq 0 ]; then

      if [ ${BATCH} -eq 1 ]; then
        ${SRUN} ${ROMS_EXE_A} ${ROMS_NLinp}
      else
        ${MPIrun} ${nPETs} ${ROMS_EXE_A} ${ROMS_NLinp} >> err
      fi

      if [ $? -ne 0 ] ; then
        echo
        echo "Error while running 4D-Var System:  Cycle = ${Cycle}" \
                                               "  Outer = ${OuterLoop}" \
                                               "  Phase = ${Phase4DVAR}"
        echo "Check ${RunDir}/log.roms for details ..."

        exit 1
      fi
    fi

## Run 4D-Var 'increment' phase.

    Phase4DVAR="increment"

    echo
    echo "Running 4D-Var System:  Cycle = ${Cycle}" \
                               "  Outer = ${OuterLoop}" \
                               "  Phase = ${Phase4DVAR}"

    My4DVarScript ${DataDir} ${SUBSTITUTE} ${OuterLoop} ${Phase4DVAR} \
                  ${OBSname} ${Fprefix} ${Fsuffix} ${Inp4DVAR}

    echo "   ${EXECUTE_B}"

    if [ ${DRYRUN} -eq 0 ]; then

      if [ ${BATCH} -eq 1 ]; then
        ${SRUN} ${ROMS_EXE_B} ${ROMS_DAinp}
      else
        ${MPIrun} ${nPETs} ${ROMS_EXE_B} ${ROMS_DAinp} >> err
      fi

      if [ $? -ne 0 ] ; then
        echo
        echo "Error while running 4D-Var System:  Cycle = ${Cycle}" \
                                               "  Outer = ${OuterLoop}" \
			                       "  Phase = ${Phase4DVAR}"
        echo "Check ${RunDir}/log.roms for details ..."
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
                  ${OBSname} ${Fprefix} ${Fsuffix} ${Inp4DVAR}

    echo "   ${EXECUTE_A}"

    if [ ${DRYRUN} -eq 0 ]; then

      if [ ${BATCH} -eq 1 ]; then
        ${SRUN} ${ROMS_EXE_A} ${ROMS_NLinp}
      else
        ${MPIrun} ${nPETs} ${ROMS_EXE_A} ${ROMS_NLinp} >> err
      fi

      if [ $? -ne 0 ] ; then
        echo
        echo "Error while running 4D-Var System:  Cycle = ${Cycle}" \
                                               "  Outer = ${OuterLoop}" \
			                       "  Phase = ${Phase4DVAR}"
        echo "Check ${RunDir}/log.roms for details ..."
        exit 1
      fi
    fi

## End of outer loops :::::::::::::::::::::::::::::::::::::::::::::::::

  done

  echo
  echo "Finished 4D-Var outer loops iterations"

## Run nonlinear model to compute posterior analysis ..................
##

  Phase4DVAR="post_analysis"

  echo
  echo "Running 4D-Var System:  Cycle = ${Cycle}" \
                             "  Phase = ${Phase4DVAR}"

  My4DVarScript ${DataDir} ${SUBSTITUTE} ${OuterLoop} ${Phase4DVAR} \
                ${OBSname} ${Fprefix} ${Fsuffix} ${Inp4DVAR}

  echo "   ${EXECUTE_A}"

  if [ ${DRYRUN} -eq 0 ]; then

    if [ ${BATCH} -eq 1 ]; then
      ${SRUN} ${ROMS_EXE_A} ${ROMS_NLinp}
    else
      ${MPIrun} ${nPETs} ${ROMS_EXE_A} ${ROMS_NLinp} >> err
    fi

    if [ $? -ne 0 ] ; then
      echo
      echo "Error while running 4D-Var System:  Cycle = ${Cycle}" \
		                             "  Phase = ${Phase4DVAR}"
      echo "Check ${RunDir}/log.roms for details ..."
      exit 1
    fi
  fi

##---------------------------------------------------------------------
## Advance to the 4D-Var cycle, if any.
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

  echo "Finished I4D-Var Cycle ${Cycle},  Elapsed time = $time_diff"

  echo
  echo "Changing to directory: ${HereDir}"

  cd ../                                          # go back to start directory

  SDAY=$(( $SDAY + $INTERVAL ))

  ROMSini="${Fprefix}_roms_dai_${Fsuffix}.nc"     # new ROMS IC (DAI file)
  ROMSiniDir="../${RunDir}"                       # new ROMS IC directory

  echo

## End of 4D-Var cycle.

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
