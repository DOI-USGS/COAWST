@echo off
rem
rem    swanrun.bat
rem
rem    Run the SWAN program by means of the SWAN input file
rem    Note: it is assumed that the extension of the input file is '.swn'
rem
rem    Usage: swanrun inputfile [nprocs]
rem

set nprocs=1

if not "%1"=="" goto OK1
  echo.
  echo Usage: swanrun inputfile [nprocs]
  goto END
:OK1

set inputfile=%1
shift

if exist %inputfile%.swn goto OK2
  echo.
  echo Error: file %inputfile%.swn does not exist
  goto END
:OK2

if "%1"=="" goto OK3
  set nprocs=%1
:OK3

copy %inputfile%.swn INPUT >> nul

if not %nprocs%==1 goto PARALLEL1
swan.exe
goto OK4

:PARALLEL1
mpiexec -n %nprocs% swan.exe

:OK4
if errorlevel 1 goto END

if not %nprocs%==1 goto PARALLEL2
  if exist PRINT copy PRINT %inputfile%.prt >> nul
  if exist PRINT del PRINT
  goto OK5
:PARALLEL2
if not exist PRINT-001 goto OK5
if %nprocs% GTR 9 goto RANGE1
  for /L %%i in (1,1,%nprocs%) do copy PRINT-00%%i %inputfile%.prt-00%%i >> nul
  for /L %%i in (1,1,%nprocs%) do del PRINT-00%%i
  goto OK5
:RANGE1
if %nprocs% GTR 99 goto RANGE2
  for /L %%i in (1,1,9) do copy PRINT-00%%i %inputfile%.prt-00%%i >> nul
  for /L %%i in (1,1,9) do del PRINT-00%%i
  for /L %%i in (10,1,%nprocs%) do copy PRINT-0%%i %inputfile%.prt-0%%i >> nul
  for /L %%i in (10,1,%nprocs%) do del PRINT-0%%i
  goto OK5
:RANGE2
if %nprocs% GTR 999 goto ERR
  for /L %%i in (1,1,9) do copy PRINT-00%%i %inputfile%.prt-00%%i >> nul
  for /L %%i in (1,1,9) do del PRINT-00%%i
  for /L %%i in (10,1,99) do copy PRINT-0%%i %inputfile%.prt-0%%i >> nul
  for /L %%i in (10,1,99) do del PRINT-0%%i
  for /L %%i in (100,1,%nprocs%) do copy PRINT-%%i %inputfile%.prt-%%i >> nul
  for /L %%i in (100,1,%nprocs%) do del PRINT-%%i
  goto OK5
:ERR
  echo Error: too many processes
  goto END
:OK5

if exist Errfile copy Errfile %inputfile%.erf >> nul
if exist Errfile del Errfile
if exist ERRPTS copy ERRPTS %inputfile%.erp >> nul
if exist ERRPTS del ERRPTS
del INPUT
if not exist norm_end goto END
  type norm_end

:END

set inputfile=
set nprocs=
