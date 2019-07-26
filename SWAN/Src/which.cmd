@echo off
rem found on http://superuser.com/questions/21067/windows-equivalent-of-whereis

@echo off
@set PATH=.;%PATH%

@rem 
@rem about:  something similar like the unix-alike-which, but with
@rem         within pure cmd
@rem 

if "%1" == "" (
    @echo Usage: 
    @echo.
    @echo   which 'cmd'
    @echo.
    @echo.if 'cmd' is not found, ERRORLEVEL is set to 1
    @echo.
) else (
    ( @for %%f in (%1 %1.exe %1.cmd %1.bat %1.pif) do if not "%%~$PATH:f" == "" ( @echo %%~$PATH:f ) else @set ERRORLEVEL=1)
)
