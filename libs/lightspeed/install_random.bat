@if (%1) == () goto usage
@set VSINSTALLDIR=%1
@set VCVARSOPTS=%2
@set use_sdk=no
@if "%VCVARSOPTS%"=="x86_amd64" (
if exist "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.cmd" (
  @set use_sdk=yes
))
@if %use_sdk%==yes (
call "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.cmd" /x64
) else (
if exist %VSINSTALLDIR%\VC\vcvarsall.bat (
call %VSINSTALLDIR%\VC\vcvarsall.bat %VCVARSOPTS%
) else (
if exist %VSINSTALLDIR%\Common7\Tools\vsvars32.bat (
call %VSINSTALLDIR%\Common7\Tools\vsvars32.bat
) else (
if exist %VSINSTALLDIR%\VC98\bin\vcvars32.bat (
call %VSINSTALLDIR%\VC98\bin\vcvars32.bat
) else (
@echo "Could not find either of these files:"
@echo %VSINSTALLDIR%\VC\vcvarsall.bat
@echo %VSINSTALLDIR%\Common7\Tools\vsvars32.bat
@echo %VSINSTALLDIR%\VC98\bin\vcvars32.bat
goto :eof
))))
@rem /MD : link with MSVCRT.lib
cl  /c /MD /O2 /DNDEBUG random.c
link random.obj /dll /def:random.def
goto :eof

:usage
@echo Error in script usage.  The correct usage is:
@echo    %0 (visual studio root folder) [vcvarsall.bat options]
