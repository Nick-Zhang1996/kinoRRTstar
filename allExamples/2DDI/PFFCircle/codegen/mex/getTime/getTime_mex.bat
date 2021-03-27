@echo off
set MATLAB=E:\MATLAB\R2018b
set MATLAB_ARCH=win64
set MATLAB_BIN="E:\MATLAB\R2018b\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=getTime_mex
set MEX_NAME=getTime_mex
set MEX_EXT=.mexw64
call setEnv.bat
echo # Make settings for getTime > getTime_mex.mki
echo CC=%COMPILER%>> getTime_mex.mki
echo CXX=%CXXCOMPILER%>> getTime_mex.mki
echo CFLAGS=%COMPFLAGS%>> getTime_mex.mki
echo CXXFLAGS=%CXXCOMPFLAGS%>> getTime_mex.mki
echo LINKER=%LINKER%>> getTime_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> getTime_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> getTime_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> getTime_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> getTime_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> getTime_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> getTime_mex.mki
echo OMPFLAGS= >> getTime_mex.mki
echo OMPLINKFLAGS= >> getTime_mex.mki
echo EMC_COMPILER=mingw64>> getTime_mex.mki
echo EMC_CONFIG=optim>> getTime_mex.mki
"E:\MATLAB\R2018b\bin\win64\gmake" -j 1 -B -f getTime_mex.mk
