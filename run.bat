@echo off
set MPI_BIN=C:\Program Files\Microsoft MPI\Bin
set EXECUTABLE=a2.exe  REM
mpiexec -n 10 "%MPI_BIN%\mpiexec.exe" %EXECUTABLE%
pause