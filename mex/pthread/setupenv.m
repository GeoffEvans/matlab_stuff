% prepare the environment (we could also use the -I, -L flags)
pthreadsInstallFolder = 'C:\Program Files (x86)\pthreads-win32\';  % change this as needed
setenv('PATH',   [getenv('PATH')    ';' pthreadsInstallFolder 'dll\x64']);
setenv('LIB',    [getenv('LIB')     ';' pthreadsInstallFolder 'lib\x64']);
setenv('INCLUDE',[getenv('INCLUDE') ';' pthreadsInstallFolder 'include']);
 
% create a 64-bit MEX that uses the pthreads DLL
mex expt_pt.c -lpthreadVC2
 
% copy the pthreadVC2.dll file to be accessible to the MEX file, otherwise it will not run
copyfile([pthreadsInstallFolder 'dll\x64\pthreadVC2.dll'], '.')
