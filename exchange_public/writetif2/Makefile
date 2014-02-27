# mule-utf-8
# expects tiff-win32-3.6.1-2.exe installation in c:\GnuWin32
# BE SURE TO ADD C:\GnuWin32\bin TO PATH!!
# works with MATLAB R2008a
# build with GHC 6.8.3


MATLAB=C:\matlab\r2008a

CC = ghc
.SUFFIXES : .c .hs .o

MEXFLAGS = -I$(MATLAB)\extern\include

bindir = bin\win32
objdir = $(bindir)\objs

# FFI stub files go to objdir 
CFLAGS = $(MEXFLAGS) -I$(objdir)
LFLAGS = $(MATLAB)\extern\lib\win32\microsoft\libmx.lib $(MATLAB)\extern\lib\win32\microsoft\libmex.lib $(MATLAB)\extern\lib\win32\microsoft\libmat.lib \
	       C:\GnuWin32\lib\libtiff.dll.a


# all: writetifcn.mexw32
all: $(objdir) $(bindir)\writetif2.mexw32

# create directory if missing
$(objdir) :
	mkdir $(objdir)

# Actually, Mex.o goes to $(objdir)/Matlab
$(objdir)\Wtifc.o: $(objdir)\Mex.o $(objdir)\LibTIFF.o

Wtifc_stub.h: $(objdir)\Wtifc.o

$(objdir)\dllMain.o: Wtifc_stub.h

$(bindir)\writetif2.mexw32: $(objdir)\dllMain.o
	$(CC) -shared -o "$(@)" -odir $(objdir) $(objdir)\*.o bin\win32\objs\Matlab\*.o bin\win32\objs\Tiff\*.o writetif2.def $(LFLAGS)

{src}.c{$(objdir)}.o:
	@echo CC  $(<) ===^> $(@)
	@$(CC) -c $(CFLAGS) "$(<)" -o "$(@)"

.hs{$(objdir)}.o:
	@echo GHC $(<) ===^> $(@)
	@$(CC) -c "$(<)" -fglasgow-exts -odir $(objdir) -stubdir $(objdir)

# VPATH??
# Actually, Mex.o goes to $(objdir)/Matlab
{Matlab}.hs{$(objdir)}.o:
	@echo GHC $(<) ===^> $(@)
	@$(CC) -c "$(<)" -fglasgow-exts -odir $(objdir) -stubdir $(objdir)
{Tiff}.hs{$(objdir)}.o:
	@echo GHC $(<) ===^> $(@)
	@$(CC) -c "$(<)" -fglasgow-exts -odir $(objdir) -stubdir $(objdir)


clean:
	@-del /Q /S $(objdir)
	@-del $(bindir)\writetif2.mexw32
	@-del $(bindir)\writetif2.mexw32.a


# Local Variables: 
# compile-command:"\"c:/Program Files/Microsoft Visual Studio 9.0/VC/BIN/nmake\""
# End:
