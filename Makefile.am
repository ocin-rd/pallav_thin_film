
# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS= axisym_thinfilm_dripping_faucet

#----------------------------------------------------------------------

# Sources for executable
axisym_thinfilm_dripping_faucet_SOURCES = axisym_thinfilm_dripping_faucet.cc
# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
axisym_thinfilm_dripping_faucet_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
INCLUDES = -I@includedir@  

#EXTRA_DIST=axisym_thinfilm_dripping_faucet_elements.h axisym_thinfilm_dripping_faucet_elements.cc


