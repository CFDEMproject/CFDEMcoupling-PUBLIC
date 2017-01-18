# Specify additional include and library paths, as well as libraries for the compilation
#
# CFDEM_ADD_INC =
# CFDEM_ADD_LIB_PATHS =
# CFDEM_ADD_LIBS =

# additional static libraries to be linked to lagrangian library
CFDEM_ADD_STATICLIBS = \
-lmpi_cxx \

# If you don't want VTK comment the following line and use the appropriate LIGGGHTS Makefile
# via setting CFDEM_LIGGGHTS_MAKEFILE_NAME that does not contain VTK.
include $(CFDEM_ADD_LIBS_DIR)/additionalLibs_vtk
include $(CFDEM_ADD_LIBS_DIR)/additionalLibs_superquadric

#################################################################
## SETTINGS FOR 2.4.x                                          ##
#################################################################
#----------------------------------------------------------------
# incompressible turbulence model settings
#----------------------------------------------------------------
# paths for incompressible turbulence models to use
CFDEM_ADD_INCOMPTURBMOD_PATHS = \
-I$(LIB_SRC)/turbulenceModels/incompressible/turbulenceModel \
-I$(LIB_SRC)/fvOptions/lnInclude \
-I$(LIB_SRC)/sampling/lnInclude

# libs for turbulence models to use
CFDEM_ADD_INCOMPTURBMOD_LIBS = \
-lincompressibleRASModels \
-lincompressibleLESModels \
-lfvOptions \
-lsampling


#----------------------------------------------------------------
# compressible turbulence model settings
#----------------------------------------------------------------
# paths for compressible turbulence models to use
CFDEM_ADD_COMPTURBMOD_PATHS = \
-I$(LIB_SRC)/turbulenceModels/compressible/turbulenceModel \
-I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
-I$(LIB_SRC)/thermophysicalModels/radiationModels/lnInclude \

# libs for turbulence models to use
CFDEM_ADD_COMPTURBMOD_LIBS = \
-lcompressibleRASModels \
-lcompressibleLESModels \
-lfluidThermophysicalModels \
