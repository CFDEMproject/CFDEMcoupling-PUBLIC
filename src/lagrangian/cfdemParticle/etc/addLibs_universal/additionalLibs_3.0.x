# Specify additional include and library paths, as well as libraries for the compilation
#
#CFDEM_ADD_INC = \
#CFDEM_ADD_LIB_PATHS = \
#CFDEM_ADD_LIBS = \

# additional static libraries to be linked to lagrangian library
CFDEM_ADD_STATICLIBS = \

# include flags for compiling with SQ
include $(CFDEM_ADD_LIBS_DIR)/additionalLibs_superquadric

#################################################################
## SETTINGS FOR 3.0.x                                          ##
#################################################################
#----------------------------------------------------------------
# incompressible turbulence model settings
#----------------------------------------------------------------
# paths for incompressible turbulence models to use
CFDEM_ADD_INCOMPTURBMOD_PATHS = \
-I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
-I$(LIB_SRC)/TurbulenceModels/incompressible/lnInclude \
-I$(LIB_SRC)/fvOptions/lnInclude \
###-I$(LIB_SRC)/thermophysicalModels/radiation/lnInclude \

# libs for turbulence models to use
CFDEM_ADD_INCOMPTURBMOD_LIBS = \
-lturbulenceModels \
-lincompressibleTurbulenceModels \
-lfvOptions \

CFDEM_TRI_SURF = \
-ltriSurface

CFDEM_SPRAY_LIBS = \
    -lliquidProperties \
    -lliquidMixtureProperties \
    -lsolidProperties \
    -lsolidMixtureProperties \
    -lthermophysicalFunctions 

#----------------------------------------------------------------
# compressible turbulence model settings
#----------------------------------------------------------------
# paths for compressible turbulence models to use
CFDEM_ADD_COMPTURBMOD_PATHS = \
-I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
-I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
-I$(LIB_SRC)/transportModels/compressible/lnInclude \
-I$(LIB_SRC)/thermophysicalModels/radiation/lnInclude \

# libs for turbulence models to use
CFDEM_ADD_COMPTURBMOD_LIBS = \
-lturbulenceModels \
-lcompressibleTurbulenceModels \
#################################################################
