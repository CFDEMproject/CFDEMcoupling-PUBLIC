# paths for additional libraries
CFDEM_ADD_LIB_PATHS = \

# additional libraries to be linked to solvers
CFDEM_ADD_LIBS = \

# additional static libraries to be linked to lagrangian library
CFDEM_ADD_STATICLIBS = \
-lmpi_cxx \

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

# libs for turbulence models to use
CFDEM_ADD_INCOMPTURBMOD_LIBS = \
-lturbulenceModels \
-lincompressibleTurbulenceModels \
-lfvOptions \

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
