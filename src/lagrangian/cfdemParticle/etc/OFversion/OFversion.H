//#define versionExt32 // being tested
//#define version30 // being tested
//#define version24Dev // uses Dev (e.g. new turbulence model structure) code within 24x
#define version24  // 2.4.x with old turbulence model structure  // currently being used
//#define version23
//#define version22
//#define version21
//#define version16ext
//#define version15

//define anisotropicRotation cloud models
//#define anisotropicRotation

// features of 16ext work also in extend 3.2
#if defined(versionExt32)
    #define version16ext
#endif

// features of 2.4Dev work also in Dev
#if defined(version30)
    #define version24Dev
#endif

// basically use 24x settings + some dev features (e.g. new turbulence model structure)
#if defined(version24Dev)
    #define version24
#endif

// features of 2.4 work also in 2.3
#if defined(version24)
    #define version23
#endif

// features of 2.1 work also in 2.3
#if defined(version23)
    #define version21
    #define version221
#endif

// features of 2.1 work also in 2.2
#if defined(version22)
    #define version21
    #define version221
#endif
