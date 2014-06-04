#define version23  // currently being tested
//#define version22  // currently being used
//#define version21
//#define version16ext
//#define version15

// choose comp/incomp
//#define comp   // if comp is on - you must use Make/options_comp!

//define multi sphere
#define multisphere

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
