#if defined(MS)
    #if defined(superquadrics_flag)
        #include "cfdemCloudRotationSuperquadric.H"
    #else
        #include "cfdemCloudMS.H"
    #endif
#else
    #include "cfdemCloud.H"
#endif
