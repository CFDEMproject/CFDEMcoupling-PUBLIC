#if defined(MS)
    #if defined(superquadrics_flag)
        cfdemCloudRotationSuperquadric particleCloud(mesh);
    #else
        cfdemCloudMS particleCloud(mesh);
    #endif
#else
    cfdemCloud particleCloud(mesh);
#endif
