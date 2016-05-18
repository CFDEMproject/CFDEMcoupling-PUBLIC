phiP2Interpolator_.reset(interpolation<scalar>::New(propsDict_.lookupOrDefault("phiP2InterpolationType",word("cellPoint")),phiP2_).ptr());
