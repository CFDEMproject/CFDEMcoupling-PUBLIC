phiP1Interpolator_.reset(interpolation<scalar>::New(propsDict_.lookupOrDefault("phiP1InterpolationType",word("cellPoint")),phiP1_).ptr());
