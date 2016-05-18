dSauterInterpolator_.reset(interpolation<scalar>::New(propsDict_.lookupOrDefault("dSauterInterpolationType",word("cellPoint")),dSauter_).ptr());
