TInterpolator_.reset(interpolation<scalar>::New(propsDict_.lookupOrDefault("TInterpolationType",word("cellPoint")),T_).ptr());
