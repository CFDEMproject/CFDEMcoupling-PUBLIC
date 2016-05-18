alphaInterpolator_.reset(interpolation<scalar>::New(propsDict_.lookupOrDefault("alphaInterpolationType",word("cellPoint")),alpha_).ptr());
