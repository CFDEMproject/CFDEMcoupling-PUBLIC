GInterpolator_.reset(interpolation<scalar>::New(propsDict_.lookupOrDefault("GInterpolationType",word("cellPoint")),G_).ptr());
