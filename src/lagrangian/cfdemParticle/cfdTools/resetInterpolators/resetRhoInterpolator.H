RhoInterpolator_.clear();
RhoInterpolator_.reset(interpolation<scalar>::New(propsDict_.lookupOrDefault("RhoInterpolationType",word("cellPoint")),rho_).ptr());
