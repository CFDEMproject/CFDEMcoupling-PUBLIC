voidfractionInterpolator_.reset(interpolation<scalar>::New(propsDict_.lookupOrDefault("voidfractionInterpolationType",word("cellPoint")),voidfraction_).ptr());
