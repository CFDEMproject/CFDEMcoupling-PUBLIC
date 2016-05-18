fluidScalarFieldInterpolator_.reset(interpolation<scalar>::New(propsDict_.lookupOrDefault("fluidScalarFieldInterpolationType",word("cellPoint")),fluidScalarField_).ptr());
