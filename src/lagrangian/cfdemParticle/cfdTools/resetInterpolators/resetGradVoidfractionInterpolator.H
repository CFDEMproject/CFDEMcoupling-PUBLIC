gradVoidfractionInterpolator_.reset(interpolation<vector>::New(propsDict_.lookupOrDefault("gradVoidfractionInterpolationType",word("cellPointFace")),gradVoidfraction_).ptr());
