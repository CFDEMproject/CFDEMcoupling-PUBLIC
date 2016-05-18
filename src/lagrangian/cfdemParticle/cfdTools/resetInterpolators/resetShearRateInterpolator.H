shearRateInterpolator_.reset(interpolation<symmTensor>::New(propsDict_.lookupOrDefault("shearRateInterpolationType",word("cellPointFace")),shearRate_).ptr());
