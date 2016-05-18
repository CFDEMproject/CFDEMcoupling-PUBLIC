gradAlphaInterpolator_.reset(interpolation<vector>::New(propsDict_.lookupOrDefault("gradAlphaInterpolationType",word("cellPointFace")),gradAlpha_).ptr());
