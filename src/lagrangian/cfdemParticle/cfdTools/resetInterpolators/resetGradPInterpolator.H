gradPInterpolator_.reset(interpolation<vector>::New(propsDict_.lookupOrDefault("gradPInterpolationType",word("cellPointFace")),gradP_).ptr());
