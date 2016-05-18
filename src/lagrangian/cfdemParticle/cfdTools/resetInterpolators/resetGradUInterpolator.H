gradUInterpolator_.reset(interpolation<vector>::New(propsDict_.lookupOrDefault("gradUInterpolationType",word("cellPointFace")),gradU_).ptr());
