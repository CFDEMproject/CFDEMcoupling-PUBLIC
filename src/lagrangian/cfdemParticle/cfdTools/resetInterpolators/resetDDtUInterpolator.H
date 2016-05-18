DDtUInterpolator_.reset(interpolation<vector>::New(propsDict_.lookupOrDefault("DDtUInterpolationType",word("cellPointFace")),DDtU_).ptr());
