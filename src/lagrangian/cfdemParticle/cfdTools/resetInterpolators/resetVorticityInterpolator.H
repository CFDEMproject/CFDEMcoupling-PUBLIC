vorticityInterpolator_.reset(interpolation<vector>::New(propsDict_.lookupOrDefault("vorticityInterpolationType",word("cellPointFace")),vorticity_).ptr());
