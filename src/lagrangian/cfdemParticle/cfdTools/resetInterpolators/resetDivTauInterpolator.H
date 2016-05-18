divTauInterpolator_.reset(interpolation<vector>::New(propsDict_.lookupOrDefault("divTauInterpolationType",word("cellPointFace")),divTau_).ptr());
