gradPsolidInterpolator_.reset(interpolation<vector>::New(propsDict_.lookupOrDefault("gradPsolidInterpolationType",word("cellPoint")),gradPsolid_).ptr());
