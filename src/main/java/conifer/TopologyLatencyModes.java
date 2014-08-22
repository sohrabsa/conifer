package conifer;

public enum TopologyLatencyModes {
	DEFAULT (false, false),	// both topology and branch lengths are latent variables 
	FIXED_TOPOLOGY_LATENT_BRANCH_LENGTHS (true, false), 
	FIXED_TOPOLOGY_FIXED_BRANCH_LENGHTS (true, true);
	
	private final boolean fixedTopology;
	private final boolean fixedBranchLength;
	
	TopologyLatencyModes(boolean fixedTopology, boolean fixedBranchLength) {
		this.fixedTopology = fixedTopology;
		this.fixedBranchLength = fixedBranchLength;
	}

	public boolean isFixedTopology() {
		return fixedTopology;
	}

	public boolean isFixedBranchLength() {
		return fixedBranchLength;
	}
	
	public static TopologyLatencyModes initFromBooleans(boolean fixedTopology, boolean fixedBranchLength) {
		if (!fixedBranchLength && !fixedTopology) {
			return TopologyLatencyModes.DEFAULT;
		} else if (!fixedBranchLength && fixedTopology) {
			return TopologyLatencyModes.FIXED_TOPOLOGY_LATENT_BRANCH_LENGTHS;
		} else if (fixedBranchLength && fixedTopology) {
			return TopologyLatencyModes.FIXED_TOPOLOGY_FIXED_BRANCH_LENGHTS;
		} else {
			throw new RuntimeException("Unsupported topology and branch length configuration!");
		}
	}
}
