package conifer;

import conifer.ctmc.expfam.ExpFamMixture;

public enum StandardEvolutionaryModels {
	KIMURA1980 (ExpFamMixture.kimura1980()),
	DNAGTR (ExpFamMixture.dnaGTR());
	
	private final ExpFamMixture expFamily;
	
	public ExpFamMixture getExpFamily() {
		return expFamily;
	}
	
	StandardEvolutionaryModels(ExpFamMixture expFamMixture) {
		this.expFamily = expFamMixture;
	}
}
