# SCRIPT FOR THE GENERATION OF THE DATASET NEEDED TO
# PLOT THE PHYLOGENETIC TREE OF THE GENES INVOLVED IN
# OUR RECONSTRUCTION OF FLAVONOID BIOSYNTHESIS:

import pandas as pd
import numpy as np

def main():
	#IMPORT THE MODEL
	path_to_dir = ''
	file_name = 'modeloFR_2020_v4.0.xlsx'
	path_to_file = ''.join((path_to_dir, file_name))
	#READ SHEETS COTAINING INFO OF GENES AND REACTIONS
	genes_df = pd.read_excel(path_to_file, sheet_name="gene sequence", engine='openpyxl')
	modules_df = pd.read_excel(path_to_file, sheet_name="Reaction List", engine='openpyxl')
	genes_df = genes_df.dropna(how='all')
	moodules_df = genes_df.dropna()
	print(genes_df)
	print(modules_df)
	species = genes_df['Organism'].values.tolist()
	genes = genes_df['gene id'].values.tolist()
	ec = genes_df['CLASS'].values.tolist()

	modules = modules_df['Abbreviation'].values.tolist()
	module_genes = modules_df['Genes'].values.tolist()
	merged_modules = []
	index = 0

	for module in modules:
		merged_modules.append(str(module)+str(module_genes[index]))
		index += 1

	precursor = []
	assembly_chasis = []
	assembly = []
	diversification = []
	decoration_chasis = []
	decoration = []

	#FOR EACH GENE COUNT ITS PARTICIPATION IN EACH OF THE MODULES BY
	#A SUM OF BOOLEAN (0 OR 1) APPARITIONS IN THE VARIABLE 'merged_modules'
	#WHICH IS A LIST OF STRINGS WITH THE FOLLOWING PATTERN :
	#
	#					{RXN_ABBREVIATION}{RXN_GPR}
	#
	#EVENTUALLY SAVES THE VARIABLES IN A DATAFRAME
	for gene in genes:
		precursor.append(sum(gene in s for s in merged_modules if 'PR_' in s))
		assembly_chasis.append(sum(gene in s for s in merged_modules if 'AS_C_' in s))
		assembly.append(sum(gene in s for s in merged_modules if 'AS_' in s) - sum(gene in s for s in merged_modules if 'AS_C_' in s))
		diversification.append(sum(gene in s for s in merged_modules if 'DI_' in s))
		decoration_chasis.append(sum(gene in s for s in merged_modules if 'DE_C_' in s))
		decoration.append(sum(gene in s for s in merged_modules if 'DE_' in s) - sum(gene in s for s in merged_modules if 'DE_C_' in s))

	df = pd.DataFrame({"organism":species,
					   "precursor":precursor,
					   "assembly_chasis":assembly_chasis,
					   "assembly":assembly,
					   "diversification":diversification,
					   "decoration_chasis":decoration_chasis,
					   "decoration":decoration})

	precursor = []
	assembly_chasis = []
	assembly = []
	diversification = []
	decoration_chasis = []
	decoration = []

	oxidorreductases = np.zeros(len(np.unique(species)))
	transferases = np.zeros(len(np.unique(species)))
	hydrolases = np.zeros(len(np.unique(species)))
	lyases = np.zeros(len(np.unique(species)))
	isomerases = np.zeros(len(np.unique(species)))
	ligases = np.zeros(len(np.unique(species)))

	index = 0

	#COMPUTE THE PRESENCE OR ABSENCE OF ECs WITHIN A GIVEN ORGANISM BY
	#PERFORMING A BILEVEL FILTERING OF THE DATAFRAME CONTAINING GENE INFO
	for organism in np.unique(species):
		precursor.append(df.loc[df['organism'] == organism, 'precursor'].sum())
		assembly_chasis.append(df.loc[df['organism'] == organism, 'assembly_chasis'].sum())
		assembly.append(df.loc[df['organism'] == organism, 'assembly'].sum())
		diversification.append(df.loc[df['organism'] == organism, 'diversification'].sum())
		decoration_chasis.append(df.loc[df['organism'] == organism, 'decoration_chasis'].sum())
		decoration.append(df.loc[df['organism']== organism, 'decoration'].sum())
		if genes_df.loc[(genes_df['CLASS'] == 1) & (genes_df['Organism']== organism), 'CLASS'].sum() > 0:
			oxidorreductases[index] = 1

		if genes_df.loc[(genes_df['CLASS'] == 2) & (genes_df['Organism']== organism), 'CLASS'].sum() > 0:
			transferases[index] = 1

		if genes_df.loc[(genes_df['CLASS'] == 3) & (genes_df['Organism']== organism), 'CLASS'].sum() > 0:
			hydrolases[index] = 1

		if genes_df.loc[(genes_df['CLASS'] == 4) & (genes_df['Organism']== organism), 'CLASS'].sum() > 0:
			lyases[index] = 1

		if genes_df.loc[(genes_df['CLASS'] == 5) & (genes_df['Organism']== organism), 'CLASS'].sum() > 0:
			isomerases[index] = 1

		if genes_df.loc[(genes_df['CLASS'] == 6) & (genes_df['Organism']== organism), 'CLASS'].sum() > 0:
			ligases[index] = 1

		index += 1

	#GENERATE THE FINAL DATAFRAME CONTAINING ALL VARIABLES OF INTEREST
	organism_df = pd.DataFrame({"organism":np.unique(species),
								"precursor":precursor,
								"assembly_chasis":assembly_chasis,
								"assembly":assembly,
								"diversification":diversification,
								"decoration_chasis":decoration_chasis,
								"decoration":decoration,
								"oxidorreductases":oxidorreductases,
								"transferases":transferases,
								"hydrolases":hydrolases,
								"lyases":lyases,
								"isomerases":isomerases,
								"ligases":ligases})

	organism_df.to_csv('rawdata_for_phylogenetic_tree.csv', sep ='\t')


if __name__ == '__main__':
	main()
