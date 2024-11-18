SynGO data downloaded on 2024-07-06 17:57 @ https://www.syngoportal.org
SynGO dataset version/release: 20231201

Summary of your gene list
D1R-only
115 / 425 genes from your gene list were mapped to 115 unique SynGO annotated genes.
experimental-evidence filtering is not enabled (default setting)
107 genes have a Cellular Component annotation, 87 for Biological Processes. Note: a gene may have multiple annotations.
186 annotations against a Cellular Component, 151 for Biological Processes.
25 Cellular Component terms are significantly enriched at 1% FDR (testing terms with at least three matching input genes), 39 for Biological Processes.
"your custom background list" background set was selected, contains 34065 unique genes in total of which 1591 overlap with SynGO annotated genes.

"user_genelist_map_to_syngo_genes.xlsx" contains your input genelist (first column) together with the matching gene ID available in SynGO. If your input could not be matched to any SynGO annotated gene, the latter will be empty.
"syngo_annotations_matching_user_input.xlsx" contains all available SynGO annotations that match your input genelist.
"syngo_ontologies_with_annotations_matching_user_input.xlsx" contains all SynGO ontology terms, together with the Gene Set Enrichment Analysis (GSEA) results and all genes from these terms that match your input genelist.

The JSON folder contains data in a format convenient for bioinformatic analysis.
The SynGO_geneset_CC/BP.json files contain all SynGO annotated genes for each synaptic ontology term. On the first level of this nested list you can find 'direct' and 'aggregate' annotations, you will want to use the latter for geneset analyses as these contain (for each term) both the annotated genes for an ontology term and all genes annotated against (recursive) child terms. For instance, in the aggregate dataset the term 'presynapse' also contains all genes annotated against its child term 'active zone'. Further, note that we originally mapped the annotated proteins from various species to HGNC identifiers. The ensembl/entrez mappings (from HGNC ID) in these JSON files were provided by the genenames.org webservice, so if you focus on Ensembl genes consider mapping the HGNC IDs to the exact Ensembl build that you are using.