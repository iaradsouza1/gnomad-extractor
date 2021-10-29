# gnomad-extractor
Tool to (1) parse the VCF file from gnomad project and (2) specifically extract the loss-of-functions SNPs. It also extracts the functional annotation of each LoF SNP.

3 Steps:
- Create index with `index`;
- Extract LoF frequencies with `extract`;
- Extract LoF annotations with `annotation`

Example:
`python gnomad-extractor.py index --help `
