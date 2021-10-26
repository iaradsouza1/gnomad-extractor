# gnomad-extractor
Tool to parse VCF from gnomad project and extract LoF SNPs.
From the VCF file processed and provided by gnomad Project, extract functional annotation. 

3 Steps:
- Create index with `index`;
- Extract LoF frequencies with `extract`;
- Extract LoF annotations with `annotation`

Example:
`python gnomad-extractor.py index --help `
