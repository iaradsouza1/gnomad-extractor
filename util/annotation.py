import pickle
from gzip import open
import pandas as pd
import os
from collections import defaultdict

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def get_annotation(index_snp, vcf_file, directory, file_name):

    """
    Extract annotation from LoF SNPs

    index_snp: VCF index file created previously
    vcf_file: gnomad VCF file (gzipped)
    directory: where to save results
    file_name: file name
    """

    # Create 'annotation' directory
    save_path = str(directory + '/annotation')
    access_rights = 0o755
    if not os.path.exists(save_path):
        os.makedirs(save_path, access_rights)

    # Load index
    index = pickle.load(open(index_snp, 'rb'))

    with open(vcf_file, 'rb') as file:

        # Define a list of dictionaries:
        # Each key represents a column in the final dataframe, and they hold a list of values.
        # This dictionary contains only the transcripts annotated for loss-of-functions mutations.
        # Each line on the final dataframe is a trancript with a LoF mutation.
        # For a given position, many different annotations may exist for
        aux_dict = defaultdict(list)

        for line_num, line in enumerate(file, 1):

            line = str(line, 'utf-8')

            if line.startswith('#') or line.startswith('##'):
                continue

            elif line_num in index.keys():

                line = line.strip()
                line = line.split('\t')

                chrom = str(line[0].strip())
                pos = str(line[1].strip())
                ref = str(line[3].strip())
                alt = str(line[4].strip())
                info = str(line[7].strip())

                info_fields = info.split(';')


                for i in info_fields:

                    f = i.strip().split('=')

                    if f[0] == 'vep':
                        annotations = f[1].split(',')

                        for transcripts in annotations:
                            ann = transcripts.split('|')

                            aux_dict['line_num'].append(line_num)
                            aux_dict['chrom'].append(chrom)
                            aux_dict['pos'].append(pos)
                            aux_dict['ref'].append(ref)
                            aux_dict['alt'].append(alt)

                            aux_dict['allele'].append(ann[0])
                            aux_dict['consequence'].append(ann[1])
                            aux_dict['impact'].append(ann[2])
                            aux_dict['symbol'].append(ann[3])
                            aux_dict['feature_type'].append(ann[5])
                            aux_dict['feature'].append(ann[6])
                            aux_dict['biotype'].append(ann[7])
                            aux_dict['cDNA_position'].append(ann[12])
                            aux_dict['cds_position'].append(ann[13])
                            aux_dict['protein_postition'].append(ann[14])
                            aux_dict['codons'].append(ann[16])
                            aux_dict['existing_variation'].append(ann[17])
                            aux_dict['allele_number'].append(ann[18])
                            aux_dict['distance'].append(ann[19])
                            aux_dict['variant_class'].append(ann[22])
                            aux_dict['lof'].append(ann[64])
                            aux_dict['lof_filter'].append(ann[65])
                            aux_dict['lof_flags'].append(ann[66])


         # Create final dataframe
        cols = ['line_num', 'chrom', 'pos', 'ref', 'alt', 'allele', 'consequence',
                'impact', 'symbol', 'feature_type', 'feature', 'biotype', 'cDNA_position',
                'cds_position', 'protein_postition', 'codons', 'existing_variation',
                'allele_number', 'distance', 'variant_class', 'lof', 'lof_filter', 'lof_flags', 'singleton']
        df = pd.DataFrame(aux_dict, columns=cols)

        # Save dataframe
        df.to_csv(str(save_path + '/' + file_name + '_annotation.csv'), mode='w', index=False)
        print(f"Table counts {str(file_name + '_annotation')} was successfully created!")













