import pickle
from gzip import open
import pandas as pd
import os

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def extract_freq(index_snp, vcf_file, subset, directory, file_name):

    """
    Extracts the LoF frequencies from the gnomad VCF
    index_snp: VCF index file created previously
    vcf_file: gnomad VCF file (gzipped)
    subset: from which gnomad group (gnomad_all, non_pop, non_cancer, non_topmed, non_neuro)
    extract SNPs frequencies
    directory: where to save results
    file_name: file name 
    """

    # Create 'chrom' directory
    save_path = str(directory + '/chrom')
    access_rights = 0o755
    if not os.path.exists(save_path):
        os.makedirs(save_path, access_rights)

    # Define fields to get from info field
    pop = ['afr', 'amr', 'asj', 'eas', 'fin', 'nfe', 'oth', 'sas']
    group = subset
    count_values = ['AC', 'AN', 'AF', 'nhomalt']

    print("Getting count values columns: ", *count_values)
    print("From the following subsets: ", *subset)

    # Create combinations
    fields_to_get = []
    for g in group:
        for count in count_values:
            for p in pop:
                fields_to_get.append(str(g + '_' + count + '_' + p))

    # Add male and female information
    for i in range(len(fields_to_get)):
        fields_to_get.append(fields_to_get[i] + '_' + 'female')
        fields_to_get.append(fields_to_get[i] + '_' + 'male')

    # Add other important fields
    add_fields = ['n_alt_alleles', 'popmax', 'variant_type']
    for i in range(len(add_fields)):
        fields_to_get.append(add_fields[i])

    # Load index
    index = pickle.load(open(index_snp, 'rb'))

    with open(vcf_file, 'rb') as file:

        # Define a list of dictionaries:
        # Each key represents a column in the final dataframe, and they hold a list of values.
        # This dictionary contains only the transcripts annotated for loss-of-functions mutations.
        # Each line on the final dataframe is a trancript with a LoF mutation.
        # For a given position, many different annotations may exist for a given allele.
        
        ls_of_dicts = []

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
                filter = str(line[6].strip())
                info = str(line[7].strip())

                info_fields = info.split(';')
                aux_dict = {}

                for i in info_fields:
                    f = i.strip().split('=')

                    if f[0]  == 'vep':
                        annotations = f[1].split(',')

                        for transcripts in annotations:
                            ann = transcripts.split('|')
                            if ann[64] != 'HC':
                                continue
                            else:
                                aux_dict['allele'] = ann[0]
                                aux_dict['consequence'] = ann[1]
                                aux_dict['impact'] = ann[2]
                                aux_dict['symbol'] = ann[3]
                                aux_dict['feature_type'] = ann[5]
                                aux_dict['feature'] = ann[6]
                                aux_dict['lof'] = ann[64]
                                aux_dict['lof_filter'] = ann[65]
                                aux_dict['lof_flags'] = ann[66]

                                for i in info_fields:
                                    f = i.strip().split('=')
                                    if len(f) == 1:
                                        aux_dict[f[0]] = True
                                    else:
                                        if any(f[0] in s for s in fields_to_get):
                                            if is_float(f[1]):
                                                aux_dict[f[0]] = float(f[1])
                                            else:
                                                aux_dict[f[0]] = f[1]
                                                aux_dict['chrom'] = chrom
                                                aux_dict['pos'] = pos
                                                aux_dict['ref'] = ref
                                                aux_dict['alt'] = alt
                                                aux_dict['filter'] = filter
                                                aux_dict['line_num'] = line_num

                                ls_of_dicts.append(aux_dict)

        # Create final dataframe
        df = pd.DataFrame(ls_of_dicts)

        # Save dataframe
        df.to_csv(str(save_path + '/' + file_name + '.csv'), mode='w', index=False)
        print(f"Table counts {file_name} was successfully created!")













