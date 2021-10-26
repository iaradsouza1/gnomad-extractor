import pickle
import re
from gzip import open

def create_index(vcf_file):
    '''
    Creates an index to the gnomad vcf file, considering the position of loss-of-function mutations.
    Here, we index based on the 65th field from CSQ field, from the INFO column, which contains such information.

    :param vcf_file: the VCF file from gnomAD
    :return: It will create an index binary file '.exaci'
    '''

    index_snp = {}

    with open(vcf_file, 'rb') as file:

        for line_num, line in enumerate(file, 1):

            line = str(line, 'utf-8').strip()

            if line.startswith('##') or line.startswith('#'):
                continue

            else:

                fields = line.split('\t')
                filters = fields[6].strip()

                # Filter lines with the 'PASS' tag
                if filters == 'PASS':

                    vep_capture = re.search('(;vep=)(.+)', fields[7])
                    vep_capture = vep_capture.group(2)

                    line_annotation = vep_capture.split(',')

                    annot_list = []

                    for annot in line_annotation:

                        annot = annot.split('|')

                        # Gnomad new version (2.1.1) still has 68 fields on VEP annotation field, just like the previous version
                        if len(annot) == 68:

                            variant_number = int(annot[18])
                            lof = annot[64]

                            if lof == 'HC':
                                if variant_number in annot_list:
                                    pass
                                else:
                                    annot_list.append(variant_number)

                            if len(annot_list) > 0:
                                index_snp[line_num] = annot_list


        index_file_name = vcf_file + '.exaci'

        pickle.dump(index_snp, open(index_file_name, 'wb'))

        print(f'Index file {index_file_name} was successfully created!')

