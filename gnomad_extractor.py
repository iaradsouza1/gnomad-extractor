import argparse
from util.index import *
from util.extract import *
from util.annotation import *

def main():
    parser = argparse.ArgumentParser(description='A tool to filter LoF mutations from gnomAD VCF file')

    subparser = parser.add_subparsers(title='commands',
                                      description='subcommands',
                                      dest='commands',
                                      help='Add a {command} -h for help')

    # Index
    parser_index = subparser.add_parser('index', help='Create an index from gnomAD VCF file')
    parser_index.add_argument('-f', '--file', dest='gnomad_vcf', action='store', default=None, required=False)

    # Extract counts
    parser_extract = subparser.add_parser('extract', help='Extract variants')
    parser_extract.add_argument('-i', '--index', dest='index', default=None, required=True, action='store',
                                help='Index file created with the index command')
    parser_extract.add_argument('-f', '--file', dest='file', default=None, required=True, action='store',
                                help='gnomAD VCF file')
    parser_extract.add_argument('-s', '--subset', dest='subsets', nargs='+', default=None, required=True, action='store',
                                help='one or more of the following subsets: ["controls", "non_topmed", "non_neuro", "non_cancer"]')
    parser_extract.add_argument('-d', '--directory', dest='directory', default=None, required=True, action='store',
                                help='path to create "chrom/" directory')
    parser_extract.add_argument('-o', '--output', dest='fname', default=None, required=True, action='store',
                                help='result file name')

    # Print lines
    parser_print= subparser.add_parser('print', help='Print VCf position')
    parser_print.add_argument('-f', '--file', dest='file', action='store', default=None, required=False)
    parser_print.add_argument('-i', '--index', dest='index', default=None, required=True, action='store',
                                help='Index file created with the  index command')
    parser_print.add_argument('-p', '--position', dest='print_pos', action='store', default=None, required=False)

    # Extract annotation
    parser_annotation = subparser.add_parser('annotation', help='get annotations from LoF mutations')
    parser_annotation.add_argument('-f', '--file', dest='file', action='store', default=None, required=False)
    parser_annotation.add_argument('-i', '--index', dest='index', default=None, required=True, action='store',
                                help='Index file created with the  index command')
    parser_annotation.add_argument('-d', '--directory', dest='directory', default=None, required=True, action='store',
                                help='path to create "annotation/" directory')
    parser_annotation.add_argument('-o', '--output', dest='fname', default=None, required=True, action='store',
                                help='result file name')


    args = parser.parse_args()

    if args.commands == 'index':
        create_index(args.gnomad_vcf)
    if args.commands == 'extract':
        extract_freq(args.index, args.file, args.subsets, args.directory, args.fname)
    if args.commands == 'annotation':
        get_annotation(args.index, args.file, args.directory, args.fname)


if __name__ == '__main__':
    main()
