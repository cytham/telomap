import os
import sys
import argparse
import distutils.spawn
from telomap import __version__

# Parse input arguments
def get_args(args=sys.argv[1:]):

    parser = argparse.ArgumentParser(description="Telomap is a tool for analyzing telomeres from telobait-captured or WGS long-read sequencing data",
                                     formatter_class=argparse.RawTextHelpFormatter, usage=msg(), add_help=False)

    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    
    required.add_argument("reads", type=str,
                          metavar="[READS]",
                          help="path to input reads (fasta/fastq/bam/pacbio-bam")

    required.add_argument("dir", type=str,
                          metavar="[WORK_DIRECTORY]",
                          help="path to work directory")

    optional.add_argument("-c", "--capoligo", type=str, metavar="path",
                          help="path to capture oligo fasta file")

    optional.add_argument("-b", "--barcodes", type=str, metavar="path",
                          help="path to barcodes fasta file")
    
    optional.add_argument("-o", "--outprefix", type=str, metavar="path",
                          default="sample",
                          help="prefix of output files [sample]")
    
    def restrict_threads(t):
        t = int(t)
        if t < 1:
            raise argparse.ArgumentTypeError("Number of threads specified < 1, minimum requirement is 1 thread.")
        return t
    
    optional.add_argument("-t", "--threads", type=restrict_threads, metavar="int",
                          default=1,
                          help="specify number of threads [1]")

    optional.add_argument("-v", "--version", action='version',
                          version=__version__,
                          help="print version")
    
    optional.add_argument("-q", "--quiet", action='store_true',
                          help="hide verbose")
  
    optional.add_argument("-h", "--help", action="help",
                          default=argparse.SUPPRESS,
                          help="show this help message and exit")
  
    args = parser.parse_args(args)
    
    # args.insfa, args.suptsv = check_files(args.insfa, args.suptsv, args.dir)
    return args

# Custom usage message
def msg():
    return "telomap [options] [READS] [WORK_DIRECTORY]"

# # Check paths and executables
# def check_exe(path, exe):
#     if path is None:
#         if distutils.spawn.find_executable(exe):
#             return exe
#         else:
#             raise Exception("Error: %s executable is not in PATH" % exe)
#     else:
#         if distutils.spawn.find_executable(path):
#             return path
#         else:
#             raise Exception("Error: %s path do not exist" % path)

# # Check file paths
# def check_files(insfa, suptsv, wk_dir):
#     if insfa is None:
#         insfa = os.path.join(wk_dir, 'ins_seq.fa')
#     if suptsv is None:
#         suptsv = os.path.join(wk_dir, 'sv_support_reads.tsv')
#     if not os.path.isfile(insfa):
#         raise Exception("Error: ins_seq.fa file is not found in %s." % insfa)  
#     if not os.path.isfile(suptsv):
#         raise Exception("Error: sv_support_reads.tsv file is not found in %s." % suptsv)
#     return insfa, suptsv
