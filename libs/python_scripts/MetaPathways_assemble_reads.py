#!/usr/bin/python
# File created on 27 Jan 2012.
from __future__ import division

__author__ = "Kishori M Konwar, Niels W Hanson"
__copyright__ = "Copyright 2014, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import os
    import sys
    import subprocess
    import re
    import gzip
    from optparse import OptionParser
#from libs.python_modules.parsers.fastareader import FastaRecord, FastaReader 

except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     sys.exit(3)

PATHDELIM = os.pathsep 

usage= sys.argv[0] + " -1 <R1_file.fastq.gz> -2 <R2_file.fastq.gz> --12 <interleaved.fastq.gz> --algos <spades>" +\
       "-o <output_dir> --log_file logfile.log " 

parser = None
def createParser():
    global parser
    epilog = """ This script assembles with the assembly algorithms specified. """

    parser = OptionParser(usage=usage, epilog=epilog)

    parser.add_option("-1", "--R1", dest="left_reads", default=None,
                      help='left_hand paired-end read file (i.e., R1)')
    parser.add_option("-2", "--R2", dest="right_reads", default=None,
                      help='right-hand paired-end read file (i.e., R2)')
    parser.add_option("--12", dest="both_reads", default=None,
                      help='interleaved (both left-hand and right-hand) read file')
    parser.add_option("-i", "--input", dest="input_reads", default=None,
                      help='simple non paired-end read file')
    parser.add_option("-l", "--log_file", dest="log_file",  
                      help='file name to write the statsitics and log into')
    parser.add_option("--algos", dest="algos", default="spades", action="append",
                      help="algorithm to assemble with")
    parser.add_option("-m", "--mem-limit", dest="mem_limit", default=None,
                      help="memory limit for assembly algorithms (if applicable)")
    parser.add_option("-o", "--output-dir", dest="output_dir", default=None,
                      help="output file directory to put assemblies") 


def valid_arguments(opts, args):
    state = True

    if ((opts.left_reads != None) != (opts.right_reads != None)) :
        print 'Need to supply both left and right read files'
        state = False
    
    if (opts.both_reads and (opts.left_reads or opts.right_reads)):
        print "Can't supply both interleaved and isolated read sets"
        state = False

    if (opts.input_reads and (opts.left_reads or opts.right_reads) or opts.both_reads):
        print "Can't specify both non paired-end and paired-end inputs"
        state = False
    
    if not opts.output_dir:
        print "Need to specify an output directory"
        state = False
    
    return state


def main(argv, errorlogger = None, runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)

    if not valid_arguments(opts, args):
       print usage
       sys.exit(0)



    return (0,'')

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

