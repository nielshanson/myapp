#!/usr/bin/python

"""This script run the pathologic """

try:
   import optparse, sys, re, csv, traceback
   from os import path, _exit
   import logging.handlers
   from glob import glob

   from libs.python_modules.utils.sysutil import pathDelim
   from libs.python_modules.utils.metapathways_utils  import fprintf, printf, eprintf,  exit_process
   from libs.python_modules.utils.sysutil import getstatusoutput

   from libs.python_modules.utils.pathwaytoolsutils import *

except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     print traceback.print_exc(10)
     sys.exit(3)


PATHDELIM= pathDelim()



def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)

def files_exist( files , errorlogger = None):
    status = True    
    for file in files:
       if not path.exists(file):
          if errorlogger:
             errorlogger.write( 'ERROR\tCould not find ptools input  file : ' +  file )
          status = False
    return not status



usage = sys.argv[0] + """ -c <contigs> -o <output> -r <reads>  -O <orfgff> --rpkmExec <rpkmexec> """
parser = None
def createParser():
    global parser

    epilog = """This script computes the RPKM values for each ORF, from the BWA 
                recruits. 
             """

    epilog = re.sub(r'\s+', ' ', epilog)

    parser = optparse.OptionParser(usage=usage, epilog = epilog)

    # Input options


    parser.add_option('-c', '--contigs', dest='contigs', default=None,
                           help='the contigs file')

    parser.add_option('-o', '--output', dest='output', default=None,
                           help='orfwise RPKM file')

    parser.add_option('--stats', dest='stats', default=None,
                           help='output stats for ORFs  into file')

    parser.add_option('-r', '--rpkmdir', dest='rpkmdir', default=None,
                           help='list of sam files that contains the read recruitments')

    parser.add_option('-O', '--orfgff', dest='orfgff', default=None,
                           help='folder of the PGDB')

    parser.add_option('-s', '--sample_name', dest='sample_name', default=None,
                           help='name of the sample')

    parser.add_option('--rpkmExec', dest='rpkmExec', default=None,
                           help='RPKM Executable')


def main(argv, errorlogger = None, runcommand = None, runstatslogger = None):
    global parser

    options, args = parser.parse_args(argv)

    parser.add_option('-r', '--rpkmdir', dest='rpkmdir', default=None,
                           help='list of sam files that contains the read recruitments')

    parser.add_option('-O', '--orfgff', dest='orfgff', default=None,
                           help='folder of the PGDB')

    parser.add_option('-s', '--sample_name', dest='sample_name', default=None,
                           help='name of the sample')

    parser.add_option('--rpkmExec', dest='rpkmExec', default=None,
                           help='RPKM Executable')



def main(argv, errorlogger = None, runcommand = None, runstatslogger = None):
    global parser

    options, args = parser.parse_args(argv)
    if options.contigs ==None:
       parser.error('ERROR\tThe contigs file is missing')

    if options.rpkmExec ==None:
       parser.error('ERROR\tThe RPKM executable is missing')

    if options.rpkmdir ==None:
       parser.error('ERROR\tThe RPKM directory')

    if options.sample_name ==None:
       parser.error('ERROR\tThe sample name is missing')


    # is there a pathwaytools executable installed
    if not path.exists(options.rpkmExec):
       eprintf("ERROR\tRPKM executable %s not found!\n", options.rpkmExec)
       if errorlogger:
          errorlogger.printf("ERROR\tRPKM executable %s not found!\n",  options.rpkmExec)
       exit_process("ERROR\tRPKM executable %s not found!\n" %(options.rpkmExec))


    # command to build the ePGDB
    command = "%s -c %s"  %(options.rpkmExec, options.contigs)
    command += " --multireads --format sam-2" 
    if options.output:
       command += " --ORF-RPKM %s" %(options.output)
       command += " --stats %s" %(options.stats)

    if options.orfgff:
       command += " -O %s" %(options.orfgff)


    samfiles = []
    if path.exists(options.rpkmdir):
       samfiles = glob(options.rpkmdir + PATHDELIM + options.sample_name + '*.sam')
    
    if not samfiles:
       return 

    for samfile in samfiles:
        command += " -r " + samfile


    try:
       status  = runRPKMCommand(runcommand = command) 
    except:
       status = 1
       pass

    if status!=0:
       eprintf("ERROR\tFailed to run RPKM \n")
       exit_process("ERROR\tFailed to run RPKM" )


def runRPKMCommand(runcommand = None):
    if runcommand == None:
      return False

    result = getstatusoutput(runcommand)
    return result[0]


# this is the portion of the code that fixes the name

def split_attributes(str, attributes):
     rawattributes = re.split(';', str)
     for attribStr in rawattributes:
        insert_attribute(attributes, attribStr)

     return attributes


# this is the function that fixes the name
def  fix_pgdb_input_files(pgdb_folder, pgdbs = []):
     pgdb_list = glob(pgdb_folder + '/*/input/organism.dat')     

     for pgdb_organism_file in pgdb_list:
        process_organism_file(pgdb_organism_file)


def fixLine(line, id):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[0]+'\t' + id
     

def getID(line):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[1]
     
def process_organism_file(filel):
     patternsToFix = [ re.compile(r'NAME\tunclassified sequences'), re.compile(r'ABBREV-NAME\tu. sequences') ]
     patternID =  re.compile(r'^ID\t.*')
     try:
         orgfile = open(filel,'r')
     except IOError:
         print "ERROR : Cannot open organism file" + str(filel)
         return 

     lines = orgfile.readlines()
     newlines = []
     needsFixing = False

     id = None
     for line in lines:
         line = line.strip()
         if len(line)==0:
            continue
         flag = False

         result = patternID.search(line)
         if result:   
             id = getID(line)
          
         for patternToFix in patternsToFix:
             result = patternToFix.search(line)
             if result and id:
                 newline = fixLine(line, id)
                 newlines.append(newline)
                 flag= True
                 needsFixing = True

         if flag==False:
            newlines.append(line)

     orgfile.close()
     if needsFixing:
       write_new_file(newlines, filel)


def write_new_file(lines, output_file):
    
    print "Fixing file " + output_file 
    try:
       outputfile = open(output_file,'w')
       pass
    except IOError:
         print "ERROR :Cannot open output file "  + output_file
   
    for line in lines:
       fprintf(outputfile, "%s\n", line)

    outputfile.close()


def MetaPathways_rpkm(argv, extra_command = None, errorlogger = None, runstatslogger =None): 
    if errorlogger != None:
       errorlogger.write("#STEP\tRPKM_CALCULATION\n")
    createParser()
    main(argv, errorlogger = errorlogger, runcommand= extra_command, runstatslogger = runstatslogger)
    return (0,'')

if __name__ == '__main__':
    createParser()
    main(sys.argv[1:])

