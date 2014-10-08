import os
import subprocess

#post-processing code for RainMaker Report integration
import lib.parse_into_tables as parser

#Base class which handles download, execution and upload
#By inheriting, we only have to supply the details of the app's execution, not the plumbing
from lib.docker_app import DockerApp

CONFIG_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "config")

class MyApp( DockerApp ):
    def __init__(self) :
        pass

    #Input: A list of files that were downloaded.
    #Output: A list of directories containing the output of your app's execution.
    #        Each folder in the list is turned into an AppResult and every file
    #        inside is uploaded into that AppResult, mirroring the directory structure to the server.
    #        The AppResult's name is that of the last folder in the path

    def parseMp2Input(self, app_session):
        mp_parsed_params = {}
        if app_session:
            for i in app_session["Properties"]["Items"]:
                if "Name" in i.keys() and "Input" in i["Name"]:
                    if i["Type"] == "string":
                        name = i["Name"].lstrip("Input.")
                        mp_parsed_params[name] = i["Content"]
                    if i["Type"] == "string[]":
                        name = i["Name"].lstrip("Input.")
                        mp_parsed_params[name] = i["Items"]

        return mp_parsed_params

    def createMp2ParameterFile(self, mp2_parsed_params):
        if mp2_parsed_params:
            param_fh = open(os.path.join(CONFIG_DIR, "template_param.txt"), "w")
            preamble = ["##V.1   do not remove this line", \
                        "# MetaPathways v3.0", \
                        "# Kishori M. Konwar, Niels W. Hanson", \
                        "# Parameter File", \
                        "INPUT:format fasta"]
            # QC params
            qc_string = ["# Quality Control  parameters"]
            if "qc-min-nuc-size" in mp2_parsed_params:
                qc_string.append("quality_control:min_length " + \
                                 mp2_parsed_params["qc-min-nuc-size"])
            else:
                # default values
                qc_string.append("quality_control:min_length 180")

            if "qc-nuc-remove-dups" in mp2_parsed_params:
                if mp2_parsed_params["qc-nuc-remove-dups"] == "1":
                    qc_string.append("quality_control:delete_replicates yes")
                else:
                    qc_string.append("quality_control:delete_replicates no")

            orf_pred_string = ["# ORF prediction parameters"]

            if "orf-predict-algorithm" in mp2_parsed_params:
                orf_pred_string.append("orf_prediction:algorithm " + \
                                       mp2_parsed_params["orf-predict-algorithm"])
            else:
                orf_pred_string.append("orf_prediction:algorithm prodigal")

            if "orf-predict-min-aa-size" in mp2_parsed_params:
                orf_pred_string.append("orf_prediction:min_length " + \
                                       mp2_parsed_params["orf-predict-min-aa-size"])
            else:
                orf_pred_string.append("orf_prediction:min_length 60")

            if "orf-predict-translation-table" in mp2_parsed_params:
                orf_pred_string.append("orf_prediction:translation_table " + \
                                       mp2_parsed_params["orf-predict-translation-table"])
            else:
                orf_pred_string.append("orf_prediction:translation_table 11")

            # write out the rest of the parameters 
            post_string = ["# ORF annotation parameters",
                           "annotation:algorithm LAST",
                           "# e.g. blast or last",
                           "annotation:dbs metacyc-v4-2011-07-03,COG_2013-12-27",
                           "annotation:min_bsr 0.4",
                           "annotation:max_evalue 0.000001",
                           "annotation:min_score 20",
                           "annotation:min_length 60",
                           "annotation:max_hits 5",
                           "# rRNA annotation parameters",
                           "rRNA:refdbs LSURef_115_tax_silva",
                           "rRNA:max_evalue 0.000001",
                           "rRNA:min_identity 20",
                           "rRNA:min_bitscore 50",
                           "# pathway tools parameters",
                           "ptools_settings:taxonomic_pruning no"]

            execution_string = """# pipeline execution flags
metapaths_steps:PREPROCESS_INPUT yes
metapaths_steps:ORF_PREDICTION yes
metapaths_steps:ORF_TO_AMINO yes
metapaths_steps:FILTER_AMINOS yes
metapaths_steps:COMPUTE_REFSCORES yes
metapaths_steps:FUNC_SEARCH skip
metapaths_steps:PARSE_FUNC_SEARCH skip
metapaths_steps:SCAN_rRNA skip
metapaths_steps:SCAN_tRNA skip
metapaths_steps:ANNOTATE_ORFS skip
metapaths_steps:PATHOLOGIC_INPUT skip
metapaths_steps:GENBANK_FILE skip
metapaths_steps:CREATE_ANNOT_REPORTS skip
metapaths_steps:MLTREEMAP_CALCULATION skip
metapaths_steps:BUILD_PGDB skip
metapaths_steps:COMPUTE_RPKM skip"""

            # write to template_param.txt file
            param_fh.write("\n".join(preamble) + "\n")
            param_fh.write("\n".join(qc_string) + "\n")
            param_fh.write("\n".join(orf_pred_string)+ "\n")
            param_fh.write("\n".join(post_string) + "\n")
            param_fh.write(execution_string)
            param_fh.close()

    def createSimpleMp2Command(self, input, output, config, param):

        mp2_command = ['python', 'MetaPathways.py',\
                       '-i', input,
                       '-o', output,
                       '-c', config,
                       '-p', param,
                       '-v', '-r', 'overlay']
        exe = subprocess.call( mp2_command )
        if exe != 0 : raise Exception("MetaPathways exited abnormally")

    def runApp( self, downloaded_files, app_session = None ) :
        upload_folders = []
        mp2_parsed_params = self.parseMp2Input(app_session)
        self.createMp2ParameterFile(mp2_parsed_params)

        for file in downloaded_files :
            if not (file.endswith( ".fastq" ) or \
                    file.endswith( ".fastq.gz") or \
                    file.endswith(".fas") ) : continue

            output_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "output")
            config = os.path.join(os.path.dirname(os.path.realpath(__file__)), "config")
            template_config = os.path.join(config, "template_config.txt")
            template_param = os.path.join(config, "template_param.txt")
            self.createSimpleMp2Command(file, output_dir, template_config, template_param)
            exit()
            # check if output folder exists
            #print "making directory"
            #if not os.path.exists( output_dir ) :
            #    os.makedirs( output_dir )

            # Assemble the file
            # assemble_out = "/assembly"
            # assemble_out = os.path.join(output_dir, "assembly", "spades")
            # print "Making Spades directory:"
            #if not os.path.exists ( assemble_out ):
            #    os.makedirs( assemble_out )
            # assemble_command = [ 'python', '/home/SPAdes-3.1.1-Linux/bin/spades.py', \
            #                     '--12', file, \
            #                     '-o', assemble_out ]

            print assemble_command
            # execute Spades
            exe = subprocess.call( assemble_command )

            # Try running spaces

            #run the app
            # output_dir = "/home/apps/myapp/output"
            #            command_list = ['perl', \
            #                            '/home/bin/fastqc_v0.10.1_source/FastQC/my_fastqc/fastqc', \
            #                            '--outdir=%s' % output_dir, file]
            #            print "command_list", command_list
            #            rcode = subprocess.call( command_list )
            #            if rcode != 0 : raise Exception("fastqc process exited abnormally")

            #identify the data file for post-processing
            #            (head,tail) = os.path.split( file )
            #            sample_name = tail.replace(".gz",'').replace(".fastq",'')
            #            fastqc_output_folder = sample_name + "_fastqc"
            #            output_sample_folder = os.path.join( output_dir, fastqc_output_folder )
            #            table_file = os.path.join( output_sample_folder, "fastqc_data.txt" )

            #post-processing
            #            parsed_dir = os.path.join( output_sample_folder, "parsed" )
            #            if not os.path.exists( parsed_dir ) :
            #                os.mkdir( parsed_dir )
            #            parser.parseIntoTables( table_file, parsed_dir )

            upload_folders.append( assemble_out )

        return upload_folders

#the entry point
if __name__ == "__main__" :
    MyApp().run()

