import os
import subprocess
import re
import gzip
import shutil
from libs.python_modules.parsers.fastareader import FastaRecord, FastaReader

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
                    if i["Type"] == "sample[]":
                        if i["Name"] not in mp_parsed_params:
                            mp_parsed_params[i["Name"]] = []
                        for j in i["Items"]:
                            mp_parsed_params[i["Name"]].append(j["Name"])
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

        mp2_command = ['python', '/home/apps/myapp/MetaPathways.py',\
                       '-i', input,
                       '-o', output,
                       '-c', config,
                       '-p', param,
                       '-v', '-r', 'overlay']
        exe = subprocess.call( mp2_command )
        if exe != 0 : raise Exception("MetaPathways exited abnormally")
        print "Returned!"

    def getSequencePairs(self, downloaded_files):
        paired_samples = {}
        paired_pattern = re.compile("(.+?)\_(R[0-2])\_(.+?)[\.fastq|\.fq][\.gz]?")

        for file in downloaded_files:
            if not (file.endswith( ".fastq" ) or \
                    file.endswith( ".fastq.gz") or \
                    file.endswith(".fas") ) or \
                    file.endswith(".fq") or \
                    file.endswith(".fq.gz"): continue
            hits = paired_pattern.search(file)
            if hits:
                sample = os.path.basename(hits.group(1)) + "_" + hits.group(3)
                if sample not in paired_samples:
                    paired_samples[sample] = {"R1":None, "R2": None}
                if hits.group(2) == "R1":
                    print "Found R1 in ", file
                    paired_samples[sample]["R1"] = file
                elif hits.group(2) == "R2":
                    print "Found R2 in ", file
                    paired_samples[sample]["R2"] = file
                else:
                    print "Error: Did not find paired-end pattern for " + sample
        return paired_samples
    
    def unzip_sample_files(self, sample):
        for r in sample:
            if os.path.exists(sample[r]):
               new_file = sample[r].replace(".gz", "")
               with gzip.open(sample[r], 'rb') as fh_in:
                   with open(new_file, 'wb') as fh_out:
                       fh_out.writelines(fh_in)
               sample[r] = new_file
        return sample

    def assemble_sample(self, sample, sample_pairs, algos, paired = False):

        root = os.path.dirname(os.path.realpath(__file__))
        # folder for resulting contigs
        contigs_folder = os.path.join(root, "assemblies", sample, "assembly", "final_contigs")
        if not os.path.exists( contigs_folder ):
            os.makedirs( contigs_folder )

        if "spades" in algos:
# spades_exe = os.path.join(root, "executables", "Darwin", "spades", "bin", "spades.py")
            spades_exe = os.path.join(root, "executables", "Unix", "SPAdes-3.1.1-Linux", "bin", "spades.py")
            assemble_out = os.path.join(root, "assemblies", sample, "assembly", "spades")
            if not os.path.exists( assemble_out ) :
                os.makedirs( assemble_out )

            if paired:
                ## TODO will have to make commands more robust to handle multiple pairs of command
                ## See: http://spades.bioinf.spbau.ru/release3.0.0/manual.html
                ## TODO memory limit -m needs to be tweeked for system
                assemble_command = [ 'python', spades_exe, \
                                     '-1', sample_pairs["R1"],\
                                     '-2', sample_pairs["R2"],\
                                     '--careful',\
                                     "-m", "2",\
                                     '-o', assemble_out ]
            else:
                assemble_command = [ 'python', spades_exe,\
                                     '-12', sample_pairs["R1"],\
                                     '--careful',\
                                     "-m", "2",\
                                     '-o', assemble_out ]
            
            ret = subprocess.call(assemble_command)
            if ret != 0 : raise Exception("Assembly exited abnormally")

            # copy resulting contig file
            spades_contigs = os.path.join(assemble_out, "contigs.fasta")
            os.rename(spades_contigs, os.path.join(contigs_folder, "spades_contigs.fasta"))

        if "idba-ud" in algos:

            # executables
            exe = os.path.join(root, "executables", "Unix", "idba-1.1.1", "bin", "idba_ud")
            fq2fa_exe = os.path.join(root, "executables", "Unix", "idba-1.1.1", "bin", "fq2fa")

            # output folder for results
            temp_fasta_file = "temp.fa"
            assemble_out = os.path.join(root, "assemblies", sample, "assembly", "idba-ud")
            if not os.path.exists( assemble_out ) :
                os.makedirs( assemble_out )
            
            def assemble_idba_ud(sample_pairs, paired=False):
                """
                    Helper function for to run idba_ud assembly
                """
                if paired:
                    # run fq2fa
                    if len(sample_pairs) > 1:
                        # two paired files (left, right)
                        fq2fa_command = [fq2fa_exe, '--merge', '--filter', sample_pairs["R1"],  sample_pairs["R2"], temp_fasta_file]
                    else:
                        # interlaced files files
                        fq2fa_command = [fq2fa_exe, '--merge', '--filter', sample_pairs["R1"], temp_fasta_file]
                    ret = subprocess.call(fq2fa_command)
                    if ret != 0 : raise Exception("fq2fa_exe exited abnormally")

                    assemble_command = [ exe,\
                                        '-r', temp_fasta_file,\
                                        '--pre_correction',\
                                        '-o', assemble_out ]
                    ret = subprocess.call(assemble_command)
                    if ret != 0 : raise Exception("idba_ud exited abnormally")
                else:
                    if len(sample_pairs) > 1:
                        # two paired files (left, right)
                        fq2fa_command = [fq2fa_exe, '--merge', '--filter', sample_pairs["R1"],  sample_pairs["R2"], temp_fasta_file]
                    else:
                        # interlaced files
                        fq2fa_command = [fq2fa_exe, '--merge', '--filter', sample_pairs["R1"], temp_fasta_file]
                    ret = subprocess.call(fq2fa_command)
                    if ret != 0 : raise Exception("fq2fa_exe exited abnormally")

                    assemble_command = [ exe,\
                                        '-r', temp_fasta_file,\
                                        '--pre_correction',\
                                        '-o', assemble_out ]
                    ret = subprocess.call(assemble_command)
                    if ret != 0 : raise Exception("idba_ud exited abnormally")

            if paired:
                ## TODO will have to make commands more robust to handle multiple pairs of command
                ## See: http://i.cs.hku.hk/~alse/hkubrg/projects/idba_ud/
                
                # check to see if files are .gziped
                if sample_pairs["R1"].endswith(".gz"):
                    sample_pairs_unzip = self.unzip_sample_files(sample_pairs)
                    assemble_idba_ud(sample_pairs_unzip, paired)
                else: 
                    assemble_idba_ud(sample_pairs, paired)
            
            # remove intermediate files
            if sample_pairs_unzip:
               for s in sample_pairs_unzip:
                   if os.path.exists(sample_pairs_unzip[s]):
                       os.remove(sample_pairs_unzip[s])
            os.remove(temp_fasta_file) 
            
            # copy resulting contig file
            idba_ud_contigs = os.path.join(assemble_out, "contig.fa")
            os.rename(idba_ud_contigs, os.path.join(contigs_folder, "idba_ud_contigs.fasta"))

    def getLens(self, filename):
        """
        Parses FASTA file using screed to create a sorted list of contig lengths.
        """
        lens = []
        reader = FastaReader(filename)
        
        for record in reader:
            lens.append(len(record.sequence))
        
        return sorted(lens) 
    
    def nx(self, lens, percent):
        """
        Calculates any NXX (e.g. N50, N90) statistic.
        """
        
        lenSum = sum(lens)
        threshold = (float(percent) / 100) * lenSum
        
        runningSum = 0
        nxx = 0
        nxxLen = 0

        for i in range(len(lens)-1, -1, -1):
            myLen = lens[i]
            nxx += 1
            runningSum += myLen

            if runningSum >= threshold:
                nxxLen = myLen
                break

        return nxxLen

    def compute_assembly_statistics(self, paired_samples, min_length=20):
        
        # write the headers
        header_nx = "Sample\tAssembler\tN_x\tvalue\n"
        header = "Sample\tAssembler\tN\tN_trimmed\tTotal_Length\tMin\tMedian\tMean\tMax\tN50\tN90\tNx_AUC\n"
        root = os.path.dirname(os.path.realpath(__file__))

        for sample in paired_samples:

            # output
            contigs_folder = os.path.join(root, "assemblies", sample, "assembly", "final_contigs") 
            file_list = os.listdir(contigs_folder)
            nx_stats_out = open(os.path.join(contigs_folder, "assembly_stats_nx.txt"),"w")
            stats_out = open(os.path.join(contigs_folder, "assembly_stats.txt"),"w")
            nx_stats_out.write(header_nx)
            stats_out.write(header)

            for f in file_list:
                algo = "a1"
                hits = re.search("(.*)\_contigs\.fasta",f)

                if hits:
                   algo = hits.group(1)
                else:
                    print 'File ', f, ' not assembler file.'
                    continue
                 
                lens = self.getLens(os.path.join(contigs_folder, f))
                trimmedLens = self.trimLens(lens, min_length)
                
                if len(trimmedLens) == 0:
                    print f + ": no sequences longer than threshold"
                
                statN = len(lens)
                statTrimmedN = len(trimmedLens)
                statSum = sum(trimmedLens)
                statMax = max(trimmedLens)
                statMin = min(trimmedLens)
                statMed = trimmedLens[ (len(trimmedLens)-1) / 2 ]
                statMean = int( statSum / float(statTrimmedN) )
                statnxs = []
                statnx_lens = []
                
                for n in range(5,96,5):
                    statnx_len = self.nx(trimmedLens, n)
                    statnxs.append( str(n))
                    statnx_lens.append(statnx_len)
                    
                    if n == 50:
                        stat_n50 = statnx_len
                    if n == 90:
                        stat_n90 = statnx_len
                
                statnx_auc = sum(statnx_lens)
                
                # TODO Add analysis of ORF_lengths here
                for i in range(len(statnxs)):
                    nx_stats_out.write("\t".join(map(str, [sample, algo, statnxs[i], statnx_lens[i]])) + "\n")
                
                stats_out.write("\t".join(map(str,[sample, algo, statN, statTrimmedN, statSum, statMin, statMed, statMean, statMax, stat_n50, stat_n90, statnx_auc])) + "\n")
            
            # close files
            stats_out.close()
            nx_stats_out.close() 

    def trimLens(self, lens, minLen):
       """
       Eliminates any reads below a certain threshold.  Function assumes that input
       list lens is sorted smallest to largest.
       """
       index = 0
       for i in range(len(lens)):
           if lens[i] < minLen:
               index += 1
           else:
               break
    
       return lens[index:len(lens)]


    def assemble_samples(self, samples, mp2_parsed_params, algos=None, paired=False):

        #TODO loop over assembly types
        print "In assemble_samples"
        algos = ["spades", "idba-ud"]
        for sample in samples:
            self.assemble_sample(sample, samples[sample], algos, paired=True)
    
    def select_winning_assembly(self, sample):
        # output
        root = os.path.dirname(os.path.realpath(__file__))
        contigs_folder = os.path.join(root, "assemblies", sample, "assembly", "final_contigs") 
        stats_file = open(os.path.join(contigs_folder, "assembly_stats.txt"),"r")
        lines = stats_file.readlines()
        headers = lines[0]

        lengths = []
        aucs = []
        n50s = []
        algorithms = []
        for line in lines[1:]:
            fields = line.split("\t")
            sample = fields[0]
            algorithm = fields[1]
            n = fields[2]
            n_trimmed = fields[3]
            length = fields[4]
            stat_min = fields[5]
            stat_median = fields[6]
            stat_mean = fields[7]
            stat_max = fields[8]
            stat_n50 = fields[9]
            stat_n90 = fields[10]
            stat_auc = fields[11]
            
            lengths.append(length)
            aucs.append(stat_auc)
            n50s.append(stat_n50)
            algorithms.append(algorithm)
        
        length_winner = algorithms[lengths.index(max(lengths))]
        auc_winner = algorithms[aucs.index(max(aucs))]
        n50_winner = algorithms[n50s.index(max(n50s))]

        vote_list = [length_winner, auc_winner, n50_winner]
        max_votes = -1
        max_candidate = None
        for candidate in set(vote_list):
            votes = vote_list.count(candidate)
            if votes > max_votes:
               max_candidate = candidate
               max_votes = votes

        return max_candidate

    def runApp( self, downloaded_files, app_session = None ) :
        upload_folders = []
        mp2_parsed_params = self.parseMp2Input(app_session)
        paired_samples = self.getSequencePairs(downloaded_files)
# self.assemble_samples(paired_samples, mp2_parsed_params, paired=True)
        print "Calculating statistics"
        app_root = os.path.dirname(os.path.realpath(__file__))
        self.compute_assembly_statistics(paired_samples)
        mp_input = os.path.join(app_root, "mp_input")
        mp_input_assembly = None
        for sample in paired_samples:
            best_assembly = self.select_winning_assembly(sample)
            paired_samples[sample]["best_assembly"] = best_assembly
            
            if not os.path.exists(mp_input):
               os.mkdir(mp_input)
            
            contigs_file = os.path.join(app_root, "assemblies", sample, "assembly", "final_contigs", best_assembly + "_contigs.fasta")
            print contigs_file
            mp_input_assembly = os.path.join(mp_input, sample + "_" + best_assembly + "_contigs.fasta")
            shutil.copyfile(contigs_file, os.path.join(mp_input, sample + "_" + best_assembly + "_contigs.fasta"))
        
        if not mp_input_assembly:
            print "mp_input_assembly is None... Bail!"
            exit()

        self.createMp2ParameterFile(mp2_parsed_params)
         
         
        output_dir = os.path.join(app_root, "mp_output")
        config = os.path.join(app_root, "config")
        template_config = os.path.join(config, "template_config.txt")
        template_param = os.path.join(config, "template_param.txt")
        
        # create output folder
        appresults_dir = "/data/output/appresults/"
        if not os.path.exists(appresults_dir):
            os.makedirs(appresults_dir)
        
        project_id_dir = os.path.join(appresults_dir, self.project_id)
        for sample in paired_samples:
            sample_name = sample.replace(" ","_")
            sample_dir = os.path.join( project_id_dir, sample_name ) 
            r_output_dir = os.path.join(sample_dir , "r_output")
            
            if not os.path.exists(r_output_dir):
                os.makedirs(r_output_dir)
            final_contigs_dir = os.path.join(app_root, "assemblies", sample_name, "/assembly/final_contigs/")
            final_contigs_dir = "/".join([app_root, "assemblies", sample_name, "/assembly/final_contigs/"])
            # compute nx graph
            r_command = ["Rscript", \
                        os.path.join(app_root, "lib/nx_plot.R"), \
                        os.path.join(final_contigs_dir, "assembly_stats_nx.txt"),
                        r_output_dir + "/nx_contigs.png" ]
            print r_command 
            res = subprocess.call(r_command)
            if res == 0:
                print "Rcommand Success!"
            
            
            # move contig stats table
            res = shutil.move(os.path.join(final_contigs_dir, "assembly_stats.txt"),\
                            os.path.join(r_output_dir, "assembly_stats.tsv"))

            
        self.createSimpleMp2Command(mp_input, output_dir, template_config, template_param)
        mp_sample = os.path.basename(mp_input_assembly).replace(".fasta","")

        file_list = os.listdir(os.path.join(output_dir, mp_sample, "run_statistics"))
        orf_lengths_file = None
        contig_lengths_file = None
        print file_list
        for file in file_list:
            if re.match(".*\.orf\.lengths\.txt", file):
                orf_lengths_file = file
            elif re.match(".*\.contig\.lengths\.txt", file):
                contig_lengths_file = file
    
        print orf_lengths_file
        print contig_lengths_file
        # run new R command
        r_command = ["Rscript", \
                        os.path.join(app_root, "lib/histogram_plot.R"), \
                        os.path.join(output_dir, mp_sample,  "run_statistics", contig_lengths_file),\
                        r_output_dir + "/contig_lengths.png" ]
        res = subprocess.call(r_command)
        
        r_command = ["Rscript", \
                        os.path.join(app_root, "lib/histogram_plot.R"), \
                        os.path.join(output_dir, mp_sample,  "run_statistics", orf_lengths_file),\
                        r_output_dir + "/orf_lengths.png" ]
        res = subprocess.call(r_command)
        
          
        # print "Calling MetaPathways!"
        # self.createSimpleMp2Command(mp_input, output_dir, template_config, template_param)
         
        







        exit()


        for file in downloaded_files :
            if not (file.endswith( ".fastq" ) or \
                    file.endswith( ".fastq.gz") or \
                    file.endswith(".fas") ) : continue



            output_dir = os.path.join(os.path.dirname(me), "output")
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

