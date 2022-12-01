import java.io.*;
import java.nio.file.FileSystems;


public class Main_Launch_Exp {

    // Flag for distinguish exact from approximate solution.
    private static int Approx_Solution = 0;


    private static final int KB_ALLOC = 1024;
    private static final int size_Block_ESS_Dataset = 128 * KB_ALLOC * KB_ALLOC; // 134.217.728 bytes --> (128 MB)   // 268.435.456 bytes --> (256 MB)    // 536.870.912 bytes --> (512 MB)
    private static final int size_Block_DSK_Dictionary = ((2 * KB_ALLOC * KB_ALLOC * KB_ALLOC) - 1); // 2.147.483.647 bytes --> (2GB)

    private static final int MAX_KMER_ABUNDANCE = 2147483647;
    private static final int MIN_KMER_ABUNDANCE = 1;

    private static final String VERBATIM_Suffix = "-verbatim";
    private static final String LEXICOGR_Suffix = "-lexicogr";
    private static final String FREQUENCY_Suffix = "-frequency";
    private static final String BIJECTION_Suffix = "-bijection";

    private static final String FASTA_EXT = ".fasta";
    private static final String FASTQ_EXT = ".fastq";

    private static final String DSK_EXT = ".txt";
    private static final String KMC3_EXT = ".txt";
    private static final String KMC3_PREFIX_DB = ".kmc_pre";
    private static final String KMC3_SUFFIX_DB = ".kmc_suf";

    private static final String nfo_EXT = ".nfo";
    private static final String ESS_COMPR_EXT = ".essc";
    private static final String ESS_DECOMPR_EXT = ".essd";

    // compressor method extensions
    private static final String bzip2_c_EXT = ".bz2";
    private static final String lz4_c_EXT = ".lz4";
    private static final String zstd_c_EXT = ".zst";
    private static final String mfc_c_EXT = ".mfc";

    private static final String spring_c_EXT = ".spr";
    private static final String spring_d_EXT = ".dec";
    private static final String BIC_c_EXT = ".bic";
    private static final String OptPFOR_c_EXT = ".opt";
    private static final String FMind_c_EXT = ".fm9";

    // decompressor method extensions
    private static final String mfc_d_EXT = ".d";
    private static final String BIC_d_EXT = ".d";
    private static final String OptPFOR_d_EXT = ".d";
    private static final String FMind_d_EXT = ".dec";

    //BCSF extension
    private static final String BCSF_backup = "-Backup";
    private static final String BCSF_CSF_EXT = ".csf.bin";
    private static final String BCSF_BLOOM_EXT = ".0.bloom.txt";


    // DP size and time Output Files
    private static final String DP_EXT = ".txt";

    // file System File Separator
    private static final String FS_SEPARATOR = FileSystems.getDefault().getSeparator();

    // k values used in this research
    private static final int[] k_values = {
            4,
            8,
            16,
            32,
            48,
            64,
    };


    public static void main(String[] args) throws IOException, InterruptedException {

        //Dataset Directory, example: /home/nome-utente/TradeOff/data
        String dataset_directory = args[0];

        //Dataset Name, example: StaphAU
        String exp = args[1];

        String exp_path = args[2] + FS_SEPARATOR;

        //DEBUG flag 0/1: O for no verbose output; 1 for verbose output
        int DEBUG = Integer.parseInt(args[3]);


        // create folder experiment
        create_folder(exp_path, exp, DEBUG);

        System.out.println("Create experiment: " + exp);

        // absolute input dataset pathname (FASTA format)
        // for example: /home/nome-utente/IdeaProjects/TradeOff/data/StaphAU.fasta
        String dataset_pathname = dataset_directory + FS_SEPARATOR + exp + FASTA_EXT;

        // absolute input dataset pathname (FASTA format)
        // for example: /home/nome-utente/IdeaProjects/TradeOff/data/StaphAU.fasta
        String dataset_pathname_FASTQ = dataset_directory + FS_SEPARATOR + exp + FASTQ_EXT;

        // Experiment Folder
        // for example: /home/nome-utente/IdeaProjects/TradeOff/exp-automazione/StaphAU
        String exp_folder = exp_path + exp;


        // FIXME: Base_Case scenario
        System.out.println("Applying Base Case scenario on entire dataset ...");
        Base_Case(dataset_pathname, exp_folder, exp, DEBUG);

        // FIXME: Clean Garbage Collector
        System.gc();

        // FIXME: DP scenarios
        for (int k : k_values) {

            int b = Compute_Items_For_DSK_Fragment(k);
            System.out.println("Blocks of 2^b items for each DSK fragments, with  b = " + b);

            System.out.println("Set k value to: " + k + " !!!");

            // FIXME: create specific k folder into experiment folder experiment
            create_folder(exp_folder, "k" + k, DEBUG);
            System.out.println("Folder " + exp_folder + " created.");

            // FIXME: Extract (Dk,Fk) dictionary with DSK from dataset
            System.out.println("Extract via DSK a (Dk, Fk) dictionary from dataset ...");
            long pre_proc_DSK_start = System.currentTimeMillis();
            extract_DSK_dictionary(dataset_pathname, k, exp_folder, exp, DEBUG);
            long pre_proc_DSK_time = System.currentTimeMillis() - pre_proc_DSK_start;

            // TODO: Extract (Dk,Fk) dictionary with KMC from dataset
            long pre_proc_KMC3_start = System.currentTimeMillis();
            extract_KMC3_dictionary(dataset_pathname_FASTQ, k, exp_folder, exp, DEBUG);
            long pre_proc_KMC3_time = System.currentTimeMillis() - pre_proc_KMC3_start;
            System.out.println("pre_proc_KMC3_time = " + pre_proc_KMC3_time);               //TODO: save info into file

            // TODO: Canonize KMC Dictionary in according to DSK_Rel_Ord
            long pre_proc_KMC3_canonize_start = System.currentTimeMillis();
            DSKCanonize_KMC3_dictionary(exp_folder, exp, k, DEBUG);
            long pre_proc_KMC3_canonize_time = System.currentTimeMillis() - pre_proc_KMC3_canonize_start;
            System.out.println("pre_proc_KMC3_canonize_time = " + pre_proc_KMC3_canonize_time);               //TODO: save info into file

            //TODO: Transform [Dk Fk] into Dk and Fk separated files
            long pre_proc_KMC3_canonize_DkFk_start = System.currentTimeMillis();
            Transform_DkFk_KMC3_dictionary(exp_folder, exp, k, DEBUG);
            long pre_proc_KMC3_canonize_DkFk_time = System.currentTimeMillis() - pre_proc_KMC3_canonize_DkFk_start;
            System.out.println("pre_proc_KMC3_canonize_DkFk_time = " + pre_proc_KMC3_canonize_DkFk_time);               //TODO: save info into file


            // FIXME: Clean Garbage Collector
            System.gc();

            // FIXME: create scenario DP0 and put results into DP0_file
            System.out.println("Create DP0 scenario");
            DP_scenario_012(exp_folder, exp, k, pre_proc_DSK_time, 0, b, DEBUG);
            System.out.println("Create Histogram from DP0 scenario");
            Histogram_exp.execute(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR, exp + "-DSK-k" + k + VERBATIM_Suffix, DEBUG);

            // FIXME: Clean Garbage Collector
            System.gc();

            // FIXME: create scenario DP3 and put results into DP3_file
            System.out.println("Create DP3 scenario");

            // FIXME: Extract Spectrum Preserving String Set S from dataset
            System.out.println("Extract via ESS a Spectrum Preserving String Set S from dataset ...");
            long pre_proc_ESS_start = System.currentTimeMillis();
            extract_SPSS(dataset_pathname, k, exp_folder, exp, DEBUG);
            long pre_proc_ESS_time = System.currentTimeMillis() - pre_proc_ESS_start;

            DP_scenario_3(exp_folder, exp, 3, k, pre_proc_DSK_time, pre_proc_ESS_time, b, DEBUG);

            // TODO: Remove Dk, Fk, Gk fragments files from DP3 and DP0 scenario
            System.out.println("Remove Temporary DP3 and DP0 files");
            Remove_Temporary_files(exp_folder, exp, k, 3, DEBUG);
            Remove_Temporary_files(exp_folder, exp, k, 0, DEBUG);

            // FIXME: Clean Garbage Collector
            System.gc();

            // FIXME: create scenario DP1 and put results into DP1_file
            System.out.println("Create DP1 scenario");
            DP_scenario_012(exp_folder, exp, k, pre_proc_DSK_time, 1, b, DEBUG);

            // TODO: Remove Dk, Fk and Gk fragments files from DP1 scenario
            System.out.println("Remove Temporary DP1 files");
            Remove_Temporary_files(exp_folder, exp, k, 1, DEBUG);

            //FIXME: Clean Garbage Collector
            System.gc();

            //TODO: create scenario DP2 and put results into DP2_file
            System.out.println("Create DP2 scenario");
            DP_scenario_012(exp_folder, exp, k, pre_proc_DSK_time, 2, b, DEBUG);

            //TODO: Remove Dk, Fk and Gk fragments files from DP2 scenario
            System.out.println("Remove Temporary DP2 files");
            Remove_Temporary_files(exp_folder, exp, k, 2, DEBUG);

            //FIXME: Clean Garbage Collector
            System.gc();
        }

        // test jar usage: java -jar TradeOff.jar /home/ismartlab/IdeaProjects/TradeOff/data dataset /home/ismartlab/IdeaProjects/TradeOff/exp-automazione DEBUG[0/1]
    }

    private static int Compute_Items_For_DSK_Fragment(int k) {
        int b = 20;
        while (((Math.pow(2, b) - 1) * (k + 5)) <= (size_Block_DSK_Dictionary / 2)) {

            int aux = (int) (Math.pow(2, b) - 1) * (k + 5);

            // update b value for next comparison
            b++;
        }
        return b;
    }

    private static long Compute_size_Dataset(String dataset_pathname) {

        long size_D = 0;

        // compute size of fasta Dataset D
        File file = new File(dataset_pathname);
        if (file.exists() && file.isFile()) {
            size_D += file.length();
        } else {
            System.out.println(dataset_pathname + " not a file or non-existent file !!!.");
        }

        return size_D;
    }

    private static void watch(final Process process) {
        // FIXME: watch the process: This thread will finish in the calling process when the called process ends.

        new Thread() {
            public void run() {
                BufferedReader input = new BufferedReader(new InputStreamReader(process.getInputStream()));
                String line = null;
                try {
                    while ((line = input.readLine()) != null) {
                        System.out.println(line);
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }.start();
    }

    private static void create_folder(String path, String exp, int debug) throws InterruptedException, IOException {

        String cmd = "mkdir " + path + FS_SEPARATOR + exp;
        if (debug == 1) System.out.println(cmd);

        ProcessBuilder builder = new ProcessBuilder("mkdir", path + FS_SEPARATOR + exp);
        builder.redirectErrorStream(true);
        final Process process = builder.start();
        // Watch the process
        watch(process);
        process.waitFor();
    }

    private static void remove_folder(String path, int debug) throws InterruptedException, IOException {

        if (debug==1) System.out.println("rm -R " + path + FS_SEPARATOR);

        ProcessBuilder builder = new ProcessBuilder("rm", "-R", path + FS_SEPARATOR);
        builder.redirectErrorStream(true);
        final Process process = builder.start();
        // Watch the process
        watch(process);
        process.waitFor();
    }


    private static void remove_file(String path_file, int debug) throws InterruptedException, IOException {

        if (debug==1) System.out.println("rm -f " + path_file);

        ProcessBuilder builder = new ProcessBuilder("rm", "-f", path_file);
        builder.redirectErrorStream(true);
        final Process process = builder.start();
        // Watch the process
        watch(process);
        process.waitFor();
    }


    private static void move_file_to_path(String file, String path_in, String path_out, int debug) throws IOException, InterruptedException {

        if (debug==1) System.out.println("mv " + path_in + FS_SEPARATOR + file + " " + path_out + FS_SEPARATOR + file);

        ProcessBuilder builder = new ProcessBuilder("mv", path_in + FS_SEPARATOR + file, path_out + FS_SEPARATOR + file);
        builder.redirectErrorStream(true);
        final Process process = builder.start();
        // Watch the process
        watch(process);
        process.waitFor();
    }


    private static void extract_DSK_dictionary(String dataset_pathname, int k, String exp_folder, String exp_name, int debug) throws IOException, InterruptedException {

        // create .h5 file
        if (debug == 1) System.out.println("./essAuxDsk -file " + dataset_pathname + " -kmer-size " + k + " -abundance-min 1 -out " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-k" + k + ".h5" + " -verbose 0 -max-memory 5000");

        ProcessBuilder builder = new ProcessBuilder("./essAuxDsk", "-file", dataset_pathname, "-kmer-size", String.valueOf(k), "-abundance-min", String.valueOf(1), "-out", exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-k" + k + ".h5", "-verbose", String.valueOf(0), "-max-memory", String.valueOf(5000));
        builder.redirectErrorStream(true);
        final Process process = builder.start();
        // Watch the process
        watch(process);
        process.waitFor();

        // conversion .h5 file to .txt file
        if (debug ==1) System.out.println("./essAuxDsk2ascii -file " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-k" + k + ".h5" + " -out " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-DSK-k" + k + DSK_EXT + " -verbose 0");

        builder = new ProcessBuilder("./essAuxDsk2ascii", "-file", exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-k" + k + ".h5", "-out", exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-DSK-k" + k + DSK_EXT, "-verbose", String.valueOf(0));
        builder.redirectErrorStream(true);
        final Process process2 = builder.start();
        // Watch the process
        watch(process2);
        process2.waitFor();

        // remove .h5 file
        if (debug==1) System.out.println("rm -f " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-k" + k + ".h5");

        remove_file(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-k" + k + ".h5", debug);

        //Cleanup Garbage Collector
        System.gc();
    }

    private static void extract_KMC3_dictionary(String dataset_pathname_FASTQ, int k, String exp_folder, String exp_name, int debug) throws IOException, InterruptedException {

        if (debug==1) System.out.println("./kmc -v -ci" + MIN_KMER_ABUNDANCE + " " + "-cs" + MAX_KMER_ABUNDANCE + " -k" + k + " " + dataset_pathname_FASTQ + " " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + " " + exp_folder);

        ProcessBuilder builder = new ProcessBuilder("./kmc", "-v", "-ci" + MIN_KMER_ABUNDANCE, "-cs" + MAX_KMER_ABUNDANCE, "-k" + k, dataset_pathname_FASTQ, exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k, exp_folder);
        builder.redirectErrorStream(true);
        final Process process = builder.start();
        // Watch the process
        watch(process);
        process.waitFor();


        if (debug==1) System.out.println("./kmc_dump " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + " " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + KMC3_EXT);

        builder = new ProcessBuilder("./kmc_dump", exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k, exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + KMC3_EXT);
        builder.redirectErrorStream(true);
        final Process process2 = builder.start();
        // Watch the process
        watch(process2);
        process2.waitFor();


        // remove .kmc_pre file
        if (debug ==1) System.out.println("rm -f " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + KMC3_PREFIX_DB);
        remove_file(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + KMC3_PREFIX_DB, debug);


        // remove .kmc_suf file
        if (debug==1) System.out.println("rm -f " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + KMC3_SUFFIX_DB);
        remove_file(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + KMC3_SUFFIX_DB, debug);

        //Cleanup Garbage Collector
        System.gc();
    }

    private static void DSKCanonize_KMC3_dictionary(String exp_folder, String exp_name, int k, int debug) throws IOException, InterruptedException {

        // Input KMC3 Dictionary into <Kmer_can_lexicogr> <frequenza> format
        BufferedReader dict_in = new BufferedReader(new FileReader(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + KMC3_EXT));

        // Output KMC3 Dictionary file in Canonical DSK Relation Order into <Kmer_can_DSK_rel_ord> <frequenza> format
        PrintWriter dict_out  = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + "-DSKcan" + KMC3_EXT));

        int block_size = 50000000;
        if (debug == 1) System.out.println("Canonize KMC3 Dictionary. Input Data Loading...");

        String line;
        String[] aux;

        int count = 0;
        while ((line = dict_in.readLine()) != null) {

            aux = line.split("\\s+");

            //storage su files di output
            dict_out.println(StringUtils.Canonical_DSK_RelOrd_speedup(aux[0]) + " " + aux[1]);
            count++;

            if (debug == 1) {
                if (count % block_size == 0) {
                    System.out.println("Canonize first " + count + " rows from input KMC3 dictionary file");
                }
            }
        }
        if (debug==1) System.out.println("Canonize KMC3 input Dictionary File Terminated !!!");

        /* Closing files */
        dict_out.close();
        dict_in.close();

        //Cleanup Garbage Collector
        System.gc();
    }

    private static void Transform_DkFk_KMC3_dictionary(String exp_folder, String exp_name, int k, int debug) throws IOException, InterruptedException {

        // Input D in <Kmer> <frequency> format
        BufferedReader dict_in = new BufferedReader(new FileReader(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + "-DSKcan" + KMC3_EXT));

        // Output Kmers file in fasta format
        PrintWriter kmers_out  = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + "-DSKcan_Dk" + FASTA_EXT));

        // Output Frequency file in <frequency> \n format
        PrintWriter frequency_out  = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-KMC-k" + k + "-DSKcan_Fk" + DSK_EXT));

        int block_size = 50000000;

        if (debug == 1) System.out.println("Transform KMC3 Dictionary into Dk and Fk. Input Data Loading...");

        String line;
        String[] aux;

        // write the character '>' in the first line of the k-mer output file
        kmers_out.println(">");

        int count = 0;
        while ((line = dict_in.readLine()) != null) {

            aux = line.split("\\s+");

            kmers_out.println(aux[0]);
            frequency_out.println(aux[1]);

            count++;
            if (debug == 1) {
                if (count % block_size == 0) {
                    System.out.println("Save first " + count + " rows from input file");
                }
            }
        }

        if (debug==1) System.out.println("Created K-mer and Frequencies files !!!");

        /* Closing files */
        kmers_out.close();
        frequency_out.close();
        dict_in.close();
    }

    private static void extract_SPSS_single_D(String dataset_pathname, int k, String exp_folder, String exp_name, int debug) throws IOException, InterruptedException {
        // create .essc file
        if (debug==1) System.out.println("./essCompress -k " + k + " -i " + dataset_pathname + " -o " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR);

        ProcessBuilder builder = new ProcessBuilder("./essCompress", "-k", String.valueOf(k), "-i", dataset_pathname, "-o", exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR);
        builder.redirectErrorStream(true);
        final Process process = builder.start();
        // Watch the process
        watch(process);
        process.waitFor();

        // decompression .essc file to .essd decompressed textual file .fasta.essd
        if (debug==1) System.out.println("./essDecompress " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + FASTA_EXT + ESS_COMPR_EXT);

        builder = new ProcessBuilder("./essDecompress", exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + FASTA_EXT + ESS_COMPR_EXT);
        builder.redirectErrorStream(true);
        final Process process2 = builder.start();
        // Watch the process
        watch(process2);
        process2.waitFor();

        // renaming .fasta.essd file to fasta ESS output file
        if (debug==1) System.out.println("mv " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + FASTA_EXT + ESS_DECOMPR_EXT + " " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-ESS-k" + k + FASTA_EXT);

        builder = new ProcessBuilder("mv", exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + FASTA_EXT + ESS_DECOMPR_EXT, exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-ESS-k" + k + FASTA_EXT);
        builder.redirectErrorStream(true);
        final Process process3 = builder.start();
        // Watch the process
        watch(process3);
        process3.waitFor();

        // remove .essc file
        if (debug==1) System.out.println("rm -f " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + FASTA_EXT + ESS_COMPR_EXT);
        remove_file(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + FASTA_EXT + ESS_COMPR_EXT, debug);
    }

    private static void extract_SPSS(String dataset_pathname, int k, String exp_folder, String exp_name, int debug) throws IOException, InterruptedException {

        //Check Dataset dimension
        long size_D = Compute_size_Dataset(dataset_pathname);
        if (debug==1) System.out.println("Dataset D size [bytes] = " + size_D);

        // setup Dataset Block .nfo file for ESS compute
        int num_Blocks = 1;     //default value
        PrintWriter Blocks_info_out = null;
        try {
            Blocks_info_out = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-ESS-k" + k + nfo_EXT));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Check condition for divide genomic dataset
        // If size(Dataset) > size_Block_ESS_Dataset ==> I divide the dataset into blocks of approximately size_Block bytes and compute ESS on the various blocks
        if ( (size_D > (size_Block_ESS_Dataset)) && (k > 8) ) {    //We are in Approximate Solution Context

            // In k value set {4,8,16,32,48,64}, when k>13 in commodity hardware ESS failed, so we consider the context of an approximate solution

            System.out.println("We are in Approximate solution Context ... ");
            System.out.println("Please Note: Possible loss of k-mers.");
            Approx_Solution = 1;

            // compute Blocks number
            double aux_div = size_D / size_Block_ESS_Dataset;
            double aux_rest = size_D % size_Block_ESS_Dataset;
            num_Blocks = (int) ((aux_rest == 0) ? Math.floor(aux_div) : Math.floor(aux_div) + 1);
            System.out.println("num_Blocks = " + num_Blocks);

            // storage num_Blocks information
            Blocks_info_out.println(num_Blocks);
            Blocks_info_out.close();

            Divide_Dataset_into_Blocks(dataset_pathname, k, exp_folder, exp_name, debug);

            // TODO: Apply ESS to every single Block
            Extract_ESS_SPSS(exp_folder, exp_name, k, num_Blocks, debug);

        } else {    //We are in Exact Solution Context
            // FIXME: Non occorre dividere in più blocchi, per cui procedo in maniera classica
            //  bisogna comunque memorizzare l'informazione da qualche parte, esempio un file .nfo
            //  vedi (exp_name)-ESS-k().nfo

            // storage num_Blocks information
            Blocks_info_out.println(num_Blocks);
            Blocks_info_out.close();

            // TODO: computa ESS sull'intero Dataset genomico (old-style ma con nomenclatura nuova)

            //We are in Exact Solution Context
            Approx_Solution = 0;
            System.out.println("We are in Exact Solution Context ... ");
            extract_SPSS_single_D(dataset_pathname, k, exp_folder, exp_name, debug);
        }


        // TODO: concatena i blocchi e rinomina come nel caso a singolo blocco
        String block = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-ESS-k" + k + "_Block_";
        String ESS_output = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-ESS-k" + k + FASTA_EXT;
        Concatenate_ESS_Blocks (block, num_Blocks, ESS_output, debug);
    }

    private static void Extract_ESS_SPSS(String exp_folder, String exp_name, int k, int num_blocks, int debug) throws IOException, InterruptedException {

        for (int i = 0; i < num_blocks; i++) {       //0 .. <num_blocks

            //FIXME: apply ESS software on each block
            String block = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-ESS-k" + k + "_Block_" + i;

            // FIXME: compression .essc file
            if (debug==1) System.out.println("./essCompress -k " + k + " -i " + block + FASTA_EXT + " -o " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR);

            ProcessBuilder builder = new ProcessBuilder("./essCompress", "-k", String.valueOf(k), "-i", block + FASTA_EXT, "-o", exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR);
            builder.redirectErrorStream(true);
            final Process process = builder.start();
            // Watch the process
            watch(process);
            process.waitFor();

            // FIXME: decompression .essc file to .essd decompressed textual file .fasta.essd
            if (debug==1) System.out.println("./essDecompress " + block + FASTA_EXT + ESS_COMPR_EXT);

            builder = new ProcessBuilder("./essDecompress", block + FASTA_EXT + ESS_COMPR_EXT);
            builder.redirectErrorStream(true);
            Process process2 = builder.start();
            // Watch the process
            watch(process2);
            process2.waitFor();


            // FIXME: remove .essc file
            if (debug==1) System.out.println("rm -f " + block + FASTA_EXT + ESS_COMPR_EXT);
            remove_file(block + FASTA_EXT + ESS_COMPR_EXT, debug);


            // FIXME: copy .fasta genomic block of dataset D
            if (debug==1) System.out.println("cp " + block + FASTA_EXT + " " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "_B_" + i + FASTA_EXT);

            builder = new ProcessBuilder("cp", block + FASTA_EXT, exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "_B_" + i + FASTA_EXT);
            builder.redirectErrorStream(true);
            process2 = builder.start();
            // Watch the process
            watch(process2);
            process2.waitFor();


            // FIXME: move .fasta.essd decompressed file into .fasta file
            if (debug==1) System.out.println("mv -f " + block + FASTA_EXT + ESS_DECOMPR_EXT + " " + block + FASTA_EXT);

            builder = new ProcessBuilder("mv", "-f", block + FASTA_EXT + ESS_DECOMPR_EXT, block + FASTA_EXT);
            builder.redirectErrorStream(true);
            process2 = builder.start();
            // Watch the process
            watch(process2);

            // Clean Garbace Collector
            System.gc();
        }
    }

    private static void Concatenate_ESS_Blocks(String block, int num_blocks, String ess_output, int debug) {

        PrintWriter S_out = null;
        try {
            S_out = new PrintWriter(new FileWriter(ess_output));

            for (int i = 0; i < num_blocks; i++) {

                String l;

                BufferedReader S_i = null;
                try {
                    S_i = new BufferedReader(new FileReader(block + i + FASTA_EXT));
                    if (debug==1) System.out.println("ESS Block: " + block + i + FASTA_EXT);

                    while ((l = S_i.readLine()) != null) {
                        S_out.println(l);
                    }

                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                }

                // close current ESS input block
                S_i.close();
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // close ESS output file
        S_out.close();

    }

    private static void Divide_Dataset_into_Blocks(String dataset_pathname, int k, String exp_folder, String exp_name, int debug) {

        // Input from Dataset
        String string_D = null;

        //String last_row_current e first_row_next
        String last_prev = null;
        String first_next = null;

        //String last_row_read = null;
        long bytes_readed = 0;
        int count_block = 0;

        BufferedReader Dataset_in = null;
        try {
            Dataset_in = new BufferedReader(new FileReader(dataset_pathname));

            //FIXME: write Block B_i in fasta format
            PrintWriter Block_i_out = null;
            try {
                Block_i_out = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-ESS-k" + k + "_Block_" + count_block + FASTA_EXT));
                if (debug==1) System.out.println("Extract Block " + count_block + " from Genomic Dataset D");
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            //TODO: write Block B_residues for recovery k-mer lost.
            PrintWriter Block_residues_out = null;
            try {
                Block_residues_out = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-ESS-k" + k + "_Block_" + "residues" + FASTA_EXT));
                if (debug==1) System.out.println("Prepare Extract Block " + count_block + " from Genomic Dataset");
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            while (((string_D = Dataset_in.readLine()) != null)) {
                //FIXME: read a line from Dataset D

                if (bytes_readed < size_Block_ESS_Dataset) {   // If I am in the range of block size allowed

                    // save data into current block
                    Block_i_out.println(string_D);

                    //update bytes readed counter
                    bytes_readed += string_D.length() + 1;

                    // update last_row current block
                    last_prev = string_D;

                } else {  // If, on the other hand, I am not in the allowed block size range then I close the current block and create a new one

                    // Current Block Closure
                    Block_i_out.close();

                    // update block counter
                    count_block++;

                    // start new block
                    //TODO: scrittura su Blocco Bi in formato fasta
                    Block_i_out = null;
                    try {
                        Block_i_out = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-ESS-k" + k + "_Block_" + count_block + FASTA_EXT));
                        if (debug==1) System.out.println("Extract Block " + count_block + " from Genomic Dataset");
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }

                    // if the string is a header
                    if (string_D.startsWith(">")) {

                        // start new block with this data
                        Block_i_out.println(string_D);

                        // set bytes readed to last line acquired
                        bytes_readed = string_D.length() + 1;

                    } else {

                        // start new block with header followed by last data
                        Block_i_out.println(">header"); //add 8 characters

                        // refactor with last_row_read data
                        Block_i_out.println(string_D);

                        // set bytes readed to last line acquired plus header line characters
                        bytes_readed = string_D.length() + 9;  // add "\n" + length(header) + "\n" and remove addiction of last_row_read.length() characters

                        // update first_row current block
                        first_next = string_D;

                        // Save last_prev and first_next to Block residues
                        int label_prec_block = count_block - 1;

                        // Save boundary reads into residues block file
                        Block_residues_out.println(">" + label_prec_block + "-" + count_block);

                        //TODO:  instead of storing all the previous read and the next read, I store the last (k-1) nucleotides of last_prev,
                        // combined with the first (k-1) nucleotides of first_next, which I will feed to DSK, when I need to extract the
                        // dictionary of missing k-meres or for which to cumulate the respective frequency ...

                        //Block_residues_out.println(last_prev);
                        //Block_residues_out.println(first_next);
                        String aux_read = last_prev.substring(last_prev.length() - k + 1, last_prev.length()) + first_next.substring(0, k);
                        Block_residues_out.println(aux_read);
                    }
                }

                // last_row_read update
                //last_row_read = string_D;
            }

            // close last block
            Block_i_out.close();

            //close Block of residues
            Block_residues_out.close();

            // close input dataset
            Dataset_in.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException ioException) {
            ioException.printStackTrace();
        }

        // Cleanup Garbage Collector
        System.gc();
    }

    private static void Remove_Temporary_files(String exp_folder, String exp_name, int k, int dp_scenario, int debug) throws IOException, InterruptedException {

        String mode = "";
        switch (dp_scenario) {
            case 0: {
                mode = VERBATIM_Suffix;
                break;
            }

            case 1: {
                mode = LEXICOGR_Suffix;
                break;
            }

            case 2: {
                mode = FREQUENCY_Suffix;
                break;
            }

            case 3: {
                mode = BIJECTION_Suffix;
                break;
            }
        }

        // remove path folder (es. /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-automazione/StaphAU/k4/StaphAU-DSK-k4-frequency/ )
        String path = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-DSK-k" + k + mode;

        remove_folder(path, debug);

        // remove "sorted" DSK output file (es. /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-automazione/StaphAU/k4/StaphAU-DSK-k4-frequency.txt )
        String file = path + DSK_EXT;
        remove_file(file, debug);
    }


    private static void Base_Case(String dataset_pathname, String exp_folder, String exp, int debug) throws IOException, InterruptedException {
        // create DP0 folder and put following file
        String DP_path = exp_folder + FS_SEPARATOR;
        String DP_folder = "Base";
        create_folder(DP_path, DP_folder, debug);

        // Create and Open DP files depends on: exp_folder, exp_name, k, DP_scenario
        PrintWriter DPBase_size_out = null;
        try {
            DPBase_size_out = new PrintWriter(new FileWriter(DP_path + FS_SEPARATOR + DP_folder + FS_SEPARATOR + "DP_" + DP_folder + "-size" + DP_EXT));
            if (debug == 1) System.out.println("Write DP file to: " + DP_path + FS_SEPARATOR + DP_folder + FS_SEPARATOR + "DP_" + DP_folder + "-size" + DP_EXT);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        PrintWriter DPBase_time_out = null;
        try {
            DPBase_time_out = new PrintWriter(new FileWriter(DP_path + FS_SEPARATOR + DP_folder + FS_SEPARATOR + "DP_" + DP_folder + "-time" + DP_EXT));
            if (debug == 1) System.out.println("Write DP file to: " + DP_path + FS_SEPARATOR + DP_folder + FS_SEPARATOR + "DP_" + DP_folder + "-time" + DP_EXT);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }


        // CD-NRAM Base Case Scenario
        Base_Case_CD_NRAM(dataset_pathname, exp_folder, exp, DP_folder, debug, DPBase_size_out, DPBase_time_out);

        // SD-RAM Base Case Scenario
        Base_Case_SD_RAM(dataset_pathname, exp_folder, exp, DP_folder, debug, DPBase_size_out, DPBase_time_out);


        // close DPBase files
        DPBase_time_out.close();
        DPBase_size_out.close();
    }


    private static void Base_Case_SD_RAM(String dataset_pathname, String exp_folder, String exp, String DP_folder, int debug, PrintWriter dpBase_size_out, PrintWriter dpBase_time_out) throws IOException, InterruptedException {

        create_folder(exp_folder + FS_SEPARATOR + DP_folder, "SD-RAM", debug);

        String cmd;
        long uncompr_size_D = 0;

        // compute size of fasta Dataset D
        File file = new File(dataset_pathname);
        if (file.exists() && file.isFile()) {
            uncompr_size_D += file.length();

            // TODO: Compressione via FM-index (Build)
            long FM_size_D = 0;
            long FM_time_compr_D;

            long FM_start_compr_D = System.currentTimeMillis();

            if (debug == 1) System.out.println("./fm-build " + dataset_pathname);

            ProcessBuilder builderC = new ProcessBuilder("./fm-build", dataset_pathname);
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update FM_size [bytes] with size of fasta dataset D
            file = new File(dataset_pathname + FMind_c_EXT);
            if (file.exists() && file.isFile()) {
                FM_size_D += file.length();
            } else {
                System.out.println(dataset_pathname + FMind_c_EXT + " not a file or non-existent file !!!.");
            }
            FM_time_compr_D = System.currentTimeMillis() - FM_start_compr_D;


            // TODO: Decompressione via FM-index (Extract)
            long FM_time_decompr_D;

            long FM_start_decompr_FM_of_D = System.currentTimeMillis();

            if (debug == 1) System.out.println("./fm-extract " + dataset_pathname + FMind_c_EXT + " " + dataset_pathname);

            //memo locate: exp_folder + FS_SEPARATOR + DP_folder + FS_SEPARATOR + "SD-RAM" + FS_SEPARATOR + exp + FASTA_EXT
            ProcessBuilder builderD = new ProcessBuilder("./fm-extract", dataset_pathname + FMind_c_EXT, dataset_pathname);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
            FM_time_decompr_D = System.currentTimeMillis() - FM_start_decompr_FM_of_D;

            //FIXME: il file .fm9 al momento non lo cancello, mi può servire per altre cose. NOTA: rimuovo solo il .dec che è pari al Dataset D in formato .fasta
            remove_file(dataset_pathname + FMind_d_EXT, debug);

            //FIXME: sposto il .fm9 nella cartella SD-RAM
            ProcessBuilder builderM = new ProcessBuilder("mv", dataset_pathname + FMind_c_EXT, exp_folder + FS_SEPARATOR + DP_folder + FS_SEPARATOR + "SD-RAM" + FS_SEPARATOR);
            builderM.redirectErrorStream(true);
            final Process processM = builderM.start();
            // Watch the process
            watch(processM);
            processM.waitFor();

            //FIXME: write to DPBase file size and time information !!!.
            //size
            dpBase_size_out.println(uncompr_size_D + ", " + FM_size_D);
            //time
            dpBase_time_out.println(FM_time_compr_D);
            dpBase_time_out.println(FM_time_decompr_D);

        } else {
            System.out.println(dataset_pathname + " not a file or non-existent file !!!.");
        }
    }


    private static void Base_Case_CD_NRAM(String dataset_pathname, String exp_folder, String exp, String DP_folder, int debug, PrintWriter dpBase_size_out, PrintWriter dpBase_time_out) throws IOException, InterruptedException {

        create_folder(exp_folder + FS_SEPARATOR + DP_folder, "CD-NRAM", debug);

        String cmd;
        long uncompr_size_D = 0;

        // compute size of fasta Dataset D
        File file = new File(dataset_pathname);
        if (file.exists() && file.isFile()) {
            uncompr_size_D += file.length();

            // bzip2 compression and decompression

            //FIXME: bzip2 compression
            long bzip2_start_compr = System.currentTimeMillis();
            long bzip2_size = 0;

            if (debug == 1) System.out.println("bzip2 -z -k " + dataset_pathname + " --best");

            ProcessBuilder builderC = new ProcessBuilder("bzip2", "-v", "-z", "-k", dataset_pathname, "--best");
            builderC.redirectErrorStream(true);
            Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update bzip2_size [bytes] with size of current fasta fragment
            file = new File(dataset_pathname + bzip2_c_EXT);
            if (file.exists() && file.isFile()) {
                bzip2_size += file.length();
            } else {
                System.out.println(dataset_pathname + bzip2_c_EXT + " not a file or non-existent file !!!.");
            }
            long bzip2_time_compr = System.currentTimeMillis() - bzip2_start_compr;

            //FIXME: bzip2 decompression
            long bzip2_start_decompr = System.currentTimeMillis();

            if (debug == 1) System.out.println("bzip2 -d -k -f " + dataset_pathname + bzip2_c_EXT);

            ProcessBuilder builderD = new ProcessBuilder("bzip2", "-v", "-d", "-k", "-f", dataset_pathname + bzip2_c_EXT);
            builderD.redirectErrorStream(true);
            Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
            long bzip2_time_decompr = System.currentTimeMillis() - bzip2_start_decompr;

            //FIXME: remove bzip2 compression files
            remove_file(dataset_pathname + bzip2_c_EXT, debug);


            // lz4 compression and decompression
            //FIXME: lz4 decompression
            long lz4_start_compr = System.currentTimeMillis();
            long lz4_size = 0;

            if (debug == 1) System.out.println("lz4 -z -f -k -9 " + dataset_pathname + " " + dataset_pathname + lz4_c_EXT);

            builderC = new ProcessBuilder("lz4", "-z", "-f", "-k", "-9", dataset_pathname, dataset_pathname + lz4_c_EXT);
            builderC.redirectErrorStream(true);
            processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update lz4_size [bytes] with size of current fasta fragment
            file = new File(dataset_pathname + lz4_c_EXT);
            if (file.exists() && file.isFile()) {
                lz4_size += file.length();
            } else {
                System.out.println(dataset_pathname + lz4_c_EXT + " not a file or non-existent file !!!.");
            }
            long lz4_time_compr = System.currentTimeMillis() - lz4_start_compr;

            //FIXME: lz4 decompression
            long lz4_start_decompr = System.currentTimeMillis();

            if (debug == 1) System.out.println("lz4 -d -k " + dataset_pathname + lz4_c_EXT + " " + dataset_pathname);

            builderD = new ProcessBuilder("lz4", "-d", "-k", "-f", dataset_pathname + lz4_c_EXT, dataset_pathname);
            builderD.redirectErrorStream(true);
            processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
            long lz4_time_decompr = System.currentTimeMillis() - lz4_start_decompr;

            //FIXME: remove lz4 compression files
            remove_file(dataset_pathname + lz4_c_EXT, debug);


            // zstd compression and decompression
            //FIXME: zstd compression
            long zstd_start_compr = System.currentTimeMillis();
            long zstd_size = 0;

            if (debug == 1) System.out.println("zstd -19 " + dataset_pathname + " --format=zstd");

            builderC = new ProcessBuilder("zstd", "-19", dataset_pathname, "--format=zstd");
            builderC.redirectErrorStream(true);
            processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update zstd_size [bytes] with size of current fasta fragment
            file = new File(dataset_pathname + zstd_c_EXT);
            if (file.exists() && file.isFile()) {
                zstd_size += file.length();
            } else {
                System.out.println(dataset_pathname + zstd_c_EXT + " not a file or non-existent file !!!.");
            }
            long zstd_time_compr = System.currentTimeMillis() - zstd_start_compr;

            //FIXME: zstd decompression
            long zstd_start_decompr = System.currentTimeMillis();

            if (debug == 1) System.out.println("zstd -d -f " + dataset_pathname + zstd_c_EXT);

            builderD = new ProcessBuilder("zstd", "-d", "-f", dataset_pathname + zstd_c_EXT);
            builderD.redirectErrorStream(true);
            processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
            long zstd_time_decompr = System.currentTimeMillis() - zstd_start_decompr;

            //FIXME: remove zstd compression files
            remove_file(dataset_pathname + zstd_c_EXT, debug);


            // MFC compression and decompression
            //FIXME: MFC compression
            long mfc_start_compr = System.currentTimeMillis();
            long mfc_size = 0;

            if (debug == 1) System.out.println("./essAuxMFCompressC -3 -p 2 " + dataset_pathname);

            builderC = new ProcessBuilder("./essAuxMFCompressC", "-3", "p 2", "-v", dataset_pathname);
            builderC.redirectErrorStream(true);
            processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update mfc_size [bytes] with size of current fasta fragment
            file = new File(dataset_pathname + mfc_c_EXT);
            if (file.exists() && file.isFile()) {
                mfc_size += file.length();
            } else {
                System.out.println(dataset_pathname + mfc_c_EXT + " not a file or non-existent file !!!.");
            }
            long mfc_time_compr = System.currentTimeMillis() - mfc_start_compr;

            //FIXME: mfc decompression
            long mfc_start_decompr = System.currentTimeMillis();

            if (debug == 1) System.out.println("./essAuxMFCompressD " + dataset_pathname + mfc_c_EXT);

            builderD = new ProcessBuilder("./essAuxMFCompressD", "-v", dataset_pathname + mfc_c_EXT);
            builderD.redirectErrorStream(true);
            processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
            long mfc_time_decompr = System.currentTimeMillis() - mfc_start_decompr;

            //FIXME: remove mfc compression files
            remove_file(dataset_pathname + mfc_c_EXT, debug);
            remove_file(dataset_pathname + mfc_c_EXT + mfc_d_EXT, debug);


            // Spring compression and decompression
            //FIXME: Spring compression
            long spr_start_compr = System.currentTimeMillis();
            long spr_size = 0;

            if (debug == 1) System.out.println("./spring -c (--fasta-input) -l -i " + dataset_pathname + " -o " + dataset_pathname + spring_c_EXT);

            builderC = new ProcessBuilder("./spring", "-c", "--fasta-input", "-l", "-i", dataset_pathname, "-o", dataset_pathname + spring_c_EXT);
            builderC.redirectErrorStream(true);
            processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update spr_size [bytes] with size of current fasta fragment
            file = new File(dataset_pathname + spring_c_EXT);
            if (file.exists() && file.isFile()) {
                spr_size += file.length();
            } else {
                System.out.println(dataset_pathname + spring_c_EXT + " not a file or non-existent file !!!.");
            }
            long spr_time_compr = System.currentTimeMillis() - spr_start_compr;

            //FIXME: spring decompression
            long spr_start_decompr = System.currentTimeMillis();

            if (debug == 1) System.out.println("./spring -d -i " + dataset_pathname + spring_c_EXT + " -o " + exp_folder + FS_SEPARATOR + DP_folder + FS_SEPARATOR + exp + "_decompr" + FASTA_EXT);

            builderD = new ProcessBuilder("./spring", "-d", "-i", dataset_pathname + spring_c_EXT, "-o", exp_folder + FS_SEPARATOR + DP_folder + FS_SEPARATOR + exp + "_decompr" + FASTA_EXT);
            builderD.redirectErrorStream(true);
            processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
            long spr_time_decompr = System.currentTimeMillis() - spr_start_decompr;

            //FIXME: remove spring compression files
            remove_file(dataset_pathname + spring_c_EXT, debug);
            remove_file(exp_folder + FS_SEPARATOR + DP_folder + FS_SEPARATOR + exp + "_decompr" + FASTA_EXT, debug);

            //FIXME: write to DPBase file size and time information !!!.
            //size
            dpBase_size_out.println(uncompr_size_D + ", " + bzip2_size + ", " + lz4_size + ", " + zstd_size + ", " + mfc_size + ", " + spr_size);
            //time compr
            dpBase_time_out.println(bzip2_time_compr + ", " + lz4_time_compr + ", " + zstd_time_compr + ", " + mfc_time_compr + ", " + spr_time_compr);
            //time decompr
            dpBase_time_out.println(bzip2_time_decompr + ", " + lz4_time_decompr + ", " + zstd_time_decompr + ", " + mfc_time_decompr + ", " + spr_time_decompr);

        } else {
            System.out.println(dataset_pathname + " not a file or non-existent file !!!.");
        }

    }


    private static void Manage_Dk_fragments(String path_in, String exp, int DP_scenario, PrintWriter DPfile_size, PrintWriter DPfile_time, int debug) throws IOException, InterruptedException {

        //FIXME: read .nfo file
        String l;
        int num_rows_file = 0;
        long total_rows = 0;
        int num_fragments = 0;
        int k = 0;

        // Reading .nfo information to derive information: total_rows, k, num_rows_file, num_fragments
        BufferedReader Header_NFO = null;
        try {
            Header_NFO = new BufferedReader(new FileReader(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + nfo_EXT));
            if (debug==1) System.out.println("Read .nfo info from: " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + nfo_EXT);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Extract info from .nfo file
        l = Header_NFO.readLine();
        total_rows = Long.parseLong(l);

        l = Header_NFO.readLine();
        k = Integer.parseInt(l);

        l = Header_NFO.readLine();
        num_rows_file = Integer.parseInt(l);

        l = Header_NFO.readLine();
        num_fragments = Integer.parseInt(l);

        // File Closure
        Header_NFO.close();

        //FIXME: Opening .fasta fragments for reading dimensions and cumulative sum of them.
        long uncompr_size_dk = 0;
        for (int i = 0; i < num_fragments; i++) {

            // compute size of fasta fragmentDict_size_DSK
            File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + i + FASTA_EXT);
            if (file.exists() && file.isFile()) {
                uncompr_size_dk += file.length();
            } else {
                System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + i + FASTA_EXT + " not a file or non-existent file !!!.");
            }
        }

        //FIXME: compress/uncompress phase with generic and specialized Dk compression methods
        String cmd;

        //FIXME: bzip2 compression
        long bzip2_start_compr = System.currentTimeMillis();
        long bzip2_size_dk = 0;
        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("bzip2 -z -k " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + " --best");

            ProcessBuilder builderC = new ProcessBuilder("bzip2", "-v", "-z", "-k", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT, "--best");
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update bzip2_size [bytes] with size of current fasta fragment
            File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + bzip2_c_EXT);
            if (file.exists() && file.isFile()) {
                bzip2_size_dk += file.length();
            } else {
                System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + bzip2_c_EXT + " not a file or non-existent file !!!.");
            }
        }
        long bzip2_time_compr_dk = System.currentTimeMillis() - bzip2_start_compr;

        //FIXME: bzip2 decompression
        long bzip2_start_decompr = System.currentTimeMillis();
        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("bzip2 -d -k -f " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + bzip2_c_EXT);

            ProcessBuilder builderD = new ProcessBuilder("bzip2", "-v", "-d", "-k", "-f", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + bzip2_c_EXT);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
        }
        long bzip2_time_decompr_dk = System.currentTimeMillis() - bzip2_start_decompr;

        //FIXME: remove bzip2 compression files
        for (int j = 0; j < num_fragments; j++) {
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + bzip2_c_EXT, debug);
        }


        //FIXME: lz4 compression
        long lz4_start_compr = System.currentTimeMillis();
        long lz4_size_dk = 0;
        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("lz4 -z -f -k -9 " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + " " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + lz4_c_EXT);

            ProcessBuilder builderC = new ProcessBuilder("lz4", "-z", "-f", "-k", "-9", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT, path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + lz4_c_EXT);
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update lz4_size [bytes] with size of current fasta fragment
            File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + lz4_c_EXT);
            if (file.exists() && file.isFile()) {
                lz4_size_dk += file.length();
            } else {
                System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + lz4_c_EXT + " not a file or non-existent file !!!.");
            }

        }
        long lz4_time_compr_dk = System.currentTimeMillis() - lz4_start_compr;

        //FIXME: lz4 decompression
        long lz4_start_decompr = System.currentTimeMillis();
        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("lz4 -d -k " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + lz4_c_EXT + " " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT);

            ProcessBuilder builderD = new ProcessBuilder("lz4", "-d", "-k", "-f", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + lz4_c_EXT, path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
        }
        long lz4_time_decompr_dk = System.currentTimeMillis() - lz4_start_decompr;

        //FIXME: remove lz4 compression files
        for (int j = 0; j < num_fragments; j++) {
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + lz4_c_EXT, debug);
        }


        //FIXME: zstd compression
        long zstd_start_compr = System.currentTimeMillis();
        long zstd_size_dk = 0;
        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("zstd -19 " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + " --format=zstd");

            ProcessBuilder builderC = new ProcessBuilder("zstd", "-19", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT, "--format=zstd");
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update zstd_size [bytes] with size of current fasta fragment
            File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + zstd_c_EXT);
            if (file.exists() && file.isFile()) {
                zstd_size_dk += file.length();
            } else {
                System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + zstd_c_EXT + " not a file or non-existent file !!!.");
            }
        }
        long zstd_time_compr_dk = System.currentTimeMillis() - zstd_start_compr;

        //FIXME: zstd decompression
        long zstd_start_decompr = System.currentTimeMillis();
        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("zstd -d -f " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + zstd_c_EXT);

            ProcessBuilder builderD = new ProcessBuilder("zstd", "-d", "-f", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + zstd_c_EXT);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
        }
        long zstd_time_decompr_dk = System.currentTimeMillis() - zstd_start_decompr;

        //FIXME: remove zstd compression files
        for (int j = 0; j < num_fragments; j++) {
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + zstd_c_EXT, debug);
        }


        //FIXME: MFC compression
        long mfc_start_compr = System.currentTimeMillis();
        long mfc_size_dk = 0;
        for (int j = 0; j < num_fragments; j++) {


            if (debug==1) System.out.println("./essAuxMFCompressC -3 -p 2 " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT);

            ProcessBuilder builderC = new ProcessBuilder("./essAuxMFCompressC", "-3", "p 2", "-v", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT);
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update mfc_size [bytes] with size of current fasta fragment
            File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + mfc_c_EXT);
            if (file.exists() && file.isFile()) {
                mfc_size_dk += file.length();
            } else {
                System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + mfc_c_EXT + " not a file or non-existent file !!!.");
            }
        }
        long mfc_time_compr_dk = System.currentTimeMillis() - mfc_start_compr;

        //FIXME: mfc decompression
        long mfc_start_decompr = System.currentTimeMillis();
        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("./essAuxMFCompressD " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + mfc_c_EXT);

            ProcessBuilder builderD = new ProcessBuilder("./essAuxMFCompressD", "-v", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + mfc_c_EXT);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
        }
        long mfc_time_decompr_dk = System.currentTimeMillis() - mfc_start_decompr;

        //FIXME: remove mfc compression files
        for (int j = 0; j < num_fragments; j++) {
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + mfc_c_EXT, debug);
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + mfc_c_EXT + mfc_d_EXT, debug);
        }


        //FIXME: Spring compression
        long spr_start_compr = System.currentTimeMillis();
        long spr_size_dk = 0;
        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("./spring -c (--fasta-input) (-l) -i " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + " -o " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + spring_c_EXT);

            ProcessBuilder builderC = new ProcessBuilder("./spring", "-c", "--fasta-input", "-l", "-i", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT, "-o", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + spring_c_EXT);
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update spr_size [bytes] with size of current fasta fragment
            File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + spring_c_EXT);
            if (file.exists() && file.isFile()) {
                spr_size_dk += file.length();
            } else {
                System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + spring_c_EXT + " not a file or non-existent file !!!.");
            }
        }
        long spr_time_compr_dk = System.currentTimeMillis() - spr_start_compr;

        //FIXME: spring decompression
        long spr_start_decompr = System.currentTimeMillis();
        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("./spring -d -i " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + spring_c_EXT + " -o " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_decompr" + "_" + j + FASTA_EXT);

            ProcessBuilder builderD = new ProcessBuilder("./spring", "-d", "-i", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + spring_c_EXT, "-o", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_decompr" + "_" + j + FASTA_EXT);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
        }
        long spr_time_decompr_dk = System.currentTimeMillis() - spr_start_decompr;

        //FIXME: remove spring compression files
        for (int j = 0; j < num_fragments; j++) {
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + spring_c_EXT, debug);
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_decompr" + "_" + j + FASTA_EXT, debug);
        }


        //TODO: FM-index
        long FM_size_dk = 0;
        long FM_time_compr_dk = 0;
        long FM_time_decompr_dk = 0;

        if (DP_scenario == 0) { // Build FM-index only for the DP0 Scenario (verbatim)

            //FIXME: FM-index(Dk) Compression (build FM-index)
            long FM_start_compr = System.currentTimeMillis();

            for (int j = 0; j < num_fragments; j++) {

                if (debug==1) System.out.println("./fm-build " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT);

                ProcessBuilder builderC = new ProcessBuilder("./fm-build", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT);
                builderC.redirectErrorStream(true);
                final Process processC = builderC.start();
                // Watch the process
                watch(processC);
                processC.waitFor();


                //update FM_size [bytes] with size of current fasta fragment
                File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + FMind_c_EXT);
                if (file.exists() && file.isFile()) {
                    FM_size_dk += file.length();
                } else {
                    System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + FMind_c_EXT + " not a file or non-existent file !!!.");
                }
            }
            FM_time_compr_dk = System.currentTimeMillis() - FM_start_compr;


            //FIXME: FM-index(Dk) Decompression  (.fasta.fm9 & .fasta -> .fasta.dec)  NOTE: (.fasta.dec === .fasta)
            long FM_start_decompr = System.currentTimeMillis();

            for (int j = 0; j < num_fragments; j++) {

                if (debug==1) System.out.println("./fm-extract " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + FMind_c_EXT + " " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT);

                ProcessBuilder builderD = new ProcessBuilder("./fm-extract", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + FMind_c_EXT, path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT);
                builderD.redirectErrorStream(true);
                final Process processD = builderD.start();
                // Watch the process
                watch(processD);
                processD.waitFor();
            }
            FM_time_decompr_dk = System.currentTimeMillis() - FM_start_decompr;

            //FIXME: remove spring compression/decompression  files (.fasta.dec)
            for (int j = 0; j < num_fragments; j++) {
                remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_" + j + FASTA_EXT + FMind_d_EXT, debug);
            }

        }

        //FIXME: write to DPfile size and time information !!!.
        if (DP_scenario == 0) {
            //size
            DPfile_size.println(uncompr_size_dk + ", " + bzip2_size_dk + ", " + lz4_size_dk + ", " + zstd_size_dk + ", " + mfc_size_dk + ", " + spr_size_dk + " - " + uncompr_size_dk + ", " + FM_size_dk);

            //compression time
            DPfile_time.println(bzip2_time_compr_dk + ", " + lz4_time_compr_dk + ", " + zstd_time_compr_dk + ", " + mfc_time_compr_dk + ", " + spr_time_compr_dk + " - " + FM_time_compr_dk);

            //decompression time
            DPfile_time.println(bzip2_time_decompr_dk + ", " + lz4_time_decompr_dk + ", " + zstd_time_decompr_dk + ", " + mfc_time_decompr_dk + ", " + spr_time_decompr_dk + " - " + FM_time_decompr_dk);
        } else {
            //size
            DPfile_size.println(uncompr_size_dk + ", " + bzip2_size_dk + ", " + lz4_size_dk + ", " + zstd_size_dk + ", " + mfc_size_dk + ", " + spr_size_dk);

            //compression time
            DPfile_time.println(bzip2_time_compr_dk + ", " + lz4_time_compr_dk + ", " + zstd_time_compr_dk + ", " + mfc_time_compr_dk + ", " + spr_time_compr_dk);

            //decompression time
            DPfile_time.println(bzip2_time_decompr_dk + ", " + lz4_time_decompr_dk + ", " + zstd_time_decompr_dk + ", " + mfc_time_decompr_dk + ", " + spr_time_decompr_dk);
        }
    }


    private static void Manage_Sk_fragments(String path_in, String exp, int DP_scenario, int num_blocks, PrintWriter DPfile_size, PrintWriter DPfile_time, int debug) throws IOException, InterruptedException {

        //FIXME: read .nfo file
        String l;
        int n_blocks = num_blocks;  // in input we have passed 1, set to minimum value

        if (n_blocks > 1) {     // read the number of ESS blocks to be read, in case (num_blocks > 1)
            // Reading .nfo information to extract information: total_rows, k, num_rows_file, num_fragments.
            BufferedReader Header_NFO = null;
            try {
                Header_NFO = new BufferedReader(new FileReader(path_in + FS_SEPARATOR + exp + nfo_EXT));
                if (debug==1) System.out.println("Read .nfo info from: " + path_in + FS_SEPARATOR + exp + nfo_EXT);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
            // Extract info from .nfo file
            l = Header_NFO.readLine();
            n_blocks = Integer.parseInt(l);
            if (debug==1) System.out.println("ESS num_blocks = " + n_blocks + "\n");
            // Chiusura File
            Header_NFO.close();
        }

        String ess_file = null;

        //FIXME: opening .fasta fragments for reading dimensions and cumulative sum of them.
        long uncompr_size_Sk = 0;
        if (n_blocks == 1) {

            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;

            // compute size of fasta block ESS
            File file = new File(ess_file);
            if (file.exists() && file.isFile()) {
                uncompr_size_Sk += file.length();
            } else {
                System.out.println(ess_file + " not a file or non-existent file !!!.");
            }

        } else {
            for (int i = 0; i < n_blocks; i++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + i + FASTA_EXT;

                // compute size of fasta block ESS
                File file = new File(ess_file);
                if (file.exists() && file.isFile()) {
                    uncompr_size_Sk += file.length();
                } else {
                    System.out.println(ess_file + " not a file or non-existent file !!!.");
                }
            }
        }

        //FIXME: compress/uncompress phase with generic and specialized Dk compression methods
        String cmd;

        //FIXME: bzip2 compression
        long bzip2_start_compr = System.currentTimeMillis();
        long bzip2_size_Sk = 0;

        if (n_blocks == 1) {    // for S
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;

            if (debug==1) System.out.println("bzip2 -z -k " + ess_file + " --best");

            ProcessBuilder builderC = new ProcessBuilder("bzip2", "-v", "-z", "-k", ess_file, "--best");
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update bzip2_size [bytes] with size of current fasta fragment
            File file = new File(ess_file + bzip2_c_EXT);
            if (file.exists() && file.isFile()) {
                bzip2_size_Sk += file.length();
            } else {
                System.out.println(ess_file + bzip2_c_EXT + " not a file or non-existent file !!!.");
            }
        } else {  // for every S_i

            for (int j = 0; j < n_blocks; j++) {

                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;

                if (debug==1) System.out.println("bzip2 -z -k " + ess_file + " --best");

                ProcessBuilder builderC = new ProcessBuilder("bzip2", "-v", "-z", "-k", ess_file, "--best");
                builderC.redirectErrorStream(true);
                final Process processC = builderC.start();
                // Watch the process
                watch(processC);
                processC.waitFor();

                //update bzip2_size [bytes] with size of current fasta fragment
                File file = new File(ess_file + bzip2_c_EXT);
                if (file.exists() && file.isFile()) {
                    bzip2_size_Sk += file.length();
                } else {
                    System.out.println(ess_file + bzip2_c_EXT + " not a file or non-existent file !!!.");
                }
            }
        }
        long bzip2_time_compr_Sk = System.currentTimeMillis() - bzip2_start_compr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: bzip2 decompression
        long bzip2_start_decompr = System.currentTimeMillis();

        if (n_blocks == 1) {    // for S
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;

            if (debug==1) System.out.println("bzip2 -d -k -f " + ess_file + bzip2_c_EXT);

            ProcessBuilder builderD = new ProcessBuilder("bzip2", "-v", "-d", "-k", "-f", ess_file + bzip2_c_EXT);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
        } else {
            for (int j = 0; j < num_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;

                if (debug==1) System.out.println("bzip2 -d -k -f " + ess_file + bzip2_c_EXT);

                ProcessBuilder builderD = new ProcessBuilder("bzip2", "-v", "-d", "-k", "-f", ess_file + bzip2_c_EXT);
                builderD.redirectErrorStream(true);
                final Process processD = builderD.start();
                // Watch the process
                watch(processD);
                processD.waitFor();
            }
        }
        long bzip2_time_decompr_Sk = System.currentTimeMillis() - bzip2_start_decompr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: remove bzip2 compression files
        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;
            remove_file(ess_file + bzip2_c_EXT, debug);
        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;
                remove_file(ess_file + bzip2_c_EXT, debug);
            }
        }


        //FIXME: lz4 compression
        long lz4_start_compr = System.currentTimeMillis();
        long lz4_size_Sk = 0;

        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;
            ProcessBuilder builderC = new ProcessBuilder("lz4", "-z", "-f", "-k", "-9", ess_file, ess_file + lz4_c_EXT);
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update lz4_size [bytes] with size of current fasta fragment
            File file = new File(ess_file + lz4_c_EXT);
            if (file.exists() && file.isFile()) {
                lz4_size_Sk += file.length();
            } else {
                System.out.println(ess_file + lz4_c_EXT + " not a file or non-existent file !!!.");
            }
        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;

                if (debug==1) System.out.println("lz4 -z -f -k -9 " + ess_file + " " + ess_file + lz4_c_EXT);

                ProcessBuilder builderC = new ProcessBuilder("lz4", "-z", "-f", "-k", "-9", ess_file, ess_file + lz4_c_EXT);
                builderC.redirectErrorStream(true);
                final Process processC = builderC.start();
                // Watch the process
                watch(processC);
                processC.waitFor();

                //update lz4_size [bytes] with size of current fasta fragment
                File file = new File(ess_file + lz4_c_EXT);
                if (file.exists() && file.isFile()) {
                    lz4_size_Sk += file.length();
                } else {
                    System.out.println(ess_file + lz4_c_EXT + " not a file or non-existent file !!!.");
                }
            }
        }
        long lz4_time_compr_Sk = System.currentTimeMillis() - lz4_start_compr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: lz4 decompression
        long lz4_start_decompr = System.currentTimeMillis();

        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;

            if (debug==1) System.out.println("lz4 -d -k " + ess_file + lz4_c_EXT + " " + ess_file);

            ProcessBuilder builderD = new ProcessBuilder("lz4", "-d", "-k", "-f", ess_file + lz4_c_EXT, ess_file);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;

                if (debug==1) System.out.println("lz4 -d -k " + ess_file + lz4_c_EXT + " " + ess_file);

                ProcessBuilder builderD = new ProcessBuilder("lz4", "-d", "-k", "-f", ess_file + lz4_c_EXT, ess_file);
                builderD.redirectErrorStream(true);
                final Process processD = builderD.start();
                // Watch the process
                watch(processD);
                processD.waitFor();
            }
        }
        long lz4_time_decompr_Sk = System.currentTimeMillis() - lz4_start_decompr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: remove lz4 compression files
        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;
            remove_file(ess_file + lz4_c_EXT, debug);
        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;
                remove_file(ess_file + lz4_c_EXT, debug);
            }
        }


        //FIXME: zstd compression
        long zstd_start_compr = System.currentTimeMillis();
        long zstd_size_Sk = 0;

        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;

            if (debug==1) System.out.println("zstd -19 " + ess_file + " --format=zstd");

            ProcessBuilder builderC = new ProcessBuilder("zstd", "-19", ess_file, "--format=zstd");
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update zstd_size [bytes] with size of current fasta fragment
            File file = new File(ess_file + zstd_c_EXT);
            if (file.exists() && file.isFile()) {
                zstd_size_Sk += file.length();
            } else {
                System.out.println(ess_file + zstd_c_EXT + " not a file or non-existent file !!!.");
            }
        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;

                if (debug==1) System.out.println("zstd -19 " + ess_file + " --format=zstd");

                ProcessBuilder builderC = new ProcessBuilder("zstd", "-19", ess_file, "--format=zstd");
                builderC.redirectErrorStream(true);
                final Process processC = builderC.start();
                // Watch the process
                watch(processC);
                processC.waitFor();

                //update zstd_size [bytes] with size of current fasta fragment
                File file = new File(ess_file + zstd_c_EXT);
                if (file.exists() && file.isFile()) {
                    zstd_size_Sk += file.length();
                } else {
                    System.out.println(ess_file + zstd_c_EXT + " not a file or non-existent file !!!.");
                }
            }
        }
        long zstd_time_compr_Sk = System.currentTimeMillis() - zstd_start_compr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: zstd decompression
        long zstd_start_decompr = System.currentTimeMillis();

        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;

            if (debug==1) System.out.println("zstd -d -f " + ess_file + zstd_c_EXT);

            ProcessBuilder builderD = new ProcessBuilder("zstd", "-d", "-f", ess_file + zstd_c_EXT);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;

                if (debug==1) System.out.println("zstd -d -f " + ess_file + zstd_c_EXT);

                ProcessBuilder builderD = new ProcessBuilder("zstd", "-d", "-f", ess_file + zstd_c_EXT);
                builderD.redirectErrorStream(true);
                final Process processD = builderD.start();
                // Watch the process
                watch(processD);
                processD.waitFor();
            }
        }
        long zstd_time_decompr_Sk = System.currentTimeMillis() - zstd_start_decompr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: remove zstd compression files
        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;
            remove_file(ess_file + zstd_c_EXT, debug);
        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;
                remove_file(ess_file + zstd_c_EXT, debug);
            }
        }


        //FIXME: MFC compression
        long mfc_start_compr = System.currentTimeMillis();
        long mfc_size_Sk = 0;

        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;

            if (debug==1) System.out.println("./essAuxMFCompressC -3 -p 2 " + ess_file);

            ProcessBuilder builderC = new ProcessBuilder("./essAuxMFCompressC", "-3", "p 2", "-v", ess_file);
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update mfc_size [bytes] with size of current fasta fragment
            File file = new File(ess_file + mfc_c_EXT);
            if (file.exists() && file.isFile()) {
                mfc_size_Sk += file.length();
            } else {
                System.out.println(ess_file + mfc_c_EXT + " not a file or non-existent file !!!.");
            }
        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;

                if (debug==1) System.out.println("./essAuxMFCompressC -3 -p 2 " + ess_file);

                ProcessBuilder builderC = new ProcessBuilder("./essAuxMFCompressC", "-3", "p 2", "-v", ess_file);
                builderC.redirectErrorStream(true);
                final Process processC = builderC.start();
                // Watch the process
                watch(processC);
                processC.waitFor();

                //update mfc_size [bytes] with size of current fasta fragment
                File file = new File(ess_file + mfc_c_EXT);
                if (file.exists() && file.isFile()) {
                    mfc_size_Sk += file.length();
                } else {
                    System.out.println(ess_file + mfc_c_EXT + " not a file or non-existent file !!!.");
                }
            }
        }
        long mfc_time_compr_Sk = System.currentTimeMillis() - mfc_start_compr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: mfc decompression
        long mfc_start_decompr = System.currentTimeMillis();

        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;

            if (debug==1) System.out.println("./essAuxMFCompressD " + ess_file + mfc_c_EXT);

            ProcessBuilder builderD = new ProcessBuilder("./essAuxMFCompressD", "-v", ess_file + mfc_c_EXT);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();

        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;

                if (debug==1) System.out.println("./essAuxMFCompressD " + ess_file + mfc_c_EXT);

                ProcessBuilder builderD = new ProcessBuilder("./essAuxMFCompressD", "-v", ess_file + mfc_c_EXT);
                builderD.redirectErrorStream(true);
                final Process processD = builderD.start();
                // Watch the process
                watch(processD);
                processD.waitFor();
            }
        }
        long mfc_time_decompr_Sk = System.currentTimeMillis() - mfc_start_decompr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: remove mfc compression files
        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;
            remove_file(ess_file + mfc_c_EXT, debug);
            remove_file(ess_file + mfc_c_EXT + mfc_d_EXT, debug);
        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;
                remove_file(ess_file + mfc_c_EXT, debug);
                remove_file(ess_file + mfc_c_EXT + mfc_d_EXT, debug);
            }
        }


        //FIXME: Spring compression
        long spr_start_compr = System.currentTimeMillis();
        long spr_size_Sk = 0;

        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;

            if (debug==1) System.out.println("./spring -c (--fasta-input) (-l) -i " + ess_file + " -o " + ess_file + spring_c_EXT);

            ProcessBuilder builderC = new ProcessBuilder("./spring", "-c", "--fasta-input", "-l", "-i", ess_file, "-o", ess_file + spring_c_EXT);
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update spr_size [bytes] with size of current fasta fragment
            File file = new File(ess_file + spring_c_EXT);
            if (file.exists() && file.isFile()) {
                spr_size_Sk += file.length();
            } else {
                System.out.println(ess_file + spring_c_EXT + " not a file or non-existent file !!!.");
            }
        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;

                if (debug==1) System.out.println("./spring -c (--fasta-input) (-l) -i " + ess_file + " -o " + ess_file + spring_c_EXT);

                ProcessBuilder builderC = new ProcessBuilder("./spring", "-c", "--fasta-input", "-l", "-i", ess_file, "-o", ess_file + spring_c_EXT);
                builderC.redirectErrorStream(true);
                final Process processC = builderC.start();
                // Watch the process
                watch(processC);
                processC.waitFor();

                //update spr_size [bytes] with size of current fasta fragment
                File file = new File(ess_file + spring_c_EXT);
                if (file.exists() && file.isFile()) {
                    spr_size_Sk += file.length();
                } else {
                    System.out.println(ess_file + spring_c_EXT + " not a file or non-existent file !!!.");
                }
            }
        }
        long spr_time_compr_Sk = System.currentTimeMillis() - spr_start_compr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: spring decompression
        long spr_start_decompr = System.currentTimeMillis();

        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;

            if (debug==1) System.out.println("./spring -d -i " + ess_file + spring_c_EXT + " -o " + ess_file + spring_d_EXT);

            ProcessBuilder builderD = new ProcessBuilder("./spring", "-d", "-i", ess_file + spring_c_EXT, "-o", ess_file + spring_d_EXT);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;

                if (debug==1) System.out.println("./spring -d -i " + ess_file + spring_c_EXT + " -o " + ess_file + spring_d_EXT);

                ProcessBuilder builderD = new ProcessBuilder("./spring", "-d", "-i", ess_file + spring_c_EXT, "-o", ess_file + spring_d_EXT);
                builderD.redirectErrorStream(true);
                final Process processD = builderD.start();
                // Watch the process
                watch(processD);
                processD.waitFor();
            }
        }
        long spr_time_decompr_Sk = System.currentTimeMillis() - spr_start_decompr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: remove spring compression files
        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;
            remove_file(ess_file + spring_c_EXT, debug);
            remove_file(ess_file + spring_d_EXT, debug);
        } else {
            for (int j = 0; j < n_blocks; j++) {
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;
                remove_file(ess_file + spring_c_EXT, debug);
                remove_file(ess_file + spring_d_EXT, debug);
            }
        }


        // Cleanup Garbage Collector
        System.gc();

        /*

        //FIXME: remove FM-index compression/decompression  files .fasta.dec
        if (n_blocks == 1) {
            ess_file = path_in + FS_SEPARATOR + exp + FASTA_EXT;
            remove_file(ess_file + FMind_d_EXT, debug);
        } else {
            for (int j = 0; j < n_blocks; j++) {
                //FIXME: remove file .fm9
                ess_file = path_in + FS_SEPARATOR + exp + "_Block_" + j + FASTA_EXT;
                remove_file(ess_file + FMind_d_EXT, debug);
            }
        }
        */


        //FIXME: write to DP3 file size and time information !!!.
        //size
        DPfile_size.println(uncompr_size_Sk + ", " + bzip2_size_Sk + ", " + lz4_size_Sk + ", " + zstd_size_Sk + ", " + mfc_size_Sk + ", " + spr_size_Sk /* + " - " + uncompr_size_Sk + ", " + FM_size_Sk */);

        //compression time
        DPfile_time.println(bzip2_time_compr_Sk + ", " + lz4_time_compr_Sk + ", " + zstd_time_compr_Sk + ", " + mfc_time_compr_Sk + ", " + spr_time_compr_Sk /* + " - " + FM_time_compr_Sk */);

        //decompression time
        DPfile_time.println(bzip2_time_decompr_Sk + ", " + lz4_time_decompr_Sk + ", " + zstd_time_decompr_Sk + ", " + mfc_time_decompr_Sk + ", " + spr_time_decompr_Sk /* + " - " + FM_time_decompr_Sk */);

    }


    private static void Manage_Fa_fragments(String path_in, String exp, int DP_scenario, int approx, PrintWriter DPfile_size, PrintWriter DPfile_time, int gap, int debug) throws IOException, InterruptedException {

        //FIXME: read .nfo file
        String l;
        int num_rows_file = 0;
        long total_rows = 0;
        int num_fragments = 0;
        int k = 0;

        // Reading .nfo information to extract information: total_rows, k, num_rows_file, num_fragments
        String header_file = null;
        String approx_Suffix = "";

        if (approx == 0) {
            header_file = path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + nfo_EXT;
        } else {
            header_file = path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + "_approx" + nfo_EXT;
            approx_Suffix = "_approx";
        }

        BufferedReader Header_NFO = null;
        try {
            Header_NFO = new BufferedReader(new FileReader(header_file));
            System.out.println("MAIN" + "\t" + "Read .nfo info from: " + header_file);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Extract info from .nfo file
        l = Header_NFO.readLine();
        total_rows = Long.parseLong(l);

        l = Header_NFO.readLine();
        k = Integer.parseInt(l);

        l = Header_NFO.readLine();
        num_rows_file = Integer.parseInt(l);

        l = Header_NFO.readLine();
        num_fragments = Integer.parseInt(l);

        // Header File Closure
        Header_NFO.close();

        // Auxiliary Variables for distinguish case Fk and case Gk
        String gap_Suffix = "";

        String scenario_Lemire = "a";

        if (gap == 1) {
            gap_Suffix = "-Gk";
            scenario_Lemire = "b";
        }

        //FIXME: opening .txt fragments for size reading and cumulative sum of dimensions.
        //compute Fk fragments cumulative size
        long uncompr_size_fk = 0;
        for (int i = 0; i < num_fragments; i++) {

            // compute size of txt fragment
            File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + i + DSK_EXT);
            if (file.exists() && file.isFile()) {
                uncompr_size_fk += file.length();
            } else {
                System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + i + DSK_EXT + " not a file or non-existent file !!!.");
            }
        }

        //FIXME: compression/decompression phase with generic and specialized Fk compression methods
        String cmd;

        //FIXME: bzip2 compression
        long bzip2_start_compr = System.currentTimeMillis();
        long bzip2_size_fk = 0;

        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("bzip2 -z -k " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + " --best");

            ProcessBuilder builderC = new ProcessBuilder("bzip2", "-v", "-z", "-k", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT, "--best");
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update bzip2_size [bytes] with size of current txt fragment
            File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + bzip2_c_EXT);
            if (file.exists() && file.isFile()) {
                bzip2_size_fk += file.length();
            } else {
                System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + bzip2_c_EXT + " not a file or non-existent file !!!.");
            }
        }
        long bzip2_time_compr_fk = System.currentTimeMillis() - bzip2_start_compr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: bzip2 decompression
        long bzip2_start_decompr = System.currentTimeMillis();

        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("bzip2 -d -k -f " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + bzip2_c_EXT);

            ProcessBuilder builderD = new ProcessBuilder("bzip2", "-v", "-d", "-k", "-f", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + bzip2_c_EXT);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();
        }
        long bzip2_time_decompr_fk = System.currentTimeMillis() - bzip2_start_decompr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: remove bzip2 compression files
        for (int j = 0; j < num_fragments; j++) {
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + bzip2_c_EXT, debug);
        }


        //TODO: lz4 compression
        long lz4_start_compr = System.currentTimeMillis();
        long lz4_size_fk = 0;

        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("lz4 -z -f -k -9 --quiet " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + " " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + lz4_c_EXT);

            ProcessBuilder builderC = new ProcessBuilder("lz4", "-v", "-z", "-f", "-k", "-9", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT, path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + lz4_c_EXT);
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update lz4_size [bytes] with size of current txt fragment
            File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + lz4_c_EXT);
            if (file.exists() && file.isFile()) {
                lz4_size_fk += file.length();
            } else {
                System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + lz4_c_EXT + " not a file or non-existent file !!!.");
            }

        }
        long lz4_time_compr_fk = System.currentTimeMillis() - lz4_start_compr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: lz4 decompression
        long lz4_start_decompr = System.currentTimeMillis();

        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("lz4 -d -k -f --quiet " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + lz4_c_EXT + " " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT);

            ProcessBuilder builderD = new ProcessBuilder("lz4", "-v", "-d", "-k", "-f", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + lz4_c_EXT, path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();

        }

        long lz4_time_decompr_fk = System.currentTimeMillis() - lz4_start_decompr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: remove bzip2 compression files
        for (int j = 0; j < num_fragments; j++) {
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + lz4_c_EXT, debug);
        }


        //FIXME: zstd compression
        long zstd_start_compr = System.currentTimeMillis();
        long zstd_size_fk = 0;

        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("zstd -19 " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + " --format=zstd");

            ProcessBuilder builderC = new ProcessBuilder("zstd", "-19", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT, "--format=zstd");
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update zstd_size [bytes] with size of current fasta fragment
            File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + zstd_c_EXT);
            if (file.exists() && file.isFile()) {
                zstd_size_fk += file.length();
            } else {
                System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + zstd_c_EXT + " not a file or non-existent file !!!.");
            }

        }
        long zstd_time_compr_fk = System.currentTimeMillis() - zstd_start_compr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: zstd decompression
        long zstd_start_decompr = System.currentTimeMillis();

        for (int j = 0; j < num_fragments; j++) {

            if (debug==1) System.out.println("zstd -d -f " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + zstd_c_EXT);

            ProcessBuilder builderD = new ProcessBuilder("zstd", "-d", "-f", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + zstd_c_EXT);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();

        }
        long zstd_time_decompr_fk = System.currentTimeMillis() - zstd_start_decompr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: remove zstd compression files
        for (int j = 0; j < num_fragments; j++) {
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + zstd_c_EXT, debug);
        }


        //FIXME: Lemire BIC compression
        long BIC_start_compr = System.currentTimeMillis();
        long BIC_size_fk = 0;
        Lemire_Fk_BIC_Compr_exp.execute(path_in + FS_SEPARATOR, exp + gap_Suffix, scenario_Lemire, approx_Suffix);
        if (debug==1) System.out.println("Apply Lemire BIC compression to: " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + " fragments.");

        //update BIC_size [bytes] with size of current txt fragment
        for (int j = 0; j < num_fragments; j++) {
            File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + BIC_c_EXT);
            if (file.exists() && file.isFile()) {
                BIC_size_fk += file.length();
            } else {
                System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + BIC_c_EXT + " not a file or non-existent file !!!.");
            }
        }
        long BIC_time_compr_fk = System.currentTimeMillis() - BIC_start_compr;

        // Cleanup Garbage Collector
        System.gc();


        //FIXME: Lemire BIC Decompression
        long BIC_start_decompr = System.currentTimeMillis();
        Lemire_Fk_BIC_Decompr_exp.execute(path_in + FS_SEPARATOR, exp + gap_Suffix, scenario_Lemire, approx_Suffix);
        if (debug==1) System.out.println("Apply Lemire BIC decompression to: " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + " fragments.");
        long BIC_time_decompr_fk = System.currentTimeMillis() - BIC_start_decompr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: remove BIC compression .bic and decompression .d files
        for (int j = 0; j < num_fragments; j++) {
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + BIC_c_EXT, debug);
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + BIC_c_EXT + BIC_d_EXT, debug);
        }


        //FIXME: Lemire OptPFOR compression
        long OptPFOR_start_compr = System.currentTimeMillis();
        long OptPFOR_size_fk = 0;

        Lemire_Fk_OptPFOR_Compr_exp.execute(path_in + FS_SEPARATOR, exp + gap_Suffix, scenario_Lemire, approx_Suffix);
        if (debug==1) System.out.println("Apply Lemire OptPFOR compression to: " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + " fragments.");

        //update OptPFOR_size [bytes] with size of current txt fragment
        for (int j = 0; j < num_fragments; j++) {
            File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + OptPFOR_c_EXT);
            if (file.exists() && file.isFile()) {
                OptPFOR_size_fk += file.length();
            } else {
                System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + OptPFOR_c_EXT + " not a file or non-existent file !!!.");
            }
        }
        long OptPFOR_time_compr_fk = System.currentTimeMillis() - OptPFOR_start_compr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: Lemire OptPFOR Decompression
        long OptPFOR_start_decompr = System.currentTimeMillis();
        Lemire_Fk_OptPFOR_Decompr_exp.execute(path_in + FS_SEPARATOR, exp + gap_Suffix, scenario_Lemire, approx_Suffix);
        if (debug==1) System.out.println("Apply Lemire OptPFOR decompression to: " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + " fragments.");
        long OptPFOR_time_decompr_fk = System.currentTimeMillis() - OptPFOR_start_decompr;

        // Cleanup Garbage Collector
        System.gc();

        //FIXME: remove OptPFOR compression .bic and decompression .d files
        for (int j = 0; j < num_fragments; j++) {
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + OptPFOR_c_EXT, debug);
            remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + DSK_EXT + OptPFOR_c_EXT + OptPFOR_d_EXT, debug);
        }


        if (DP_scenario == 0) { // Build BCSF only for the DP0 Sccenario, and without gaps

            if (gap == 0) {
                //FIXME: BCSF - input dictionary to locom: without gaps.
                long Dict_size_DSK = 0;

                // BCSF size and  construction time of the succinct data structure.
                long BCSF_size_build_CSF_fk = 0, BCSF_size_build_BLOOM_fk = 0, BCSF_size_build_fk = 0, BCSF_time_build_fk = 0;

                // Size of BCSF compressions using standard BZIP2, LZ4, and ZSTD methods, along with compression and decompression times.
                long BZIP2_BCSF_size_fk = 0, LZ4_BCSF_size_fk = 0, ZSTD_BCSF_size_fk = 0;
                long BZIP2_BCSF_time_compr_fk = 0, BZIP2_BCSF_time_startC_fk = 0, LZ4_BCSF_time_compr_fk = 0, LZ4_BCSF_time_startC_fk = 0, ZSTD_BCSF_time_compr_fk = 0, ZSTD_BCSF_time_startC_fk = 0;
                long BZIP2_BCSF_time_decompr_fk = 0, BZIP2_BCSF_time_startD_fk = 0, LZ4_BCSF_time_decompr_fk = 0, LZ4_BCSF_time_startD_fk = 0, ZSTD_BCSF_time_decompr_fk = 0, ZSTD_BCSF_time_startD_fk = 0;


                //FIXME: build BCSF data structure.
                long BCSF_time_startbuild_fk = System.currentTimeMillis();

                for (int j = 0; j < num_fragments; j++) {

                    //FIXME: create DSK_i with Dk_i and Fk_i fragments
                    Create_DSK_Fragment(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + "_" + j);

                    if (debug==1) System.out.println("python3 locom.py build -i " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + "_DSK" + DSK_EXT + " -o " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j);

                    ProcessBuilder builderC = new ProcessBuilder("python3", "locom.py", "build", "-i", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + "_DSK" + DSK_EXT, "-o", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j);
                    builderC.redirectErrorStream(true);
                    final Process processC = builderC.start();
                    // Watch the process
                    watch(processC);
                    processC.waitFor();

                    //Get size of input file to locom
                    File file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + "_DSK" + DSK_EXT);
                    if (file.exists() && file.isFile()) {
                        Dict_size_DSK += file.length();
                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + "_DSK" + DSK_EXT + " not a file or non-existent file !!!.");
                    }

                    //update BCSF_size_CSF_fk [bytes] with size of current txt fragment
                    file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT);
                    if (file.exists() && file.isFile()) {
                        BCSF_size_build_CSF_fk += file.length();
                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + " not a file or non-existent file !!!.");
                    }

                    //update BCSF_size_BLOOM_fk [bytes] with size of current txt fragment
                    file = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT);
                    if (file.exists() && file.isFile()) {
                        BCSF_size_build_BLOOM_fk += file.length();
                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + " not a file or non-existent file !!!.");
                    }
                }
                BCSF_size_build_fk += BCSF_size_build_CSF_fk + BCSF_size_build_BLOOM_fk;
                BCSF_time_build_fk = System.currentTimeMillis() - BCSF_time_startbuild_fk;

                // Cleanup Garbage Collector
                System.gc();

                //FIXME: BCSF(Fk) compression with BZIP2
                BZIP2_BCSF_time_startC_fk = System.currentTimeMillis();
                for (int j = 0; j < num_fragments; j++) {

                    File fileCSF = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT);
                    if (fileCSF.exists() && fileCSF.isFile()) {

                        //compressione CSF part
                        if (debug==1) System.out.println("bzip2 -z -k " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + " --best");

                        ProcessBuilder builderC = new ProcessBuilder("bzip2", "-v", "-z", "-k", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT, "--best");
                        builderC.redirectErrorStream(true);
                        final Process processC = builderC.start();
                        // Watch the process
                        watch(processC);
                        processC.waitFor();

                        File fileBZ2_CSF = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + bzip2_c_EXT);
                        if (fileBZ2_CSF.exists() && fileBZ2_CSF.isFile()) {
                            //update FM_size [bytes] with size of current fragment
                            BZIP2_BCSF_size_fk += fileBZ2_CSF.length();
                        }
                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + " not a file or non-existent file !!!.");
                    }

                    File fileBLOOM = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT);
                    if (fileBLOOM.exists() && fileBLOOM.isFile()) {

                        //compressione BLOOM part
                        if (debug==1) System.out.println("bzip2 -z -k " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + " --best");

                        ProcessBuilder builderC = new ProcessBuilder("bzip2", "-v", "-z", "-k", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT, "--best");
                        builderC.redirectErrorStream(true);
                        final Process processC = builderC.start();
                        // Watch the process
                        watch(processC);
                        processC.waitFor();

                        File fileBZ2_BLOOM = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + bzip2_c_EXT);
                        if (fileBZ2_BLOOM.exists() && fileBZ2_BLOOM.isFile()) {
                            //update FM_size [bytes] with size of current fragment
                            BZIP2_BCSF_size_fk += fileBZ2_BLOOM.length();
                        }
                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + " not a file or non-existent file !!!.");
                    }
                }
                BZIP2_BCSF_time_compr_fk = System.currentTimeMillis() - BZIP2_BCSF_time_startC_fk;

                // Cleanup Garbage Collector
                System.gc();

                //FIXME: BCSF(Fk) decompression with BZIP2
                BZIP2_BCSF_time_startD_fk = System.currentTimeMillis();

                for (int j = 0; j < num_fragments; j++) {

                    File fileCSF = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + bzip2_c_EXT);
                    if (fileCSF.exists() && fileCSF.isFile()) {

                        //decompressione CSF part
                        if (debug==1) System.out.println("bzip2 -d -k -f " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + bzip2_c_EXT);
                        ProcessBuilder builderD = new ProcessBuilder("bzip2", "-v", "-d", "-k", "-f", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + bzip2_c_EXT);

                        builderD.redirectErrorStream(true);
                        final Process processD = builderD.start();
                        // Watch the process
                        watch(processD);
                        processD.waitFor();

                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + bzip2_c_EXT + " not a file or non-existent file !!!.");
                    }


                    File fileBLOOM = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + bzip2_c_EXT);
                    if (fileBLOOM.exists() && fileBLOOM.isFile()) {

                        //decompressione BLOOM part
                        if (debug==1) System.out.println("bzip2 -d -k -f " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + bzip2_c_EXT);

                        ProcessBuilder builderD = new ProcessBuilder("bzip2", "-v", "-d", "-k", "-f", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + bzip2_c_EXT);
                        builderD.redirectErrorStream(true);
                        final Process processD = builderD.start();
                        // Watch the process
                        watch(processD);
                        processD.waitFor();
                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + bzip2_c_EXT + " not a file or non-existent file !!!.");
                    }
                }
                BZIP2_BCSF_time_decompr_fk = System.currentTimeMillis() - BZIP2_BCSF_time_startD_fk;

                // Cleanup Garbage Collector
                System.gc();

                //FIXME: remove BZIP2 BCSF compression/decompression  temporary files
                for (int j = 0; j < num_fragments; j++) {
                    remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + bzip2_c_EXT, debug);
                    remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + bzip2_c_EXT, debug);
                }


                //FIXME: BCSF(Fk) compression with LZ4
                LZ4_BCSF_time_startC_fk = System.currentTimeMillis();
                for (int j = 0; j < num_fragments; j++) {

                    File fileCSF = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT);
                    if (fileCSF.exists() && fileCSF.isFile()) {

                        //compressione CSF part
                        if (debug==1) System.out.println("lz4 -z -f -k -9 --quiet " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + " " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + lz4_c_EXT);

                        ProcessBuilder builderC = new ProcessBuilder("lz4", "-v", "-z", "-f", "-k", "-9", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT, path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + lz4_c_EXT);
                        builderC.redirectErrorStream(true);
                        final Process processC = builderC.start();
                        // Watch the process
                        watch(processC);
                        processC.waitFor();

                        File fileLZ4_CSF = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + lz4_c_EXT);
                        if (fileLZ4_CSF.exists() && fileLZ4_CSF.isFile()) {
                            //update FM_size [bytes] with size of current fragment
                            LZ4_BCSF_size_fk += fileLZ4_CSF.length();
                        }
                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + " not a file or non-existent file !!!.");
                    }

                    File fileBLOOM = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT);
                    if (fileBLOOM.exists() && fileBLOOM.isFile()) {

                        //compressione BLOOM part
                        if (debug==1) System.out.println("lz4 -z -f -k -9 --quiet " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + " " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + lz4_c_EXT);

                        ProcessBuilder builderC = new ProcessBuilder("lz4", "-v", "-z", "-f", "-k", "-9", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT, path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + lz4_c_EXT);
                        builderC.redirectErrorStream(true);
                        final Process processC = builderC.start();
                        // Watch the process
                        watch(processC);
                        processC.waitFor();

                        File fileLZ4_BLOOM = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + lz4_c_EXT);
                        if (fileLZ4_BLOOM.exists() && fileLZ4_BLOOM.isFile()) {
                            //update FM_size [bytes] with size of current fragment
                            LZ4_BCSF_size_fk += fileLZ4_BLOOM.length();
                        }
                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + " not a file or non-existent file !!!.");
                    }
                }
                LZ4_BCSF_time_compr_fk = System.currentTimeMillis() - LZ4_BCSF_time_startC_fk;

                // Cleanup Garbage Collector
                System.gc();

                //FIXME: BCSF(Fk) decompression with LZ4
                LZ4_BCSF_time_startD_fk = System.currentTimeMillis();

                for (int j = 0; j < num_fragments; j++) {

                    File fileCSF = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + lz4_c_EXT);
                    if (fileCSF.exists() && fileCSF.isFile()) {

                        //decompressione CSF part
                        if (debug==1) System.out.println("lz4 -d -k -f --quiet " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + lz4_c_EXT + " " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT);

                        ProcessBuilder builderD = new ProcessBuilder("lz4", "-v", "-d", "-k", "-f", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + lz4_c_EXT, path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT);
                        builderD.redirectErrorStream(true);
                        final Process processD = builderD.start();
                        // Watch the process
                        watch(processD);
                        processD.waitFor();

                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + lz4_c_EXT + " not a file or non-existent file !!!.");
                    }


                    File fileBLOOM = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + lz4_c_EXT);
                    if (fileBLOOM.exists() && fileBLOOM.isFile()) {

                        //decompressione BLOOM part
                        if (debug==1) System.out.println("lz4 -d -k -f --quiet " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + lz4_c_EXT + " " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT);

                        ProcessBuilder builderD = new ProcessBuilder("lz4", "-v", "-d", "-k", "-f", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + lz4_c_EXT, path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT);
                        builderD.redirectErrorStream(true);
                        final Process processD = builderD.start();
                        // Watch the process
                        watch(processD);
                        processD.waitFor();
                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + lz4_c_EXT + " not a file or non-existent file !!!.");
                    }
                }
                LZ4_BCSF_time_decompr_fk = System.currentTimeMillis() - LZ4_BCSF_time_startD_fk;

                // Cleanup Garbage Collector
                System.gc();

                //FIXME: remove LZ4 BCSF compression/decompression  temporary files
                for (int j = 0; j < num_fragments; j++) {
                    remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + lz4_c_EXT, debug);
                    remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + lz4_c_EXT, debug);
                }


                //FIXME: BCSF(Fk) compression with ZSTD
                ZSTD_BCSF_time_startC_fk = System.currentTimeMillis();
                for (int j = 0; j < num_fragments; j++) {

                    File fileCSF = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT);
                    if (fileCSF.exists() && fileCSF.isFile()) {

                        //compressione CSF part
                        if (debug==1) System.out.println("zstd -19 " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + " --format=zstd");

                        ProcessBuilder builderC = new ProcessBuilder("zstd", "-19", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT, "--format=zstd");
                        builderC.redirectErrorStream(true);
                        final Process processC = builderC.start();
                        // Watch the process
                        watch(processC);
                        processC.waitFor();

                        File fileZSTD_CSF = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + zstd_c_EXT);
                        if (fileZSTD_CSF.exists() && fileZSTD_CSF.isFile()) {
                            //update FM_size [bytes] with size of current fragment
                            ZSTD_BCSF_size_fk += fileZSTD_CSF.length();
                        }
                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + " not a file or non-existent file !!!.");
                    }

                    File fileBLOOM = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT);
                    if (fileBLOOM.exists() && fileBLOOM.isFile()) {

                        //compressione BLOOM part
                        if (debug==1) System.out.println("zstd -19 " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + " --format=zstd");

                        ProcessBuilder builderC = new ProcessBuilder("zstd", "-19", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT, "--format=zstd");
                        builderC.redirectErrorStream(true);
                        final Process processC = builderC.start();
                        // Watch the process
                        watch(processC);
                        processC.waitFor();

                        File fileZSTD_BLOOM = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + zstd_c_EXT);
                        if (fileZSTD_BLOOM.exists() && fileZSTD_BLOOM.isFile()) {
                            //update FM_size [bytes] with size of current fragment
                            ZSTD_BCSF_size_fk += fileZSTD_BLOOM.length();
                        }
                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + " not a file or non-existent file !!!.");
                    }
                }
                ZSTD_BCSF_time_compr_fk = System.currentTimeMillis() - ZSTD_BCSF_time_startC_fk;

                // Cleanup Garbage Collector
                System.gc();

                //FIXME: BCSF(Fk) decompression with ZSTD
                ZSTD_BCSF_time_startD_fk = System.currentTimeMillis();

                for (int j = 0; j < num_fragments; j++) {

                    File fileCSF = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + zstd_c_EXT);
                    if (fileCSF.exists() && fileCSF.isFile()) {

                        //decompressione CSF part
                        if (debug==1) System.out.println("zstd -d -f " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + zstd_c_EXT);

                        ProcessBuilder builderD = new ProcessBuilder("zstd", "-d", "-f", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + zstd_c_EXT);
                        builderD.redirectErrorStream(true);
                        final Process processD = builderD.start();
                        // Watch the process
                        watch(processD);
                        processD.waitFor();

                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + zstd_c_EXT + " not a file or non-existent file !!!.");
                    }


                    File fileBLOOM = new File(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + zstd_c_EXT);
                    if (fileBLOOM.exists() && fileBLOOM.isFile()) {

                        //decompressione BLOOM part
                        if (debug==1) System.out.println("zstd -d -f " + path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + zstd_c_EXT);

                        ProcessBuilder builderD = new ProcessBuilder("zstd", "-d", "-f", path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + zstd_c_EXT);
                        builderD.redirectErrorStream(true);
                        final Process processD = builderD.start();
                        // Watch the process
                        watch(processD);
                        processD.waitFor();
                    } else {
                        System.out.println(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + zstd_c_EXT + " not a file or non-existent file !!!.");
                    }
                }
                ZSTD_BCSF_time_decompr_fk = System.currentTimeMillis() - ZSTD_BCSF_time_startD_fk;

                // Cleanup Garbage Collector
                System.gc();

                //FIXME: remove ZSTD BCSF compression/decompression  temporary files
                for (int j = 0; j < num_fragments; j++) {
                    remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_CSF_EXT + zstd_c_EXT, debug);
                    remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + BCSF_BLOOM_EXT + zstd_c_EXT, debug);
                }

                //FIXME: remove DSK Blocks temporary files
                for (int j = 0; j < num_fragments; j++) {
                    remove_file(path_in + FS_SEPARATOR + exp + FS_SEPARATOR + exp + approx_Suffix + gap_Suffix + "_" + j + "_DSK" + DSK_EXT, debug);
                }


                //FIXME: write to DPfile size and time information on DP0 scenario !!!.
                //size
                DPfile_size.println(uncompr_size_fk + ", " + bzip2_size_fk + ", " + lz4_size_fk + ", " + zstd_size_fk + ", " + BIC_size_fk + ", " + OptPFOR_size_fk + " - " + Dict_size_DSK + ", " + BCSF_size_build_fk + ", " + BZIP2_BCSF_size_fk + ", " + LZ4_BCSF_size_fk + ", " + ZSTD_BCSF_size_fk);

                //compression time
                DPfile_time.println(bzip2_time_compr_fk + ", " + lz4_time_compr_fk + ", " + zstd_time_compr_fk + ", " + BIC_time_compr_fk + ", " + OptPFOR_time_compr_fk + " - " + BCSF_time_build_fk + ", " + BZIP2_BCSF_time_compr_fk + ", " + LZ4_BCSF_time_compr_fk + ", " + ZSTD_BCSF_time_compr_fk);

                //decompression time
                DPfile_time.println(bzip2_time_decompr_fk + ", " + lz4_time_decompr_fk + ", " + zstd_time_decompr_fk + ", " + BIC_time_decompr_fk + ", " + OptPFOR_time_decompr_fk + " - " + BZIP2_BCSF_time_decompr_fk + ", " + LZ4_BCSF_time_decompr_fk + ", " + ZSTD_BCSF_time_decompr_fk);
            } else {
                //stampa informazioni senza i gaps
                //FIXME: write to DPfile size and time information on DP0 scenario !!!.
                //size
                DPfile_size.println(uncompr_size_fk + ", " + bzip2_size_fk + ", " + lz4_size_fk + ", " + zstd_size_fk + ", " + BIC_size_fk + ", " + OptPFOR_size_fk);

                //compression time
                DPfile_time.println(bzip2_time_compr_fk + ", " + lz4_time_compr_fk + ", " + zstd_time_compr_fk + ", " + BIC_time_compr_fk + ", " + OptPFOR_time_compr_fk);

                //decompression time
                DPfile_time.println(bzip2_time_decompr_fk + ", " + lz4_time_decompr_fk + ", " + zstd_time_decompr_fk + ", " + BIC_time_decompr_fk + ", " + OptPFOR_time_decompr_fk);

            }

        }

        if (DP_scenario > 0) {
            //FIXME: write to DPfile size and time information on DP1 and DP2 scenario !!!.
            //size
            DPfile_size.println(uncompr_size_fk + ", " + bzip2_size_fk + ", " + lz4_size_fk + ", " + zstd_size_fk + ", " + BIC_size_fk + ", " + OptPFOR_size_fk);

            //compression time
            DPfile_time.println(bzip2_time_compr_fk + ", " + lz4_time_compr_fk + ", " + zstd_time_compr_fk + ", " + BIC_time_compr_fk + ", " + OptPFOR_time_compr_fk);

            //decompression time
            DPfile_time.println(bzip2_time_decompr_fk + ", " + lz4_time_decompr_fk + ", " + zstd_time_decompr_fk + ", " + BIC_time_decompr_fk + ", " + OptPFOR_time_decompr_fk);
        }
    }


    private static void Manage_Fb_fragments(String path_in, String exp, int DP_scenario, int approx, PrintWriter DPfile_size, PrintWriter DPfile_time, int debug) throws IOException, InterruptedException {

        //FIXME: Transform Fk to Gk (and Gk0) script after Manage_Fb_fragments per scrivere sui DPfile files
        long Fk_To_Gk_Gk0_start_time = System.currentTimeMillis();
        From_Fk_To_Gk_and_Gk0_exp.execute(path_in + FS_SEPARATOR, exp, approx, debug);
        long Fk_To_Gk_Gk0_time = System.currentTimeMillis() - Fk_To_Gk_Gk0_start_time;
        DPfile_time.println(Fk_To_Gk_Gk0_time);

        //FIXME: apply scenario script to Gk files
        Manage_Fa_fragments(path_in + FS_SEPARATOR, exp, DP_scenario, approx, DPfile_size, DPfile_time, 1, debug);

        //FIXME: Transform from Gk and Gk0 To Fk-Post for post-processing time (all'interno di Manage_Fb_fragments) per scrivere sui files
        long Gk_and_Gk0_To_Fk_start_time = System.currentTimeMillis();
        From_Gk_and_Gk0_To_Fk_exp.execute(path_in + FS_SEPARATOR, exp, approx, debug);
        long Gk_and_Gk0_To_Fk_time = System.currentTimeMillis() - Gk_and_Gk0_To_Fk_start_time;
        DPfile_time.println(Gk_and_Gk0_To_Fk_time);
    }


    private static void Create_DSK_Fragment(String path) throws IOException {
        //input
        BufferedReader Dk_fragment_in = new BufferedReader(new FileReader(path + FASTA_EXT));

        BufferedReader Fk_fragment_in = new BufferedReader(new FileReader(path + DSK_EXT));

        //output
        PrintWriter DSK_fragment_out = new PrintWriter(new FileWriter(path + "_DSK" + DSK_EXT));

        String l_Dk;
        String l_Fk;

        // read the first header line from Dk
        l_Dk = Dk_fragment_in.readLine();

        // until it reaches the end of the file (e.g. of k-meri), I read from the two files and write to the output file
        while ((l_Dk = Dk_fragment_in.readLine()) != null) {
            l_Fk = Fk_fragment_in.readLine();

            int freq_aux = Integer.parseInt(l_Fk);

            DSK_fragment_out.println(l_Dk + " " + freq_aux);
        }

        //close files
        DSK_fragment_out.close();
        Fk_fragment_in.close();
        Dk_fragment_in.close();
    }


    private static void DP_scenario_012(String exp_folder, String exp_name, int k, long DSK_time, int DP_scenario, int b, int debug) throws IOException, InterruptedException {
        //TODO: DP0, DP1 and DP2 scenario

        //FIXME: create DP0 folder and put following file
        String DP_path = exp_folder + FS_SEPARATOR + "k" + k;
        String DP_folder = "DP" + DP_scenario;
        create_folder(DP_path, DP_folder, debug);

        //FIXME: Create and Open DP files depends on: exp_folder, exp_name, k, DP_scenario
        PrintWriter DPscenario_size_out = null;
        try {
            DPscenario_size_out = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP" + DP_scenario + FS_SEPARATOR + "DP" + DP_scenario + "-k" + k + "-size" + DP_EXT));
            if (debug==1) System.out.println("Write DP file to: " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP" + DP_scenario + FS_SEPARATOR + "DP" + DP_scenario + "-k" + k + "-size" + DP_EXT);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        PrintWriter DPscenario_time_out = null;
        try {
            DPscenario_time_out = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP" + DP_scenario + FS_SEPARATOR + "DP" + DP_scenario + "-k" + k + "-time" + DP_EXT));
            if (debug==1) System.out.println("Write DP file to: " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP" + DP_scenario + FS_SEPARATOR + "DP" + DP_scenario + "-k" + k + "-time" + DP_EXT);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        //FIXME: KmerDSKUtility script
        String mode = "";
        switch (DP_scenario) {
            case 0: {
                mode = VERBATIM_Suffix;
                break;
            }

            case 1: {
                mode = LEXICOGR_Suffix;
                break;
            }

            case 2: {
                mode = FREQUENCY_Suffix;
                break;
            }
        }

        long extract_fragment_start_time = System.currentTimeMillis();
        KmerDSKUtility_exp.execute(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR, exp_name + "-DSK-k" + k, DP_scenario, exp_name + "-DSK-k" + k + mode, b, debug);
        long extract_fragment_time = System.currentTimeMillis() - extract_fragment_start_time;

        DPscenario_time_out.println(DSK_time);
        DPscenario_time_out.println(extract_fragment_time);

        //FIXME: compute size of uncompressed, compressed and an other uncompressed Dk fragments files
        Manage_Dk_fragments(exp_folder + FS_SEPARATOR + "k" + k, exp_name + "-DSK-k" + k + mode, DP_scenario, DPscenario_size_out, DPscenario_time_out, debug);

        //FIXME: compute size of uncompressed, compressed and an other uncompressed Fk fragments files (no Gaps)
        Manage_Fa_fragments(exp_folder + FS_SEPARATOR + "k" + k, exp_name + "-DSK-k" + k + mode, DP_scenario, 0, DPscenario_size_out, DPscenario_time_out, 0, debug);

        //TODO: create DSK Dictionary with Gaps before managing Fb_fragments and storing into parent directory (i.e. where is StaphAU-DSK-k4-verbatim.txt)
        // NOTE: We actually use only for DP0 scenario, but we can use also for DP1 and DP2 scenarios.
        if (DP_scenario == 0) {
            Create_DSK_Dictionary_with_Gaps(exp_folder, exp_name, k, mode, debug);
            //TODO: no longer used because BCSF by construction works with non-negative occurrences.
        }

        //FIXME: compute size of uncompressed, compressed and an other uncompressed Gk fragments files (with Gaps)
        Manage_Fb_fragments(exp_folder + FS_SEPARATOR + "k" + k, exp_name + "-DSK-k" + k + mode, DP_scenario, 0, DPscenario_size_out, DPscenario_time_out, debug);

        //FIXME: close DP files
        DPscenario_time_out.close();
        DPscenario_size_out.close();
    }


    private static void Create_DSK_Dictionary_with_Gaps(String exp_folder, String exp_name, int k, String mode, int debug) throws IOException {

        String pathfile_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-DSK-k" + k + mode + DSK_EXT;
        if (debug==1) System.out.println("pathfile_in = " + pathfile_in);

        String pathfile_out = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp_name + "-DSK-k" + k + mode + "-Gk" + DSK_EXT;
        if (debug==1) System.out.println("pathfile_out = " + pathfile_out);

        // Input reading of the DSK verbatim dictionary for the construction of the Gk version
        BufferedReader dict_Dk = new BufferedReader(new FileReader(pathfile_in));

        // Preparation of DSK dictionary output file with gaps
        PrintWriter dict_Gk = new PrintWriter(new FileWriter(pathfile_out));

        // Support string for reading the lines of the file
        String line;

        // Reading from input file and contextual construction of gaps file
        String[] aux;
        String kmer, freq;
        int aux_freq, gap = 0;

        // The first time, subtract zero from the first frequency, i.e. assign freq[0] to the first gap, while the other times,
        // subtract freq[i] from partial updated to freq[i-1] .
        int partial = 0;

        while ((line = dict_Dk.readLine()) != null) {
            aux = line.split("\\s+");

            kmer = aux[0];
            freq = aux[1];

            aux_freq = Integer.parseInt(aux[1]);
            gap = aux_freq - partial;

            //Storage (kmer, gap) into dict_Gk file
            dict_Gk.println(kmer + " " + gap);

            //update partial
            partial = aux_freq;
        }

        //closing input and output files
        dict_Gk.close();
        dict_Dk.close();
    }


    private static void DP_scenario_3(String exp_folder, String exp, int DP_scenario, int k, long pre_proc_dsk_time, long pre_proc_ess_time, int b, int debug) throws IOException, InterruptedException {


        //FIXME: create DP0 folder and put following file
        String DP_path = exp_folder + FS_SEPARATOR + "k" + k;
        String DP_folder = "DP3";
        create_folder(DP_path, DP_folder, debug);

        // Compute DP3 Case on CD-NRAM Scenario
        //FIXME: Create and Open DP files depends on: exp_folder, exp_name, k
        PrintWriter DPscenario_size_out = null;
        try {
            DPscenario_size_out = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP3" + FS_SEPARATOR + "DP3" + "-k" + k + "-size" + "-CD_NRAM" + DP_EXT));
            if (debug==1) System.out.println("Write DP file to: " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP3" + FS_SEPARATOR + "DP3" + "-k" + k + "-size" + "-CD_NRAM" + DP_EXT);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        PrintWriter DPscenario_time_out = null;
        try {
            DPscenario_time_out = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP3" + FS_SEPARATOR + "DP3" + "-k" + k + "-time" + "-CD_NRAM" + DP_EXT));
            if (debug==1) System.out.println("Write DP file to: " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP3" + FS_SEPARATOR + "DP3" + "-k" + k + "-time" + "-CD_NRAM" + DP_EXT);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }


        if (Approx_Solution == 0) { // We are in Exact Solution Context - CD-NRAM

            System.out.println("DP3 - We are in Exact Solution Context");

            // Compute DP3 Case on CD-NRAM Scenario
            //FIXME: Compute F'k bijection from SPSS S and DSK sorted output via BijectionFkFromDSKToSPSSKmer script
            long bijection_fk_fragment_start_time = System.currentTimeMillis();
            BijectionFkFromDSKToSPSSkmer_exact.execute(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR, exp + "-ESS-k" + k, k, exp + "-DSK-k" + k + VERBATIM_Suffix, exp + "-DSK-k" + k + BIJECTION_Suffix, b, debug);
            long bijection_fk_fragment_time = System.currentTimeMillis() - bijection_fk_fragment_start_time;

            DPscenario_time_out.println(pre_proc_dsk_time);     // first row in DP3-time
            DPscenario_time_out.println(pre_proc_ess_time);     // second row in DP3-time
            DPscenario_time_out.println(bijection_fk_fragment_time);    // (F'k Bijection) third row in DP3-time


            //FIXME: compute size of uncompressed, compressed and an other uncompressed F'k fragments files (no Gaps)
            //4-th and 5-th rows in DP3-time
            Manage_Fa_fragments(exp_folder + FS_SEPARATOR + "k" + k, exp + "-DSK-k" + k + BIJECTION_Suffix, DP_scenario, 0, DPscenario_size_out, DPscenario_time_out, 0, debug);

            //FIXME: compute size of uncompressed, compressed and an other uncompressed G'k fragments files (with Gaps)
            //6-th row in DP3-time - Transform Fk to Gk (and Gk0)
            //7-th and 8-th rows in DP3-time
            Manage_Fb_fragments(exp_folder + FS_SEPARATOR + "k" + k, exp + "-DSK-k" + k + BIJECTION_Suffix, DP_scenario, 0, DPscenario_size_out, DPscenario_time_out, debug);

            //FIXME: Compute D'k bijection from SPSS S and (F'k bijection) via From_S_and_Fk_rearranged_To_Dk_fragment script
            long bijection_dk_fragment_start_time = System.currentTimeMillis();
            From_S_and_Fk_rearranged_To_Dk_exact.execute(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR, exp + "-ESS-k" + k, k, exp + "-DSK-k" + k + BIJECTION_Suffix, b, debug);
            long bijection_dk_fragment_time = System.currentTimeMillis() - bijection_dk_fragment_start_time;
            //10-th row in DP3-time - Compute D'k bijection from SPSS S and (F'k bijection)
            DPscenario_time_out.println(bijection_dk_fragment_time);

            //FIXME: compute size of uncompressed, compressed and an other uncompressed D'k fragments files
            //11-th and 12-th rows in DP3-time
            Manage_Sk_fragments(exp_folder + FS_SEPARATOR + "k" + k, exp + "-ESS-k" + k, DP_scenario, 1, DPscenario_size_out, DPscenario_time_out, debug);

            // Cleanup Garbage Collector
            System.gc();

        } else { // We are in Approximate Solution Context - CD-NRAM    //TODO: Approx_Solution == 1

            System.out.println("DP3 - We are in Approximate Solution Context");
            DPscenario_size_out.println("Approximate Solution");
            DPscenario_time_out.println("Approximate Solution");

            // TODO: Leggi Num_Blocks dal file .nfo
            int num_blocks = 1;     // set to default value
            BufferedReader S_nfo = null;
            try {
                S_nfo = new BufferedReader(new FileReader(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-ESS-k" + k + nfo_EXT));
                if (debug==1) System.out.println("ESS Blocks info .nfo file: " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-ESS-k" + k + nfo_EXT);

                String l = S_nfo.readLine();
                num_blocks = Integer.parseInt(l);

                S_nfo.close();
                if (debug==1) System.out.println("num_blocks = " + num_blocks);

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }


            // TODO: Compute DSK for every Block od Dataset DSK for Dictionary extraction
            long pre_proc_DSK_start = System.currentTimeMillis();
            extract_DSK_blocks_dictionaries(exp_folder, exp, k, num_blocks, debug);
            long pre_proc_DSK_time = System.currentTimeMillis() - pre_proc_DSK_start;


            // FIXME: Compute Bijection F'k for every Block ESS and Block DSK for construction Dk anf Fk
            long bijection_fk_fragment_start_time = System.currentTimeMillis();
            BijectionFkFromDSKToSPSSkmer_approx.execute(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR, exp, exp + "-ESS-k" + k, k, exp + "-DSK-k" + k + VERBATIM_Suffix, exp + "-DSK-k" + k + BIJECTION_Suffix, num_blocks, debug);
            long bijection_fk_fragment_time = System.currentTimeMillis() - bijection_fk_fragment_start_time;


            //FIXME: Compute D'k bijection from SPSS S and (F'k bijection) via From_S_and_Fk_rearranged_To_Dk_fragment script
            long bijection_dk_fragment_start_time = System.currentTimeMillis();
            From_S_and_Fk_rearranged_To_Dk_approx.execute(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR, exp, exp + "-ESS-k" + k, k, exp + "-DSK-k" + k + BIJECTION_Suffix, num_blocks, debug);
            long bijection_dk_fragment_time = System.currentTimeMillis() - bijection_dk_fragment_start_time;


            //FIXME: Compute DSK Approximation
            long bijection_createDictionary_start_time = System.currentTimeMillis();
            Bijection_Create_Dictionary_approx.execute(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR, exp, exp + "-ESS-k" + k, k, exp + "-DSK-k" + k + BIJECTION_Suffix, b, num_blocks, debug);
            long bijection_createDictionary_time = System.currentTimeMillis() - bijection_createDictionary_start_time;

            //FIXME: compute size of uncompressed, compressed and an other uncompressed F'k fragments files (no Gaps)
            //4-th and 5-th rows in DP3-time
            Manage_Fa_fragments(exp_folder + FS_SEPARATOR + "k" + k, exp + "-DSK-k" + k + BIJECTION_Suffix, DP_scenario, 1, DPscenario_size_out, DPscenario_time_out, 0, debug);

            //FIXME: compute size of uncompressed, compressed and an other uncompressed G'k fragments files (with Gaps)
            //6-th row in DP3-time - Transform Fk to Gk (and Gk0)
            //7-th and 8-th rows in DP3-time
            Manage_Fb_fragments(exp_folder + FS_SEPARATOR + "k" + k, exp + "-DSK-k" + k + BIJECTION_Suffix, DP_scenario, 1, DPscenario_size_out, DPscenario_time_out, debug);

            //FIXME: compute size of uncompressed, compressed and an other uncompressed D'k fragments files
            //11-th and 12-th rows in DP3-time
            Manage_Sk_fragments(exp_folder + FS_SEPARATOR + "k" + k, exp + "-ESS-k" + k, DP_scenario, num_blocks, DPscenario_size_out, DPscenario_time_out, debug);

            DPscenario_time_out.println(pre_proc_DSK_time);                 // first row in DP3-time
            DPscenario_time_out.println(pre_proc_ess_time);                 // second row in DP3-time
            DPscenario_time_out.println(bijection_fk_fragment_time);        // (F'k Bijection) third row in DP3-time
            DPscenario_time_out.println(bijection_dk_fragment_time);        // (D'k bijection) from SPSS S and (F'k bijection)
            DPscenario_time_out.println(bijection_createDictionary_time);   // Bijection Create Dictionary from D_i and F_i

            // Cleanup Garbage Collector
            System.gc();
        }

        //FIXME: close DP files
        DPscenario_time_out.close();
        DPscenario_size_out.close();


        if (Approx_Solution == 0) { // We are in Exact Solution Context - SD-RAM

            // Compute DP3 scenario with SD-RAM Scenario
            //FIXME: Create and Open DP files depends on: exp_folder, exp_name, k for Succinct on Disk scenario
            PrintWriter DPscenarioSD_size_out = null;
            try {
                DPscenarioSD_size_out = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP3" + FS_SEPARATOR + "DP3" + "-k" + k + "-size" + "-SD_RAM" + DP_EXT));
                if (debug==1) System.out.println("Write DP file to: " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP3" + FS_SEPARATOR + "DP3" + "-k" + k + "-size" + "-SD_RAM" + DP_EXT);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            PrintWriter DPscenarioSD_time_out = null;
            try {
                DPscenarioSD_time_out = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP3" + FS_SEPARATOR + "DP3" + "-k" + k + "-time" + "-SD_RAM" + DP_EXT));
                if (debug==1) System.out.println("Write DP file to: " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP3" + FS_SEPARATOR + "DP3" + "-k" + k + "-time" + "-SD_RAM" + DP_EXT);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }


            long bijection_FMindex_Of_S_start_time = System.currentTimeMillis();
            Compute_FMindex_Of_S(exp_folder, exp, BIJECTION_Suffix, k, 1, DPscenarioSD_size_out, debug);
            long bijection_FMindex_Of_S_time = System.currentTimeMillis() - bijection_FMindex_Of_S_start_time;
            //1_th row on DPscenarioSD_time_out
            DPscenarioSD_time_out.println(bijection_FMindex_Of_S_time);

            // move FM-index of S (.fm9) into exp-Bijection directory. We don't rename with -Bijection suffix
            move_file_to_path(exp + "-ESS-k" + k + FASTA_EXT + FMind_c_EXT, exp_folder + FS_SEPARATOR + "k" + k, exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix, debug);


            //TODO: Extract from FM-index(S)
            long bijection_Extract_from_FMindex_Of_S_start_time = System.currentTimeMillis();
            Extract_from_FMindex_Of_S(exp_folder, exp, BIJECTION_Suffix, k, 1, DPscenarioSD_size_out, debug);
            long bijection_Extract_from_FMindex_Of_S_time = System.currentTimeMillis() - bijection_Extract_from_FMindex_Of_S_start_time;
            //2_th row on DPscenarioSD_time_out
            DPscenarioSD_time_out.println(bijection_Extract_from_FMindex_Of_S_time);

            //TODO: Recovery Dk from S (in .fasta format)
            long bijection_Recovery_Dk_from_S_start_time = System.currentTimeMillis();
            Recovery_Dk_from_S(exp_folder, exp, k, 1, debug);
            long bijection_Recovery_Dk_from_S_time = System.currentTimeMillis() - bijection_Recovery_Dk_from_S_start_time;
            //3_th row on DPscenarioSD_time_out
            DPscenarioSD_time_out.println(bijection_Recovery_Dk_from_S_time);


            //TODO: prepare Dk for Fk recovering
            long bijection_Prepare_Dk_for_Recovering_Fk_start_time = System.currentTimeMillis();
            bijection_Prepare_Dk_for_Recovering_Fk(exp_folder, exp, k, 1, debug);
            long bijection_Prepare_Dk_for_Recovering_Fk_time = System.currentTimeMillis() - bijection_Prepare_Dk_for_Recovering_Fk_start_time;
            //4_th row on DPscenarioSD_time_out
            DPscenarioSD_time_out.println(bijection_Prepare_Dk_for_Recovering_Fk_time);


            //TODO: Recovery Fk from Dk on BCSF data structure (in .txt format) without gaps
            for (int indexGap = 0; indexGap < 1; indexGap++) {
                long bijection_Recovery_Fk_from_Dk_and_BCSF_start_time = System.currentTimeMillis();
                Recovery_Fk_from_Dk_and_BCSF(exp_folder, exp, k, 1, indexGap, debug);
                long bijection_Recovery_Fk_from_Dk_and_BCSF_time = System.currentTimeMillis() - bijection_Recovery_Fk_from_Dk_and_BCSF_start_time;
                //5_th row on DPscenarioSD_time_out
                DPscenarioSD_time_out.println(bijection_Recovery_Fk_from_Dk_and_BCSF_time);
            }

            //FIXME: close DP files for Succinct on Disk scenario
            DPscenarioSD_time_out.close();
            DPscenarioSD_size_out.close();

            // Cleanup Garbage Collector
            System.gc();

        } else {  // We are in Approximate Solution Context - SD-RAM    //TODO: Approx_Solution == 1

            // TODO: Leggi Num_Blocks dal file .nfo
            int num_blocks = 1;     // set to default value
            BufferedReader S_nfo = null;
            try {
                S_nfo = new BufferedReader(new FileReader(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-ESS-k" + k + nfo_EXT));
                if (debug==1) System.out.println("ESS Blocks info .nfo file: " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-ESS-k" + k + nfo_EXT);

                String l = S_nfo.readLine();
                num_blocks = Integer.parseInt(l);

                S_nfo.close();
                if (debug==1) System.out.println("num_blocks = " + num_blocks);

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            // Compute DP3 scenario with SD-RAM Scenario
            //FIXME: Create and Open DP files depends on: exp_folder, exp_name, k for Succinct on Disk scenario
            PrintWriter DPscenarioSD_size_out = null;
            try {
                DPscenarioSD_size_out = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP3" + FS_SEPARATOR + "DP3" + "-k" + k + "-size" + "-SD_RAM" + DP_EXT));
                if (debug==1) System.out.println("Write DP file to: " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP3" + FS_SEPARATOR + "DP3" + "-k" + k + "-size" + "-SD_RAM" + DP_EXT);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            PrintWriter DPscenarioSD_time_out = null;
            try {
                DPscenarioSD_time_out = new PrintWriter(new FileWriter(exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP3" + FS_SEPARATOR + "DP3" + "-k" + k + "-time" + "-SD_RAM" + DP_EXT));
                if (debug==1) System.out.println("Write DP file to: " + exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + "DP3" + FS_SEPARATOR + "DP3" + "-k" + k + "-time" + "-SD_RAM" + DP_EXT);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            DPscenarioSD_size_out.println("Approximate Solution");
            DPscenarioSD_time_out.println("Approximate Solution");


            long bijection_FMindex_Of_S_start_time = System.currentTimeMillis();
            Compute_FMindex_Of_S(exp_folder, exp, BIJECTION_Suffix, k, num_blocks, DPscenarioSD_size_out, debug);
            long bijection_FMindex_Of_S_time = System.currentTimeMillis() - bijection_FMindex_Of_S_start_time;
            //1_th row on DPscenarioSD_time_out
            DPscenarioSD_time_out.println(bijection_FMindex_Of_S_time);

            // move FM-index of S_i (.fm9) into exp-Bijection directory. We don't rename with -Bijection suffix
            for (int i = 0; i < num_blocks; i++) {
                move_file_to_path(exp + "-ESS-k" + k + "_Block_" + i + FASTA_EXT + FMind_c_EXT, exp_folder + FS_SEPARATOR + "k" + k, exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix, debug);
            }


            //TODO: Extract from FM-index(S_i)
            long bijection_Extract_from_FMindex_Of_S_start_time = System.currentTimeMillis();
            Extract_from_FMindex_Of_S(exp_folder, exp, BIJECTION_Suffix, k, num_blocks, DPscenarioSD_size_out, debug);
            long bijection_Extract_from_FMindex_Of_S_time = System.currentTimeMillis() - bijection_Extract_from_FMindex_Of_S_start_time;
            //2_th row on DPscenarioSD_time_out
            DPscenarioSD_time_out.println(bijection_Extract_from_FMindex_Of_S_time);


            //TODO: Recovery Dk from S_i (in .fasta format)
            long bijection_Recovery_Dk_from_S_start_time = System.currentTimeMillis();
            Recovery_Dk_from_S(exp_folder, exp, k, num_blocks, debug);
            long bijection_Recovery_Dk_from_S_time = System.currentTimeMillis() - bijection_Recovery_Dk_from_S_start_time;
            //3_th row on DPscenarioSD_time_out
            DPscenarioSD_time_out.println(bijection_Recovery_Dk_from_S_time);


            //TODO: prepare Dk for Fk recovering
            long bijection_Prepare_Dk_for_Recovering_Fk_start_time = System.currentTimeMillis();
            bijection_Prepare_Dk_for_Recovering_Fk(exp_folder, exp, k, num_blocks, debug);
            long bijection_Prepare_Dk_for_Recovering_Fk_time = System.currentTimeMillis() - bijection_Prepare_Dk_for_Recovering_Fk_start_time;
            //4_th row on DPscenarioSD_time_out
            DPscenarioSD_time_out.println(bijection_Prepare_Dk_for_Recovering_Fk_time);


            //TODO: Recovery Fk from Dk on BCSF data structure (in .txt format) without gaps
            for (int indexGap = 0; indexGap < 1; indexGap++) {  //lo faccio eseguire solo per il caso senza gaps, infatti nel caso con gaps restituirebbe dati errati, per ratio delle occorrenze di BCSF
                long bijection_Recovery_Fk_from_Dk_and_BCSF_start_time = System.currentTimeMillis();
                Recovery_Fk_from_Dk_and_BCSF(exp_folder, exp, k, num_blocks, indexGap, debug);
                long bijection_Recovery_Fk_from_Dk_and_BCSF_time = System.currentTimeMillis() - bijection_Recovery_Fk_from_Dk_and_BCSF_start_time;
                //5_th row on DPscenarioSD_time_out
                DPscenarioSD_time_out.println(bijection_Recovery_Fk_from_Dk_and_BCSF_time);
            }

            //FIXME: close DP files for Succinct on Disk scenario
            DPscenarioSD_time_out.close();
            DPscenarioSD_size_out.close();

            // Cleanup Garbage Collector
            System.gc();
        }
    }


    private static void extract_DSK_blocks_dictionaries(String exp_folder, String exp, int k, int num_blocks, int debug) throws IOException, InterruptedException {

        for (int i = 0; i < num_blocks; i++) {

            //FIXME: define block on which apply DSK software
            String block_fasta = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "_B_" + i + FASTA_EXT;
            String block_h5 = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "_B_" + i + ".h5";
            String block_dsk = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + "-verbatim_" + i + DSK_EXT;

            // create .h5 file
            if (debug==1) System.out.println("./essAuxDsk -file " + block_fasta + " -kmer-size " + k + " -abundance-min 1 -out " + block_h5 + " -verbose 0 -max-memory 5000");
            ProcessBuilder builder = new ProcessBuilder("./essAuxDsk", "-file", block_fasta, "-kmer-size", String.valueOf(k), "-abundance-min", String.valueOf(1), "-out", block_h5, "-verbose", String.valueOf(0), "-max-memory", String.valueOf(5000));
            builder.redirectErrorStream(true);
            final Process process = builder.start();
            // Watch the process
            watch(process);
            process.waitFor();

            // conversion .h5 file to .txt file
            if (debug==1) System.out.println("./essAuxDsk2ascii -file " + block_h5 + " -out " + block_dsk + " -verbose 0");
            builder = new ProcessBuilder("./essAuxDsk2ascii", "-file", block_h5, "-out", block_dsk, "-verbose", String.valueOf(0));
            builder.redirectErrorStream(true);
            final Process process2 = builder.start();
            // Watch the process
            watch(process2);
            process2.waitFor();

            // remove .h5 file
            if (debug==1) System.out.println("rm -f " + block_h5);
            remove_file(block_h5, debug);

            // FIXME: remove block.fast? Will I need it later?
            //  Yes, I have to take the border lines of the blocks for the k-mer lost due to the cut.
            //  However, I save this information when creating blocks, so here I can remove blocks.
            remove_file(block_fasta, debug);

            // Clean Garbage Collector
            System.gc();
        }
    }


    private static void bijection_Prepare_Dk_for_Recovering_Fk(String exp_folder, String exp, int k, int num_blocks, int debug) throws IOException, InterruptedException {

        for (int i = 0; i < num_blocks; i++) {

            String pathfile_Dk_in = null;
            String pathfile_Dk_out = null;

            if (num_blocks == 1) {

                pathfile_Dk_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + "-Recovery" + FASTA_EXT;
                if (debug==1) System.out.println("pathfile_Dk_in = " + pathfile_Dk_in);

                pathfile_Dk_out = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + "-Recovery" + BCSF_backup + FASTA_EXT;
                if (debug==1) System.out.println("pathfile_Dk_out = " + pathfile_Dk_out);
            } else {

                pathfile_Dk_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + "-Recovery" + "_" + i + FASTA_EXT;
                if (debug==1) System.out.println("pathfile_Dk_in = " + pathfile_Dk_in);

                pathfile_Dk_out = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + "-Recovery" + "_" + i + BCSF_backup + FASTA_EXT;
                if (debug==1) System.out.println("pathfile_Dk_out = " + pathfile_Dk_out);
            }

            //inizialize output file
            PrintWriter Dk_out = null;
            try {
                Dk_out = new PrintWriter(new FileWriter(pathfile_Dk_out));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            // Write header into Dk_out output file
            Dk_out.println(">");

            // Initialize input file
            BufferedReader Dk_in = null;
            try {
                Dk_in = new BufferedReader(new FileReader(pathfile_Dk_in));

                //Scan S for storage k-mers dictionary
                String l;
                while ((l = Dk_in.readLine()) != null) {
                    if (!l.startsWith(">")) {

                        if (l.length() == k) {
                            Dk_out.println(l + " 0");
                        }
                    }
                }
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }

            //close input and output files
            Dk_out.close();
            Dk_in.close();
        }
    }


    private static void Recovery_Fk_from_Dk_and_BCSF(String exp_folder, String exp, int k, int num_blocks, int gap, int debug) throws IOException, InterruptedException {

        //gap Info Management
        String gap_Suffix = "";
        if (gap == 1) {
            gap_Suffix = "-Gk";
        }

        for (int i = 0; i < num_blocks; i++) {

            String pathfile_Dk_in = null;
            String pathfile_BCSF_csf_in = null;
            String pathfile_BCSF_bloom_in = null;
            String pathfile_Fk_out = null;

            if (num_blocks == 1) {

                pathfile_Dk_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + "-Recovery" + BCSF_backup + FASTA_EXT;
                if (debug==1) System.out.println("pathfile_Dk_in = " + pathfile_Dk_in);

                pathfile_BCSF_csf_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + VERBATIM_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + VERBATIM_Suffix + gap_Suffix + BCSF_CSF_EXT;
                if (debug==1) System.out.println("pathfile_BCSF_csf_in = " + pathfile_BCSF_csf_in);

                pathfile_BCSF_bloom_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + VERBATIM_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + VERBATIM_Suffix + gap_Suffix + BCSF_BLOOM_EXT;
                if (debug==1) System.out.println("pathfile_BCSF_bloom_in = " + pathfile_BCSF_bloom_in);

                pathfile_Fk_out = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + "-Recovery" + gap_Suffix + DSK_EXT;
                if (debug==1) System.out.println("pathfile_Fk_out = " + pathfile_Fk_out);

            } else {

                pathfile_Dk_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + "-Recovery" + "_" + i + BCSF_backup + FASTA_EXT;
                if (debug==1) System.out.println("pathfile_Dk_in = " + pathfile_Dk_in);

                pathfile_BCSF_csf_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + VERBATIM_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + VERBATIM_Suffix + gap_Suffix + "_" + i + BCSF_CSF_EXT;
                if (debug==1) System.out.println("pathfile_BCSF_csf_in = " + pathfile_BCSF_csf_in);

                pathfile_BCSF_bloom_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + VERBATIM_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + VERBATIM_Suffix + gap_Suffix + "_" + i + BCSF_BLOOM_EXT;
                if (debug==1) System.out.println("pathfile_BCSF_bloom_in = " + pathfile_BCSF_bloom_in);

                pathfile_Fk_out = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + "-Recovery" + gap_Suffix + "_" + i + DSK_EXT;
                if (debug==1) System.out.println("pathfile_Fk_out = " + pathfile_Fk_out);
            }

            //FIXME: BCSF query to get Fk
            String cmd;

            //TODO: distinguere i casi in cui ho o meno il bloom filter associato
            File fileB = new File(pathfile_BCSF_bloom_in);
            if (fileB.exists() && fileB.isFile()) {
                //bloom filter presente
                if (debug==1) System.out.println("./bcsf_query " + pathfile_BCSF_bloom_in + " " + pathfile_BCSF_csf_in + " < " + pathfile_Dk_in + " > " + pathfile_Fk_out);
                ProcessBuilder builderQ = new ProcessBuilder("./bcsf_query", pathfile_BCSF_bloom_in, pathfile_BCSF_csf_in)
                        .redirectInput(new File(pathfile_Dk_in))
                        .redirectOutput(new File(pathfile_Fk_out));
                builderQ.redirectErrorStream(true);
                final Process processQ = builderQ.start();
                // Watch the process
                watch(processQ);
                processQ.waitFor();
            } else {
                //bloom filter assente, sostituire con "-" placeholder
                if (debug==1) System.out.println("./bcsf_query " + "- " + pathfile_BCSF_csf_in + " < " + pathfile_Dk_in + " > " + pathfile_Fk_out);
                ProcessBuilder builderQ = new ProcessBuilder("./bcsf_query", "-", pathfile_BCSF_csf_in)
                        .redirectInput(new File(pathfile_Dk_in))
                        .redirectOutput(new File(pathfile_Fk_out));
                builderQ.redirectErrorStream(true);
                final Process processQ = builderQ.start();
                // Watch the process
                watch(processQ);
                processQ.waitFor();
            }

            //TODO: remove temporary files, but external of the function, bacause are needed for recovery Fk
            //1
            //remove_file(pathfile_Dk_in);
            //2
            //String pathfile_Dk_Recovery = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + "-Recovery" + "_" + i + FASTA_EXT;
            //remove_file(pathfile_Dk_Recovery);
            //3
            //remove_file(pathfile_Fk_out);
        }
    }


    private static void Recovery_Dk_from_S(String exp_folder, String exp, int k, int num_blocks, int debug) throws IOException {

        for (int i = 0; i < num_blocks; i++) {

            String pathfile_ESS_in = null;
            String pathfile_Dk_out = null;

            if (num_blocks == 1) {

                pathfile_ESS_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-ESS-k" + k + FASTA_EXT;
                if (debug==1) System.out.println("pathfile_ESS_in = " + pathfile_ESS_in);

                pathfile_Dk_out = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + "-Recovery" + FASTA_EXT;
                if (debug==1) System.out.println("pathfile_Dk_out = " + pathfile_Dk_out);
            } else {

                pathfile_ESS_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-ESS-k" + k + "_Block_" + i + FASTA_EXT;
                if (debug==1) System.out.println("pathfile_ESS_in = " + pathfile_ESS_in);

                pathfile_Dk_out = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + "-Recovery" + "_" + i + FASTA_EXT;
                if (debug==1) System.out.println("pathfile_Dk_out = " + pathfile_Dk_out);
            }

            //initialize output file
            PrintWriter Dk_out = null;
            try {
                Dk_out = new PrintWriter(new FileWriter(pathfile_Dk_out));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            // Write header into Dk_out output file
            Dk_out.println(">");

            // inizialize input file
            BufferedReader ESS = null;
            try {
                ESS = new BufferedReader(new FileReader(pathfile_ESS_in));

                //Scan S for storage k-mers Dictionary
                String l;

                //k-mer counter
                int i_kmer = 0;

                while ((l = ESS.readLine()) != null) {
                    // Skipping lines beginning with the '>' character
                    if (!l.startsWith(">")) {
                        // K-mer extraction from SPSS sequence
                        for (int j = 0; j < l.length() - k + 1; j++) {

                            // k-mer selection
                            String kmer = l.substring(j, j + k);

                            //update i_kmer counter
                            i_kmer++;

                            //FIXME: For each x of SPSS S, consider its canonical with respect to the DSK relation order, i.e. A<C<T<G .
                            //  see https://github.com/GATB/dsk#kmers-and-their-reverse-complements
                            Dk_out.println(StringUtils.Canonical_DSK_RelOrd(kmer));
                        }
                    }
                }

            /*
            //TODO: apply here the k-mers file correction for correct operation of the SPRING specialized compressor
            if (k > 7) {
                if ((i_kmer + 2) % 4 != 0) {

                    int i_aux = 1;
                    Dk_out.println(">");

                    while ((i_kmer + 2 + i_aux) % 4 != 0) {
                        Dk_out.println("");

                        //update i_aux
                        i_aux++;
                    }
                }
            }
            */

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }

            // closing files
            ESS.close();
            Dk_out.close();
        }
    }


    private static void Extract_from_FMindex_Of_S(String exp_folder, String exp, String bijection_suffix, int k, int num_blocks, PrintWriter dPscenarioSD_size_out, int debug) throws IOException, InterruptedException {

        for (int i = 0; i < num_blocks; i++) {

            String pathfile_S_fasta_in = null;
            String pathfile_S_fm9_in = null;

            if (num_blocks == 1) {
                pathfile_S_fasta_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-ESS-k" + k + FASTA_EXT;
                if (debug==1) System.out.println("pathfile_S_fasta_in = " + pathfile_S_fasta_in);

                pathfile_S_fm9_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-ESS-k" + k + FASTA_EXT + FMind_c_EXT;
                if (debug==1) System.out.println("pathfile_S_fm9_in = " + pathfile_S_fm9_in);
            } else {
                pathfile_S_fasta_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-ESS-k" + k + "_Block_" + i + FASTA_EXT;
                if (debug==1) System.out.println("pathfile_S_fasta_in = " + pathfile_S_fasta_in);

                pathfile_S_fm9_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-ESS-k" + k + "_Block_" + i + FASTA_EXT + FMind_c_EXT;
                if (debug==1) System.out.println("pathfile_S_fm9_in = " + pathfile_S_fm9_in);
            }

            // decompression of FM-index of S
            if (debug==1) System.out.println("./fm-extract " + pathfile_S_fm9_in + " " + pathfile_S_fasta_in);

            ProcessBuilder builderD = new ProcessBuilder("./fm-extract", pathfile_S_fm9_in, pathfile_S_fasta_in);
            builderD.redirectErrorStream(true);
            final Process processD = builderD.start();
            // Watch the process
            watch(processD);
            processD.waitFor();

            //TODO: remove temporary decompression FM-index .dec file
            remove_file(pathfile_S_fasta_in + FMind_d_EXT, debug);
        }
    }


    private static void Compute_BCSF_Of_DSK_Verbatim(String exp_folder, String exp, int k, PrintWriter DP_SD_RAM_size) throws IOException, InterruptedException {

        // Compute BCSF with verbatim version of DSK
        String pathfile_DSK_verbatim_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + VERBATIM_Suffix + DSK_EXT;
        //System.out.println("pathfile_DSK_verbatim_in = " + pathfile_DSK_verbatim_in);

        String pathfile_BCSF_out = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix + FS_SEPARATOR + exp + "-DSK-k" + k + BIJECTION_Suffix;

        String cmd;
        long uncompr_DSK_size_fk = 0, BCSF_size_fk = 0;
        long BCSF_size_build_CSF_fk = 0, BCSF_size_build_BLOOM_fk = 0, BCSF_size_build_fk = 0;

        File file = new File(pathfile_DSK_verbatim_in);

        if (file.exists() && file.isFile()) {
            uncompr_DSK_size_fk += file.length();
            DP_SD_RAM_size.print(uncompr_DSK_size_fk + ", ");

            //FIXME: BCSF(Fk) Compression (build BCSF)

            cmd = "python3 locom.py build -i " + pathfile_DSK_verbatim_in + " -o " + pathfile_BCSF_out;
            System.out.println(cmd);

            ProcessBuilder builderC = new ProcessBuilder("python3", "locom.py", "build", "-i", pathfile_DSK_verbatim_in, "-o", pathfile_BCSF_out);
            builderC.redirectErrorStream(true);
            final Process processC = builderC.start();
            // Watch the process
            watch(processC);
            processC.waitFor();

            //update BCSF_size_CSF_fk [bytes] with size of current txt fragment
            file = new File(pathfile_BCSF_out + BCSF_CSF_EXT);
            if (file.exists() && file.isFile()) {
                BCSF_size_build_CSF_fk += file.length();
            } else {
                System.out.println(pathfile_BCSF_out + BCSF_CSF_EXT + " not a file or non-existent file !!!.");
            }

            //update BCSF_size_BLOOM_fk [bytes] with size of current txt fragment
            file = new File(pathfile_BCSF_out + BCSF_BLOOM_EXT);
            if (file.exists() && file.isFile()) {
                BCSF_size_build_BLOOM_fk += file.length();
            } else {
                System.out.println(pathfile_BCSF_out + BCSF_BLOOM_EXT + " not a file or non-existent file !!!.");
            }
        }
        BCSF_size_build_fk += BCSF_size_build_CSF_fk + BCSF_size_build_BLOOM_fk;
        DP_SD_RAM_size.print(BCSF_size_build_fk + ", ");

    }


    private static void Compute_FMindex_Of_S(String exp_folder, String exp, String mode, int k, int num_blocks, PrintWriter DP_SD_RAM_size, int debug) throws IOException, InterruptedException {

        //FIXME: refactoring via num_blocks: exact Vs. approximate solution
        String cmd;
        long uncompr_ESS_size_dk = 0, FM_size_dk = 0;

        for (int i = 0; i < num_blocks; i++) {

            String pathfile_ESS_in = null;

            if (num_blocks == 1) {
                pathfile_ESS_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-ESS-k" + k + FASTA_EXT;
            } else {
                pathfile_ESS_in = exp_folder + FS_SEPARATOR + "k" + k + FS_SEPARATOR + exp + "-ESS-k" + k + "_Block_" + i + FASTA_EXT;
            }
            if (debug==1) System.out.println("pathfile_ESS_in = " + pathfile_ESS_in);

            // compute size of ESS S or Block of S in FASTA format
            File file = new File(pathfile_ESS_in);
            if (file.exists() && file.isFile()) {
                uncompr_ESS_size_dk += file.length();

                //FIXME: FM-index(Dk) Compression (build FM-index)

                if (debug==1) System.out.println("./fm-build " + pathfile_ESS_in);

                ProcessBuilder builderC = new ProcessBuilder("./fm-build", pathfile_ESS_in);
                builderC.redirectErrorStream(true);
                final Process processC = builderC.start();
                // Watch the process
                watch(processC);
                processC.waitFor();

                //update FM_size [bytes] with size of fasta file
                file = new File(pathfile_ESS_in + FMind_c_EXT);
                if (file.exists() && file.isFile()) {
                    FM_size_dk += file.length();
                } else {
                    System.out.println(pathfile_ESS_in + FMind_c_EXT + " not a file or non-existent file !!!.");
                }
            } else {
                System.out.println(pathfile_ESS_in + " not a file or non-existent file !!!.");
            }
        }

        DP_SD_RAM_size.print(uncompr_ESS_size_dk + ", ");
        //complete 1_th row on DP_SD_RAM_size
        DP_SD_RAM_size.println(FM_size_dk + "\n");
    }
}
