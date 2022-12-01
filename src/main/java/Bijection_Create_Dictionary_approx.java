import java.io.*;
import java.nio.file.FileSystems;

public class Bijection_Create_Dictionary_approx {

    private Bijection_Create_Dictionary_approx() {
    }

    private static final String FASTA_EXT = ".fasta";
    private static final String DSK_EXT = ".txt";
    private static final String DSK_TMP = ".tmp";
    private static final String D_TMP_EXT = ".d00";
    private static final String F_TMP_EXT = ".f00";
    private static final String nfo_EXT = ".nfo";

    // file System File Separator
    private static final String FS_SEPARATOR = FileSystems.getDefault().getSeparator();


    protected static void execute(String path, String exp, String ESS_in, int k, String Dk_Bijection, int num_bits, int num_blocks, int debug) throws InterruptedException, IOException {

         /*
            (D'k)
            0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-automazione/     (initial path to read and write files)
            1   ->  Staph                       (Experiment Dataset Name)
            2   ->  Staph-ESS-k4                (ESS input file in fasta format)
            3   ->  4                           k value of k-mer tokenization of reads
            4   ->  Staph-k4-Bijection          (Dk output folder and file divided into fragments)
            5   ->  5                           (number of rows per file =  Math.pow(2,num_bits) -1)
            6   ->  10                          num_blocks
            7   ->  0/1                         no/yes verbose mode

         */

        String DSK_approx_tmp = path + Dk_Bijection + FS_SEPARATOR + Dk_Bijection + "_approx" + DSK_TMP;
        PrintWriter DSK_approx_out = new PrintWriter(new FileWriter(DSK_approx_tmp));

        // Output DSK approximation output file in <k-mer> <frequency> \n format
        String DSK_approx = path + Dk_Bijection + FS_SEPARATOR + Dk_Bijection + "_approx" + DSK_EXT;

        for (int i = 0; i < num_blocks; i++) {      // iterate num_blocks iterations

            String block_Di = path + Dk_Bijection + FS_SEPARATOR + Dk_Bijection + "_" + i + D_TMP_EXT;
            String block_Fi = path + Dk_Bijection + FS_SEPARATOR + Dk_Bijection + "_" + i + F_TMP_EXT;

            // Input Di
            BufferedReader Di_in = new BufferedReader(new FileReader(block_Di));

            // Input Fi
            BufferedReader Fi_in = new BufferedReader(new FileReader(block_Fi));

            String lineD;
            String lineF;
            while ((lineD = Di_in.readLine()) != null) {
                lineF = Fi_in.readLine();

                DSK_approx_out.println(lineD + " " + lineF);
            }
            //close Di and Fi temporary files
            Di_in.close();
            Fi_in.close();

            //Cleanup Garbage Collector
            System.gc();
        }

        //close DSK_approx_out file
        DSK_approx_out.close();

        //FIXME: Sorting by k-mer and removing original .tmp
        String cmd = "sort -k1,1 " + DSK_approx_tmp + " -o " + DSK_approx_tmp + ".sort";
        if (debug==1) System.out.println(cmd);
        Runtime run = Runtime.getRuntime();
        Process pr = run.exec(cmd);
        pr.waitFor();


        //FIXME: remove DSK_approx temporary file
        cmd = "rm -f " + DSK_approx_tmp;
        if (debug==1) System.out.println(cmd);
        run = Runtime.getRuntime();
        pr = run.exec(cmd);
        pr.waitFor();

        //TODO: sum frequencies on DSK_approx.sort file and return dictionary
        CumulateFrequencies_from_DSKapprox(DSK_approx_tmp + ".sort", DSK_approx, debug);

        //FIXME: remove DSK_approx temporary file (sorted version)
        cmd = "rm -f " + DSK_approx_tmp + ".sort";
        if (debug==1) System.out.println(cmd);
        run = Runtime.getRuntime();
        pr = run.exec(cmd);
        pr.waitFor();


        //TODO: Divisione dizionario in frammenti secondo num_bits item per frammento
        Create_DSK_Fragments(DSK_approx, num_bits, path + Dk_Bijection + FS_SEPARATOR + Dk_Bijection + "_approx", debug);


        //TODO: remove temporary files: { *.d00 ; *.f00 }
        System.out.println("Deletion of *.d00 and *.f00 blocks ...");
        for (int i = 0; i < num_blocks; i++) {

            String block_Di = path + Dk_Bijection + FS_SEPARATOR + Dk_Bijection + "_" + i + D_TMP_EXT;
            String block_Fi = path + Dk_Bijection + FS_SEPARATOR + Dk_Bijection + "_" + i + F_TMP_EXT;

            //FIXME: remove Di temporary file (.d00)
            cmd = "rm -f " + block_Di;
            if (debug==1) System.out.println(cmd);
            run = Runtime.getRuntime();
            pr = run.exec(cmd);
            pr.waitFor();

            //FIXME: remove Fi temporary file (.f00)
            cmd = "rm -f " + block_Fi;
            if (debug==1) System.out.println(cmd);
            run = Runtime.getRuntime();
            pr = run.exec(cmd);
            pr.waitFor();
        }

        //TODO: remove DSK approx Dictionary: PLESE NOTE: It is best to delete it only after finding the Residue Dictionary
        cmd = "rm -f " + DSK_approx;
        if (debug==1) System.out.println(cmd);
        run = Runtime.getRuntime();
        pr = run.exec(cmd);
        pr.waitFor();
    }


    private static void Create_DSK_Fragments(String dsk_approx, int num_bits, String dsk_approx_fragment_base, int debug) throws IOException {

        // FIXME: Compute num_rows_file information
        int num_rows_file = (int) (Math.pow(2, num_bits) - 1);

        // Input DSK Approximation in <k-mer> <frequency> format
        BufferedReader dict_in = new BufferedReader(new FileReader(dsk_approx));

        String line;
        String[] aux;

        int id_fragment = 0;
        int count_file = 0;
        long cumulative_rows = 0;
        long total_rows = 0;
        int num_fragments;

        String kmer;
        int freq;

        // FIXME: Create kmer and frequency fragments file
        // Output Kmers file in FASTA format
        PrintWriter kmers_out = new PrintWriter(new FileWriter(dsk_approx_fragment_base + "_" + id_fragment + FASTA_EXT));

        // Output Frequency file in <frequency> \n format
        PrintWriter frequency_out = new PrintWriter(new FileWriter(dsk_approx_fragment_base + "_" + id_fragment + DSK_EXT));

        // write the character '>' in the first line of the output FASTA file
        kmers_out.println(">");

        //initialize k value to zero
        int k = 0;

        while ((line = dict_in.readLine()) != null) {

            if (count_file < num_rows_file) {

                aux = line.split("\\s+");

                kmer = aux[0];
                freq = Integer.parseInt(aux[1]);

                // In the first fragment, in the first line of the file, read the length of the k-mer, which is k.
                // In the other fragments it will not evaluate the expression because count_file will be initialized to 1.
                if (count_file == 0) {
                    k = kmer.length();
                }

                //storage on the output files
                kmers_out.println(kmer);
                frequency_out.println(freq);

                count_file++;
            } else {

                // num_rows_file has been saved !
                cumulative_rows += count_file;
                //System.out.println("Save first " + cumulative_rows + " rows from input file to Dk and Fk output files");

                count_file = 1; // storing first row on both files
                id_fragment++;

                /*
                //TODO: apply here k-mers file correction for proper operation of SPRING specialist compressor
                if (k > 7) {
                    if ((num_rows_file + 2) % 4 != 0) {

                        int i_aux = 1;
                        kmers_out.println(">");

                        while ((num_rows_file + 2 + i_aux) % 4 != 0) {
                            kmers_out.println("");

                            //update i_aux
                            i_aux++;
                        }
                    }
                }
                */

                // Closing Fragments
                kmers_out.close();

                frequency_out.close();

                //Cleanup Garbage Collector
                System.gc();

                // FIXME: Reopen kmer and frequency fragments file

                // Output Kmers file in FASTA format
                kmers_out = new PrintWriter(new FileWriter(dsk_approx_fragment_base + "_" + id_fragment + FASTA_EXT));

                // Output Frequency file in <frequency> \n format
                frequency_out = new PrintWriter(new FileWriter(dsk_approx_fragment_base + "_" + id_fragment + DSK_EXT));

                // write the character > in the first line of the output FASTA file
                kmers_out.println(">");

                // Storage on the current row files
                aux = line.split("\\s+");

                kmer = aux[0];
                freq = Integer.parseInt(aux[1]);

                kmers_out.println(kmer);
                frequency_out.println(freq);
            }
            total_rows = cumulative_rows + count_file;
        }
        num_fragments = id_fragment + 1;

        if (debug==1) {
            System.out.println("total_rows = " + total_rows);
            System.out.println("Created output fragment files !!!\n");
        }

        /* Closing files */
        dict_in.close();

        /*
        //TODO: apply here k-mers file correction for proper operation of SPRING specialist compressor
        if (k > 7) {
            if ((count_file + 2) % 4 != 0) {

                int i_aux = 1;
                kmers_out.println(">");

                while ((count_file + 2 + i_aux) % 4 != 0) {
                    kmers_out.println("");

                    //update i_aux
                    i_aux++;
                }
            }
        }
        */

        frequency_out.close();
        kmers_out.close();

        // FIXME: Creating, loading info and finalizing file containing information data (NFO) for whole file reconstruction.
        PrintWriter info_out = new PrintWriter(new FileWriter(dsk_approx_fragment_base + nfo_EXT));

        //storage info

        //total_rows (long)
        info_out.println(total_rows);

        // k (int)
        info_out.println(k);

        // rows for fragment (int)
        if (total_rows <= num_rows_file)
            info_out.println(total_rows);
        else
            // total_rows > num_rows_file
            info_out.println(num_rows_file);

        //num_fragments (int)
        info_out.println(num_fragments);

        info_out.close();
    }

    private static void CumulateFrequencies_from_DSKapprox(String tmp_sort, String dsk_approx, int debug) throws IOException {
        // Read from file DSK approximation Dictionary sort (with duplicates)
        BufferedReader DSK_tmp_sort = null;
        try {
            DSK_tmp_sort = new BufferedReader(new FileReader(tmp_sort));
            if (debug==1) System.out.println("Input - DSK approximation Dictionary sort (with duplicates): " + tmp_sort);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Write to file DSK approximation (without duplicates)
        PrintWriter DSK_out = null;
        try {
            DSK_out = new PrintWriter(new FileWriter(dsk_approx));
            if (debug==1) System.out.println("Output - DSK approximation Dictionary (without duplicates): " + dsk_approx);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        long count_distinct_kmer = 1;   // at the first comparison I have 2 different ones found

        String l = DSK_tmp_sort.readLine();
        String[] tokens = l.split(" ");
        String aux_kmer = tokens[0];
        int aux_freq = Integer.parseInt(tokens[1]);
        long count_rows = 1;    // first row encountered

        while ((l = DSK_tmp_sort.readLine()) != null) {

            //update count_rows
            count_rows++;

            tokens = l.split(" ");
            if (tokens[0].compareTo(aux_kmer) != 0) {   // Row with different k-mer than the previous one ...

                // distinct k-mer founded
                count_distinct_kmer++;

                // Storage (aux_kmer,aux_freq) and (tokens[0], token[1]) to DSK_out
                DSK_out.println(aux_kmer + " " + aux_freq);

                // update aux_kmer
                aux_kmer = tokens[0];
                aux_freq = Integer.parseInt(tokens[1]);

            } else { // cumulate the frequency with the current frequency partial

                // it is still the same kmer, so I update aux_kmer and cumulate the frequency
                aux_kmer = tokens[0];
                aux_freq += Integer.parseInt(tokens[1]);
            }
        }
        DSK_out.println(aux_kmer + " " + aux_freq);

        if (debug==1) {
            System.out.println("Rows scanned: " + count_rows);
            System.out.println("Distinct k-mers founded = " + count_distinct_kmer);
        }

        // close input and output files
        DSK_out.close();
        DSK_tmp_sort.close();
    }


    public static void main(String[] args) throws IOException, InterruptedException {

        /*
           0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-automazione/     (initial path to read and write files)
           1   ->  Staph                       (Experiment Dataset Name)
           2   ->  Staph-k4-ESS                (ESS input file in fasta format)
           3   ->  4                           k value of k-mer tokenization of reads
           4   ->  Staph-k4-Bijection          (Dk output folder and file divided into fragments)
           5   ->  5                           (number of rows per file =  Math.pow(2,num_bits) -1)
           6   ->  10                           num_blocks
           7   ->  0/1                         no/yes verbose mode
        */

        Bijection_Create_Dictionary_approx.execute("/home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-automazione/", "Staph", "Staph-k4-ESS", 4, "Staph-k4-Bijection", 5, 10, 1);
    }
}
