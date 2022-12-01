import java.io.*;
import java.nio.file.FileSystems;

public class KmerDSKUtility_exp {

    private KmerDSKUtility_exp() {}

    protected static void execute (String path, String DSK_in, int mode_Sorting, String DSK_out_sort_fragments, int b, int debug) throws InterruptedException, IOException {

        /*
            0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/     (initial path to read and write files)
            1   ->  Staph-k4-DSK                (DSK input file in (<Dk> <Fk>) txt format)
            2   ->  1                           (0.verbatim - 1.lexicogr - 2.frequency)
            3   ->  Staph-k4-DSK-lexicogr       (DSK output file divided into fragments sorted in according to previous parameter)
            4   ->  5                           (number of rows per file =  Math.pow(2,num_bits) -1)
            5   ->  0/1                         0 no debug, 1 verbose mode

         */

        String DICT_EXT = ".txt";
        String Dk_EXT = ".fasta";
        String Fk_EXT = ".txt";
        String nfo_EXT = ".nfo";

        String VERBATIM_Suffix = "-verbatim";
        String LEXICOGR_Suffix = "-lexicogr";
        String FREQUENCY_Suffix = "-frequency";

        String FS_SEPARATOR = FileSystems.getDefault().getSeparator();

        String aux_mode = "";

        // each fragment, except possibly the last one will contain exactly ( 2^b -1 ) items
        int num_rows_file = (int) (Math.pow(2,b) - 1);

        // Create folder for storage fragments and update into mode_Sorting following step
        String folder = "mkdir " + path + DSK_out_sort_fragments + FS_SEPARATOR;

        if (debug==1) System.out.println(folder);

        Runtime run = Runtime.getRuntime();
        Process pr = run.exec(folder);
        pr.waitFor();

        // Perform input file sorting according to the type of sorting required: verbatim (nothing to do), lexicogr, frequency
        switch (mode_Sorting) {

            case 0: {
                System.out.println("Ordering in according to DSK output (verbatim)\n");  //non faccio nulla
                aux_mode = VERBATIM_Suffix;

                String cmd = "cp " + path + DSK_in + DICT_EXT + " " + path + DSK_in + aux_mode + DICT_EXT;

                if (debug==1) System.out.println(cmd);

                run = Runtime.getRuntime();
                pr = run.exec(cmd);
                pr.waitFor();
                break;

            }

            case 1: {
                System.out.println("Ordering DSK output via K-mer (lexicogr)\n");   // Ordino per Kmero
                aux_mode = LEXICOGR_Suffix;

                String cmd = "sort -k1,1 " + path + DSK_in + DICT_EXT + " -o " + path + DSK_in + aux_mode + DICT_EXT;

                if (debug==1) System.out.println(cmd);

                run = Runtime.getRuntime();
                pr = run.exec(cmd);
                pr.waitFor();
                break;
            }

            case 2: {
                System.out.println("Ordering DSK output via Frequency (frequency)\n");   // Ordino per Frequenza
                aux_mode = FREQUENCY_Suffix;

                String cmd = "sort -k2,2 -n " + path + DSK_in + DICT_EXT + " -o " + path + DSK_in + aux_mode + DICT_EXT;

                if (debug==1) System.out.println(cmd);

                run = Runtime.getRuntime();
                pr = run.exec(cmd);
                pr.waitFor();
                break;
            }
        }

        // Input D in the format <Kmer> <frequency> rearranged according to the mode_Sorting parameter
        BufferedReader dict_in = new BufferedReader(new FileReader(path + DSK_in + aux_mode + DICT_EXT));

        String line;
        String[] aux;

        int id_fragment = 0;
        int count_file = 0;
        long cumulative_rows = 0;
        long total_rows = 0;
        int num_fragments;

        String kmer;
        int  freq;

        // FIXME: Create kmer and frequency fragments file
        // Output Kmers file in FASTA format
        PrintWriter kmers_out = new PrintWriter(new FileWriter(path + DSK_out_sort_fragments + FS_SEPARATOR + DSK_out_sort_fragments + "_" + id_fragment + Dk_EXT));

        // Output Frequency file in <frequency> \n format
        PrintWriter frequency_out = new PrintWriter(new FileWriter(path + DSK_out_sort_fragments + FS_SEPARATOR + DSK_out_sort_fragments + "_" + id_fragment + Fk_EXT));

        // I write the character '>' in the first line of the output FASTA file
        kmers_out.println(">");

        //initialize k value to zero
        int k = 0;

        while ((line = dict_in.readLine()) != null) {

            if (count_file < num_rows_file) {

                aux = line.split("\\s+");

                kmer = aux[0];
                freq = Integer.parseInt(aux[1]);

                // In the first fragment, in the first line of the file I read the length of the k-mer, which is k.
                // In the other fragments it will not evaluate the expression because count_file will be initialized to 1.
                if (count_file == 0) {
                    k =  kmer.length();
                }

                // storage on output files
                kmers_out.println(kmer);
                frequency_out.println(freq);

                count_file++;
            }
            else {

                // num_rows_file has been saved !
                cumulative_rows += count_file;
                if (debug==1) System.out.println("Save first " + cumulative_rows + " rows from input file to Dk and Fk output files");

                count_file = 1; // storing first row on both files
                id_fragment++;

                /*
                //TODO: applica qui la correzione del file dei k-mers per un corretto funzionamento del compressore specialistico SPRING
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

                // Fragment closure
                kmers_out.close();
                frequency_out.close();


                // FIXME: Reopen kmer and frequency fragments file
                // Output Kmers file in FASTA format
                kmers_out = new PrintWriter(new FileWriter(path + DSK_out_sort_fragments + FS_SEPARATOR + DSK_out_sort_fragments + "_" + id_fragment + Dk_EXT));

                // Output Frequency file in <frequency> \n format
                frequency_out = new PrintWriter(new FileWriter(path + DSK_out_sort_fragments + FS_SEPARATOR + DSK_out_sort_fragments + "_" + id_fragment + Fk_EXT));

                // I write the character '>' in the first line of the output FASTA file
                kmers_out.println(">");

                // Storage of the current row on the output files
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

        /* Input Dictionary closure */
        dict_in.close();

        /*
        //TODO: applica qui la correzione del file dei k-mers per un corretto funzionamento del compressore specialistico SPRING
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

        // Output Files Closure
        frequency_out.close();
        kmers_out.close();

        // FIXME: Creazione, caricamento info e finalizzazione file contenente i dati di informazione (NFO) per la ricostruzione dei files interi
        PrintWriter info_out = new PrintWriter(new FileWriter(path + DSK_out_sort_fragments + FS_SEPARATOR + DSK_out_sort_fragments + nfo_EXT));


        //FIXME: storage info
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

        //FIXME: PLEASE NOTE: I do not remove the sorted file created as a result of desired sorting, I use it for the 'BijectionFkFromDSKTOSPSSKmer_exp algorithm.

        // Cleanup Garbage Collector
        System.gc();

    }

    public static void main(String[] args) throws IOException, InterruptedException {

        /*
        // Example of input:
            /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/     (String)
            Staph-k4-DSK                                                    (String)
            1                                                               (int)
            Staph-k4-DSK-lexicogr                                           (String)
            5                                                               (int)
            1  (0 no debug, 1 verbose mode)                                 (int)
         */

        KmerDSKUtility_exp.execute("/home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/", "Staph-k4-DSK", 1, "Staph-k4-DSK-lexicogr", 5, 1);
    }
}
