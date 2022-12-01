import java.io.*;
import java.nio.file.FileSystems;

public class BijectionFkFromDSKToSPSSkmer_approx {

    // file System File Separator
    private static final String FS_SEPARATOR = FileSystems.getDefault().getSeparator();

    private BijectionFkFromDSKToSPSSkmer_approx() {}

    protected static void execute (String path, String exp, String ESS_in, int k, String DSK_in_sort_fragments, String Fk_out_Bijection, int num_blocks, int debug) throws InterruptedException, IOException {

        /*
            (F'k)
            0   ->  /home/user/IdeaProjects/TradeOff/experiments/Staph/k4/   (initial path to read and write files)
            1   ->  Staph                                                               (experiment dataset name)
            2   ->  Staph-ESS-k4                                                        (ESS input file in FASTA format)
            3   ->  4                                                                   (k value of k-mer tokenization of reads)
            4   ->  Staph_DSK-k4-verbatim                                               (DSK output folder and file divided into fragments sorted in according to mode_Sorting parameter of KmerDSKUtility_exp)
            5   ->  Staph-DSK-k4-bijection                                              (F'k frequency files)
            6   ->  num_blocks                                                          es. 10
            7   ->  0/1                                                                 no/yes verbose mode
        */

        String SPSS_EXT = ".fasta";
        String DSK_EXT = ".txt";
        String FkOut_EXT = ".f00";
        String nfo_EXT = ".nfo";

        for (int i = 0; i < num_blocks; i++) {      // iterate num_blocks iterations

            String block_ess = path + ESS_in + "_Block_" + i + SPSS_EXT;
            String block_dsk = path + exp + "-DSK-k" + k + "-verbatim_" + i + DSK_EXT;

            // Reading SPSS Block S_i in FASTA format.
            BufferedReader S_i = null;
            try {
                S_i = new BufferedReader(new FileReader(block_ess));
                if (debug==1) System.out.println("ESS: " + block_ess);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            // Reading Block <Dk> <Fk> from in txt format.
            BufferedReader DSK_i = null;
            try {
                DSK_i = new BufferedReader(new FileReader(block_dsk));
                if (debug==1) System.out.println("DSK: " + block_dsk);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            // Writing temporary file SPSS.txt (useful for bijection)
            PrintWriter SPSS_tmp = null;
            try {
                SPSS_tmp = new PrintWriter(new FileWriter(path + "SPSS" + DSK_EXT));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }


            // Writing temporary file SPSS_DSK_Bijection.txt (bijection)
            PrintWriter SPSS_DSK_Bijection = null;
            try {
                SPSS_DSK_Bijection = new PrintWriter(new FileWriter(path + "SPSS_DSK_Bijection" + DSK_EXT));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }


            // Scan S_i by storage: <kmer> <IdSeq> <IdPos> on SPSS.txt tmp file in txt format
            // auxiliary string variable read rows input files
            String l;

            // sequence index variable in SPSS
            long idSeq = 0;

            // count index variable k-mer encountered
            long countESS = 0;

            while ((l = S_i.readLine()) != null) {
                // skip lines that begin with the character '>'
                if (!l.startsWith(">")) {

                    // Updating Progressive sequence number
                    ++idSeq;

                    // K-mer extraction from SPSS sequence
                    for (int j = 0; j < l.length() - k + 1; j++) {
                        String kmer = l.substring(j, j + k);

                        // saving data to temporary file SPSS.txt
                        countESS++;

                        //FIXME: For each x of SPSS S, consider its canonical with respect to  DSK relation order, i.e., A<C<T<G .
                        //  see https://github.com/GATB/dsk#kmers-and-their-reverse-complements
                        SPSS_tmp.println(StringUtils.Canonical_DSK_RelOrd(kmer) + " " + idSeq + " " + j);
                    }
                }
            }

            // Closing the SPSS S and SPSS_tmp input files
            S_i.close();
            SPSS_tmp.close();


            // FIXME: sorting SPSS.txt with KEYS k-mer and IdSeq
            // sorting SPSS.txt (KEY: canonical k-mer) and returning temporary file SPSS-sort.txt
            String cmd = "sort -k1,1 -k2,2n " + path + "SPSS" + DSK_EXT + " -o " + path + "SPSS" + "-sort" + DSK_EXT;
            if (debug==1) System.out.println(cmd);

            Runtime run = Runtime.getRuntime();
            Process pr = run.exec(cmd);
            pr.waitFor();


            // FIXME: sorting DSK with KEY k-mer
            //  sorting "DSK".txt (KEY: k-mer) and returning temporary file "DSK"-sort.txt
            cmd = "sort -k1,1 " + path + DSK_in_sort_fragments + "_" + i + DSK_EXT + " -o " + path + DSK_in_sort_fragments + "-sort" + "_" + i + DSK_EXT;
            if (debug==1) System.out.println(cmd);

            run = Runtime.getRuntime();
            pr = run.exec(cmd);
            pr.waitFor();


            //FIXME: Scan D'k with Dk in sorted versions to associate x's in D'k with the Fk frequency that matches Dk
            // For each match, keep files with the matches (IdSeq.IdPos_in_S, frequency) in output, which I will sort by position
            // and later throw away the positions by storing the frequencies. So in the post-processing phase, from scanning the k-mers of S,
            // calculating for each the canon, it will match with the frequency file obtained in output in the pre-processing phase.


            //FIXME: I create a temporary file with <IdSeq.IdPos> from the SPSS_tmp file (sorted) and with <Fk> from the DSK file (sorted).
            // Can I anticipate on-line acquisition of the information? - No, because they are being sorted.

            // Read from temporary file SPSS-sort.txt
            BufferedReader SPSS_sort = null;
            try {
                SPSS_sort = new BufferedReader(new FileReader(path + "SPSS" + "-sort" + DSK_EXT));
                if (debug==1) System.out.println("SPSS-sort: " + path + "SPSS" + "-sort" + DSK_EXT);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            // Read from temporary file  "DSK"-sort.txt
            BufferedReader DSK_sort = null;
            try {
                System.out.println("DSK-sort = " + path + DSK_in_sort_fragments + "-sort" + "_" + i + DSK_EXT);
                DSK_sort = new BufferedReader(new FileReader(path + DSK_in_sort_fragments + "-sort" + "_" + i + DSK_EXT));
                if (debug==1) System.out.println("DSK-sort: " + path + DSK_in_sort_fragments + "-sort" + "_" + i + DSK_EXT);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            //start bijection
            String line_SPSS_sort;
            String line_DSK_sort;
            String[] aux_SPSS;
            String[] aux_DSK;

            long count_SPSS = 0, count_DSK = 0, count_mismatch = 0;

            for (long index = 0; index < countESS ; index++) {

                line_SPSS_sort = SPSS_sort.readLine();
                line_DSK_sort = DSK_sort.readLine();

                if (line_SPSS_sort.compareTo("") != 0) count_SPSS++;
                if (line_DSK_sort.compareTo("") != 0) count_DSK++;

                aux_SPSS = line_SPSS_sort.split("\\s+");
                aux_DSK = line_DSK_sort.split("\\s+");

                if (aux_SPSS[0].compareTo(aux_DSK[0]) != 0) {

                    count_mismatch++;
                    System.out.println("N. " + i + " mismatch at position: " + i);

                    //DEBUG
                    if (debug==1) System.out.println(i + ": " + "SPSS[0]: " + aux_SPSS[0] + "\t" + "SPSS[1]: " + aux_SPSS[1] + "\t" + "SPSS[2]: " + aux_SPSS[2] + "\t" + "DSK[0]: " + aux_DSK[0] + "\t" + "DSK[1]: " + aux_DSK[1] + "\t\t ==> \t\t" + "Bijection file: " + aux_SPSS[1] + " " + aux_SPSS[2] + " " + aux_DSK[1]);

                    // scan from keyboard
                    BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
                    System.out.println("Proceed Anyway ... Press any key !!!");
                    try {
                        String s = br.readLine();
                    } catch (Exception e) {
                        System.out.println(e);
                    }
                }
                SPSS_DSK_Bijection.println(aux_SPSS[1] + " " + aux_SPSS[2] + " " + aux_DSK[1]);
            }

            //DEBUG
            if (debug==1) {
                System.out.println("Conteggio righe SPSS sorted: " + count_SPSS);
                System.out.println("Conteggio righe DSK sorted: " + count_DSK);
                System.out.println("Conteggio mismatch: " + count_mismatch);
            }

            // Closing temporary bjection file: SPSS_DSK_Bijection.txt
            SPSS_DSK_Bijection.close();

            // Closing temporary file "DSK"-sort.txt
            DSK_sort.close();
            // Closing temporary file SPSS-sort.txt
            SPSS_sort.close();


            // Sort newly created file (SPSS_DSK_Bijection.txt) (KEY: IdSeq.IdPos)
            cmd = "sort -k1,1n -k2,2n " + path + "SPSS_DSK_Bijection" + DSK_EXT + " -o " + path + "SPSS_DSK_Bijection" + "-sort" + DSK_EXT;
            if (debug==1) System.out.println(cmd);

            run = Runtime.getRuntime();
            pr = run.exec(cmd);
            pr.waitFor();

            // From the sorted file (SPSS_DSK_Bijection-sort.txt),  generate Fk_out. From the sorted file (SPSS_DSK_Bijection-sort.txt), generate Fk_out .
            // Read from temporary file SPSS_DSK_Bijection-sort.txt
            BufferedReader SPSS_DSK_Bijection_sort = null;
            try {
                SPSS_DSK_Bijection_sort = new BufferedReader(new FileReader(path + "SPSS_DSK_Bijection" + "-sort" + DSK_EXT));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            // write to Fk_out in a dedicated folder (DP3), as in the DP0, DP1, DP2 Cases.
            // Create folder for DP3 scenario
            String folder = "mkdir " + path + Fk_out_Bijection + FS_SEPARATOR;
            run = Runtime.getRuntime();
            pr = run.exec(folder);
            pr.waitFor();


            //FIXME: Create F'k i_th Block
            PrintWriter Fk_i_out = null;
            try {
                Fk_i_out = new PrintWriter(new FileWriter(path + Fk_out_Bijection + FS_SEPARATOR + Fk_out_Bijection + "_" + i + FkOut_EXT));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            String freq;
            String line_Bijection_sort = null;
            String[] aux_Bijection_sort;
            while ( ((line_Bijection_sort = SPSS_DSK_Bijection_sort.readLine()) != null) ) {

                aux_Bijection_sort = line_Bijection_sort.split("\\s+");
                //storage on current fragment output files
                //0: IdSeq; 1:IdPos; 2: frequency
                Fk_i_out.println(aux_Bijection_sort[2]);

            }
            // Closing  file tmp SPSS_DSK_Bijection-sort.txt
            SPSS_DSK_Bijection_sort.close();
            // Closing output files Fk_i
            Fk_i_out.close();



            // FIXME: Remove Temporary Files
            // remove SPSS.txt
            String rm_file = "rm -f " + path + "SPSS" + DSK_EXT;
            if (debug==1) System.out.println(rm_file);
            run = Runtime.getRuntime();
            pr = run.exec(rm_file);
            pr.waitFor();

            // remove SPSS-sort.txt
            rm_file = "rm -f " + path + "SPSS" + "-sort" + DSK_EXT;
            if (debug==1) System.out.println(rm_file);
            run = Runtime.getRuntime();
            pr = run.exec(rm_file);
            pr.waitFor();

            // remove SPSS_DSK_Bijection.txt
            rm_file = "rm -f " + path + "SPSS_DSK_Bijection" + DSK_EXT;
            if (debug==1) System.out.println(rm_file);
            run = Runtime.getRuntime();
            pr = run.exec(rm_file);
            pr.waitFor();

            // remove SPSS_DSK_Bijection-sort.txt
            rm_file = "rm -f " + path + "SPSS_DSK_Bijection" + "-sort" + DSK_EXT;
            if (debug==1) System.out.println(rm_file);
            run = Runtime.getRuntime();
            pr = run.exec(rm_file);
            pr.waitFor();

            // remove "DSK_in_sort_fragments"-sort.txt
            rm_file = "rm -f " + path + DSK_in_sort_fragments + "-sort" + "_" + i + DSK_EXT;
            if (debug==1) System.out.println(rm_file);
            run = Runtime.getRuntime();
            pr = run.exec(rm_file);
            pr.waitFor();

            // Clean Garbage Collector
            System.gc();
        }
    }

    public static void main(String[] args) throws IOException, InterruptedException {

        // Example of input
        /*
            0   ->  /home/user/IdeaProjects/TradeOff/exp-fragments/     (initial path to read and write files)
            1   ->  Staph                       (experiment Dataset name)
            2   ->  Staph-k4-ESS                (ESS input file in FASTA format)
            3   ->  4                           k value of k-mer tokenization of reads
            4   ->  Staph-k4-DSK-verbatim       (DSK output folder and file divided into fragments sorted in according to mode_Sorting parameter of KmerDSKUtility_exp)
            5   ->  Staph-k4-bijection          (F'k frequency files)
            6   ->  num_blocks                  es. 10
            7   ->  0/1                         no/yes verbose mode
        */
        BijectionFkFromDSKToSPSSkmer_approx.execute("/home/user/IdeaProjects/TradeOff/exp-fragments/", "Staph", "Staph-k4-ESS", 4, "Staph-k4-DSK-verbatim", "Staph-k4-bijection", 10, 1);
    }
}
