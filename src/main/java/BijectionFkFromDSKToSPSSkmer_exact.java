import java.io.*;
import java.nio.file.FileSystems;

public class BijectionFkFromDSKToSPSSkmer_exact {

    private BijectionFkFromDSKToSPSSkmer_exact() {}

    protected static void execute (String path, String ESS_in, int k, String DSK_in_sort_fragments, String Fk_out_Bijection, int num_bits, int debug) throws InterruptedException, IOException {

        /*
            (F'k)
            0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/     (initial path to read and write files)
            1   ->  Staph-ESS-k4                (ESS input file in fasta format)
            2   ->  4                           k value of k-mer tokenization of reads
            3   ->  Staph-DSK-k4-verbatim       (DSK output folder and file divided into fragments sorted in according to mode_Sorting parameter of KmerDSKUtility_exp)
            4   ->  Staph-DSK-k4-bijection      (F'k frequency files)
            5   ->  5                           (number of rows per file =  Math.pow(2,num_bits) -1)
            6   ->  0/1                         no/yes verbose mode
        */

        String SPSS_EXT = ".fasta";
        String DSK_EXT = ".txt";
        String FkOut_EXT = ".txt";
        String nfo_EXT = ".nfo";
        String FS_SEPARATOR = FileSystems.getDefault().getSeparator();

        // SPSS S reading in the FASTA format
        BufferedReader S = null;
        try {
            S = new BufferedReader(new FileReader(path + ESS_in + SPSS_EXT));
            if (debug==1) System.out.println("ESS: " + path + ESS_in + SPSS_EXT);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Read <Dk> <Fk> from DSK in txt format rearranged by KmerDSKUtility_exp
        BufferedReader DSK = null;
        try {
            DSK = new BufferedReader(new FileReader(path + DSK_in_sort_fragments + DSK_EXT));
            if (debug==1) System.out.println("DSK Dictionary (rearranged): " + path + DSK_in_sort_fragments + DSK_EXT);
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


        // FIXME: Scanning S for storage: <kmer> <IdSeq> <IdPos> on SPSS.txt temporary  file in txt format

        // auxiliary string variable read rows input files
        String l;

        // sequence variable index in SPSS
        long idSeq = 0;

        // count index variable k-mer encountered
        long countESS = 0;

        while ((l = S.readLine()) != null) {
            // I skip lines that begin with the character '>'
            if (!l.startsWith(">")) {
                // Progressive sequence number updating
                ++idSeq;

                // K-mer extraction from SPSS sequence
                for (int i = 0; i < l.length() - k + 1; i++) {
                    String kmer = l.substring(i, i + k);

                    // save data to temporary file SPSS.txt
                    countESS++;

                    //FIXME: For each x of SPSS S, I consider its canonical with respect to DSK order relation, i.e., A<C<T<G
                    //  see https://github.com/GATB/dsk#kmers-and-their-reverse-complements
                    SPSS_tmp.println(StringUtils.Canonical_DSK_RelOrd(kmer) + " " + idSeq + " " + i);
                }
            }
        }

        //DEBUG
        if (debug==1) System.out.println("kmers dictionary size (possibly with duplicates): " + countESS);

        // Closing the SPSS S and SPSS_tmp input files SPSS.txt
        S.close();
        SPSS_tmp.close();


        // FIXME: sorting SPSS.txt with k-mer and IdSeq keys
        // sorting SPSS.txt (KEY: canonical k-mer) and returning temporary file SPSS-sort.txt
        String cmd = "sort -k1,1 -k2,2n " + path + "SPSS" + DSK_EXT + " -o " + path + "SPSS" + "-sort" + DSK_EXT;
        if (debug==1) System.out.println(cmd);

        Runtime run = Runtime.getRuntime();
        Process pr = run.exec(cmd);
        pr.waitFor();

        // FIXME: sorting "DSK".txt (KEY: k-mer) and returning temporary file "DSK"-sort.txt
        cmd = "sort -k1,1 " + path + DSK_in_sort_fragments + DSK_EXT + " -o " + path + DSK_in_sort_fragments + "-sort" + DSK_EXT;
        if (debug==1) System.out.println(cmd);

        run = Runtime.getRuntime();
        pr = run.exec(cmd);
        pr.waitFor();


        //FIXME: Scan D'k with Dk in sorted versions to associate x's in D'k with the Fk frequency that matches Dk
        //      For each match, I keep output files with matches (IdSeq.IdPos_in_S, frequency) that I will sort by position
        //      and later I will throw out avenue positions by memorizing the frequencies. So in post-processing, from the k-mers scan of S,
        //      calculating for each the canonical, it will match with the frequency file obtained as output in the pre-processing step.


        //FIXME: I create a temporary file with <IdSeq.IdPos> from the SPSS_tmp file (sorted) and with <Fk> from the DSK file (sorted).
        // Can I anticipate on-line acquisition of the information? - No, because they are being sorted.

        //rEAD FROM temporary file SPSS-sort.txt
        BufferedReader SPSS_sort = null;
        try {
            SPSS_sort = new BufferedReader(new FileReader(path + "SPSS" + "-sort" + DSK_EXT));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        //Read from temporary file "DSK"-sort.txt
        BufferedReader DSK_sort = null;
        try {
            DSK_sort = new BufferedReader(new FileReader(path + DSK_in_sort_fragments + "-sort" + DSK_EXT));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        String line_SPSS_sort;
        String line_DSK_sort;
        String[] aux_SPSS;
        String[] aux_DSK;

        long count_SPSS = 0, count_DSK = 0, count_mismatch = 0;

        for (long i = 0; i < countESS ; i++) {

            line_SPSS_sort = SPSS_sort.readLine();
            line_DSK_sort = DSK_sort.readLine();

            if (line_SPSS_sort.compareTo("") != 0) count_SPSS++;
            if (line_DSK_sort.compareTo("") != 0) count_DSK++;

            aux_SPSS = line_SPSS_sort.split("\\s+");
            aux_DSK = line_DSK_sort.split("\\s+");

            if (aux_SPSS[0].compareTo(aux_DSK[0]) != 0) {

                // update mismatch counter
                count_mismatch++;

                if (debug==1) {
                    System.out.println("N. " + i + " mismatch at position: " + i);
                    System.out.println(i + ": " + "SPSS[0]: " + aux_SPSS[0] + "\t" + "SPSS[1]: " + aux_SPSS[1] + "\t" + "SPSS[2]: " + aux_SPSS[2] + "\t" + "DSK[0]: " + aux_DSK[0] + "\t" + "DSK[1]: " + aux_DSK[1] + "\t\t ==> \t\t" + "Bijection file: " + aux_SPSS[1] + " " + aux_SPSS[2] + " " + aux_DSK[1]);
                }

                // scan from keyboard
                BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
                System.out.println("Proceed Anyway ... ");
                try {
                    String s = br.readLine();
                } catch (Exception e) {
                    System.out.println(e);
                }
            }

            SPSS_DSK_Bijection.println(aux_SPSS[1] + " " + aux_SPSS[2] + " " + aux_DSK[1]);
        }

        // Close temporary bijection file: SPSS_DSK_Bijection.txt
        SPSS_DSK_Bijection.close();

        // Close temporary file "DSK"-sort.txt
        DSK_sort.close();

        // Close temporary file SPSS-sort.txt
        SPSS_sort.close();


        // Sort newly created file (SPSS_DSK_Bijection.txt) (KEY: IdSeq.IdPos)
        cmd = "sort -k1,1n -k2,2n " + path + "SPSS_DSK_Bijection" + DSK_EXT + " -o " + path + "SPSS_DSK_Bijection" + "-sort" + DSK_EXT;
        if (debug==1) System.out.println(cmd);

        run = Runtime.getRuntime();
        pr = run.exec(cmd);
        pr.waitFor();

        // FIXME: From the sorted file (SPSS_DSK_Bijection-sort.txt), I generate Fk_out
        // I read from the temporary file SPSS_DSK_Bijection-sort.txt
        BufferedReader SPSS_DSK_Bijection_sort = null;
        try {
            SPSS_DSK_Bijection_sort = new BufferedReader(new FileReader(path + "SPSS_DSK_Bijection" + "-sort" + DSK_EXT));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // I write to Fk_out in a dedicated folder (DP3), as in the DP0, DP1, DP2 cases.
        // Create folder for DP3 scenario
        String folder = "mkdir " + path + Fk_out_Bijection + FS_SEPARATOR;
        run = Runtime.getRuntime();
        pr = run.exec(folder);
        pr.waitFor();

        // FIXME: Compute num_rows_file information
        int num_rows_file = (int) (Math.pow(2,num_bits) -1);

        int id_fragment = 0;
        int count_file = 0;
        int cumulative_rows = 0;
        long total_rows = 0;
        int num_fragments = 0;

        String freq;

        String line_Bijection_sort = null;
        String[] aux_Bijection_sort;

        PrintWriter Fk_out = null;
        try {
            Fk_out = new PrintWriter(new FileWriter(path + Fk_out_Bijection + FS_SEPARATOR + Fk_out_Bijection + "_" + id_fragment + FkOut_EXT));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        while ( ((line_Bijection_sort = SPSS_DSK_Bijection_sort.readLine()) != null) ) {

            aux_Bijection_sort = line_Bijection_sort.split("\\s+");
            //storage su current fragment files di output
            //0: IdSeq; 1:IdPos; 2: frequenza
            freq = aux_Bijection_sort[2];

            if (count_file < num_rows_file) {

                Fk_out.println(freq);

                count_file++;
            }
            else {
                // num_rows_file has been saved !
                cumulative_rows += count_file;

                count_file = 1; // storing first row new fragment frequency output files
                id_fragment++;

                // Close current fragment
                Fk_out.close();

                // FIXME: Reopen kmer and frequency fragments file
                Fk_out = new PrintWriter(new FileWriter(path + Fk_out_Bijection + FS_SEPARATOR + Fk_out_Bijection + "_" + id_fragment + FkOut_EXT));

                Fk_out.println(freq);

            }
            total_rows = cumulative_rows + count_file;
        }
        num_fragments = id_fragment + 1;

        // Close tmp SPSS_DSK_Bijection-sort.txt
        SPSS_DSK_Bijection_sort.close();

        // Close output file
        Fk_out.close();


        // FIXME: Creating, loading info and finalizing file containing information data (NFO) for whole file reconstruction
        PrintWriter info_out = new PrintWriter(new FileWriter(path + Fk_out_Bijection + FS_SEPARATOR + Fk_out_Bijection + nfo_EXT));

        //storage info
        info_out.println(total_rows);

        info_out.println(k);

        if (total_rows < num_rows_file)
            info_out.println(total_rows);
        else
            info_out.println(num_rows_file);

        info_out.println(num_fragments);

        info_out.close();


        // TODO: Remove Temporary Files
        System.out.println("");

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
        rm_file = "rm -f " + path + DSK_in_sort_fragments + "-sort" + DSK_EXT;
        if (debug==1) System.out.println(rm_file);
        run = Runtime.getRuntime();
        pr = run.exec(rm_file);
        pr.waitFor();

        // Cleanup Garbage Collector
        System.gc();
    }


    public static void main(String[] args) throws IOException, InterruptedException {

        // Example of input
        /*
            0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/     (initial path to read and write files)
            1   ->  Staph-k4-ESS                (ESS input file in fasta format)
            2   ->  4                           k value of k-mer tokenization of reads
            3   ->  Staph-k4-DSK-verbatim       (DSK output folder and file divided into fragments sorted in according to mode_Sorting parameter of KmerDSKUtility_exp)
            4   ->  Staph-k4-bijection          (F'k frequency files)
            5   ->  5                           (number of rows per file =  Math.pow(2,num_bits) -1)
            6   ->  0/1                         no/yes verbose mode
        */
        BijectionFkFromDSKToSPSSkmer_exact.execute("/home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/", "Staph-k4-ESS", 4, "Staph-k4-DSK-verbatim", "Staph-k4-bijection", 5, 1);
    }
}
