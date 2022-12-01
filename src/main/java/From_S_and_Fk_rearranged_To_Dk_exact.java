import java.io.*;
import java.nio.file.FileSystems;

public class From_S_and_Fk_rearranged_To_Dk_exact {

    private From_S_and_Fk_rearranged_To_Dk_exact() {
    }

    protected static void execute(String path, String ESS_in, int k, String Dk_Bijection, int num_bits, int debug) throws InterruptedException, IOException {

         /*
            (D'k)
            0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/     (initial path to read and write files)
            1   ->  Staph-k4-ESS                (ESS input file in fasta format)
            2   ->  4                           k value of k-mer tokenization of reads
            3   ->  Staph-k4-Bijection          (Dk output folder and file divided into fragments)
            4   ->  5                           (number of rows per file =  Math.pow(2,num_bits) -1)
            5   ->  0/1                         no/yes verbose mode
         */


        // TODO: PLEASE NOTE - I don't actually make use of Fk rearranged, because, by construction I sorted it by IdSeq.IdPos from 1.0 through xx.yy .
        //  It guarantees me a match with k-mers in S transformed to canonical according to the order relation used by DSK.

        String SPSS_EXT = ".fasta";
        String DSK_EXT = ".fasta";
        String nfo_EXT = ".nfo";
        String TMP_EXT = ".tmp";
        String FS_SEPARATOR = FileSystems.getDefault().getSeparator();

        int num_rows_file = 0;
        int num_blocks = 0;
        String l;

        // FIXME: Read num_Blocks input parameter into .nfo file
        BufferedReader S_nfo = null;
        try {
            S_nfo = new BufferedReader(new FileReader(path + ESS_in + nfo_EXT));
            if (debug==1) System.out.println("ESS info file: " + path + ESS_in + nfo_EXT);
            l = S_nfo.readLine();
            num_blocks = Integer.parseInt(l);
            S_nfo.close();
            if (debug==1) System.out.println("num_blocks = " + num_blocks);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // FIXME: Compute num_rows_file information
        num_rows_file = (int) (Math.pow(2, num_bits) - 1);

        // FIXME: Storage Dk k-mer temporary dictionary
        // Read SPSS S in FASTA format (SPSS_EXT))
        BufferedReader S = null;
        try {
            // Output Dk_out temporary file in .fasta format
            PrintWriter Dk_out = null;
            try {
                Dk_out = new PrintWriter(new FileWriter(path + ESS_in + TMP_EXT));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            S = new BufferedReader(new FileReader(path + ESS_in + SPSS_EXT));

            // scan S to obtain from the k-mer(S) the respective canonicals  ones (with the DSK order relation)
            while ((l = S.readLine()) != null) {

                // I skip lines that begin with the character '>'
                if (!l.startsWith(">")) {

                    // Extraction of k-mer from SPSS current sequence
                    for (int i = 0; i < l.length() - k + 1; i++) {
                        String kmer = l.substring(i, i + k);

                        //TODO: For each x of SPSS S, I consider its canonical with respect to DSK relation order, i.e., A<C<T<G .
                        //  see   https://github.com/GATB/dsk#kmers-and-their-reverse-complements

                        Dk_out.println(StringUtils.Canonical_DSK_RelOrd(kmer));
                    }
                }
            }

            // close Dk_out temporary file
            Dk_out.close();

            // close ESS input file
            S.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }


        // TODO: starting from (path + ESS_in + TMP_EXT) and Divide k-mer dictionary obtained from S into fragments
        int id_fragment = 0;
        int count_file = 0;
        int cumulative_rows = 0;
        long total_rows = 0;
        int num_fragments = 0;

        // Input Dk file in .tmp format (TMP_EXT)
        BufferedReader Dk_in = null;
        try {
            Dk_in = new BufferedReader(new FileReader(path + ESS_in + TMP_EXT));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Output Dk_out file in .fasta format
        PrintWriter Dk_out = null;
        try {
            Dk_out = new PrintWriter(new FileWriter(path + Dk_Bijection + FS_SEPARATOR + Dk_Bijection + "_" + id_fragment + DSK_EXT));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }


        // Save header D'k
        Dk_out.println(">");

        // scan Dk_in to get the various k-mer(S)
        while ((l = Dk_in.readLine()) != null) {

            // Saving canonical k-mer to current fragment file output Dk_out
            if (count_file < num_rows_file) {
                Dk_out.println(l);
                count_file++;
            } else {
                // num_rows_file has been saved !
                cumulative_rows += count_file;
                if (debug==1) System.out.println("Save first " + cumulative_rows + " rows to fragment output file Dk");

                count_file = 1; // storing first row on new fragment file
                id_fragment++;

                /*
                //TODO: applica qui la correzione del file dei k-mers per un corretto funzionamento del compressore specialistico SPRING
                if (k >= 7) {
                    if ((count_file + 2) % 4 != 0) {

                        int i_aux = 1;
                        Dk_out.println(">");

                        while ((count_file + 2 + i_aux) % 4 != 0) {
                            Dk_out.println("");

                            //update i_aux
                            i_aux++;
                        }
                    }
                }
                */

                // Chiusura frammento corrente
                Dk_out.close();

                // FIXME: Reopen kmer and frequency fragments file
                Dk_out = new PrintWriter(new FileWriter(path + Dk_Bijection + FS_SEPARATOR + Dk_Bijection + "_" + id_fragment + DSK_EXT));

                //Scrivo il carattere > nella prima riga del file fasta di output
                Dk_out.println(">");

                // save last k-mer readed into new fragment
                Dk_out.println(l);
            }
        }

        // DEBUG
        total_rows = cumulative_rows + count_file;
        if (debug==1) System.out.println("total_rows = " + total_rows);


        /*
        //TODO: applica qui la correzione del file dei k-mers per un corretto funzionamento del compressore specialistico SPRING
        if (k > 7) {
            if ((count_file + 2) % 4 != 0) {

                int i_aux = 1;
                Dk_out.println(">");

                while ((count_file + 2 + i_aux) % 4 != 0) {
                    Dk_out.println("");

                    //update i_aux
                    i_aux++;
                }
            }
        }
        */

        // Chiusura dei file di input Dk_in e dell'ultimo fragment file di output Dk_out
        Dk_in.close();
        Dk_out.close();

        //FIXME: remove temporary (path + ESS_in + TMP_EXT) file
        String cmd = "rm -f " + path + ESS_in + TMP_EXT;
        if (debug==1) System.out.println(cmd);

        Runtime run = Runtime.getRuntime();
        Process pr = run.exec(cmd);
        pr.waitFor();
    }


    public static void main(String[] args) throws IOException, InterruptedException {

        /*
           0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/     (initial path to read and write files)
           1   ->  Staph-k4-ESS                (ESS input file in fasta format)
           2   ->  4                           k value of k-mer tokenization of reads
           3   ->  Staph-k4-Bijection          (Dk output folder and file divided into fragments)
           4   ->  5                           (number of rows per file =  Math.pow(2,num_bits) -1)
           5   ->  0/1                         no/yes verbose mode
        */

        From_S_and_Fk_rearranged_To_Dk_exact.execute("/home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/", "Staph-k4-ESS", 4, "Staph-k4-Bijection", 5, 1);
    }
}
