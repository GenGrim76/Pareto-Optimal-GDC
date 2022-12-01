import java.io.*;
import java.nio.file.FileSystems;

public class From_S_and_Fk_rearranged_To_Dk_approx {

    // file System File Separator
    private static final String FS_SEPARATOR = FileSystems.getDefault().getSeparator();

    private From_S_and_Fk_rearranged_To_Dk_approx() {
    }

    protected static void execute(String path, String exp, String ESS_in, int k, String Dk_Bijection, int num_blocks, int debug) throws InterruptedException, IOException {

         /*
            (D'k)
            0   ->  /home/user/IdeaProjects/TradeOff/experiments/     (initial path to read and write files)
            1   ->  Staph                       (Experiment Dataset Name)
            2   ->  Staph-k4-ESS                (ESS input file in FASTA format)
            3   ->  4                           k value of k-mer tokenization of reads
            4   ->  Staph-k4-Bijection          (Dk output folder and file divided into fragments)
            5   ->  num_blocks                  es. 10
            6   ->  0/1                         no/yes verbose mode
         */


        // TODO: PLEASE NOTE - I don't actually make use of Fk rearranged, because,
        //  by construction I sorted it by IdSeq.IdPos from 1.0 through xx.yy.
        //  It guarante a match with k-mers in S transformed to canonical, according
        //  to the order relation used by DSK.

        String SPSS_EXT = ".fasta";
        String DSK_EXT = ".d00";

        for (int i = 0; i < num_blocks; i++) {      // iterate num_blocks iterations

            String block_ess = path + ESS_in + "_Block_" + i + SPSS_EXT;

            String l;
            // FIXME: Storage Dk k-mer temporary dictionary
            // Read SPSS S_i in FASTA format (SPSS_EXT))
            BufferedReader S_i = null;
            try {
                // Output Dk_out temporary file in .txt format  (not in .fasta . We create block in .fasta  format after create dictionary)
                PrintWriter Dk_out = null;
                try {
                    Dk_out = new PrintWriter(new FileWriter(path + Dk_Bijection + FS_SEPARATOR + Dk_Bijection + "_" + i + DSK_EXT));
                    if (debug==1) System.out.println("Dk_out = " + path + Dk_Bijection + FS_SEPARATOR + Dk_Bijection + "_" + i + DSK_EXT);

                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                }

                S_i = new BufferedReader(new FileReader(block_ess));
                if (debug==1) System.out.println("Dk_out = " + block_ess);

                // scan S to obtain from the k-mer(S) the respective canons (with the order relation of DSK)
                while ((l = S_i.readLine()) != null) {

                    // skip lines that begin with the character '>'
                    if (!l.startsWith(">")) {

                        // K-mer extraction from the SPSS current sequence
                        for (int j = 0; j < l.length() - k + 1; j++) {
                            String kmer = l.substring(j, j + k);

                            //FIXME: For each x of SPSS S, consider its canonical with respect to DSK order relation, i.e., A<C<T<G .
                            //   see   https://github.com/GATB/dsk#kmers-and-their-reverse-complements

                            Dk_out.println(StringUtils.Canonical_DSK_RelOrd(kmer));
                        }
                    }
                }
                // close Dk_out temporary file
                Dk_out.close();

                // close ESS input file
                S_i.close();

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            //Cleanup Garbage Collector
            System.gc();
        }
    }


    public static void main(String[] args) throws IOException, InterruptedException {

        /*
           0   ->  /home/user/IdeaProjects/TradeOff/experiments/     (initial path to read and write files)
           1   ->  Staph                       (Experiment Dataset Name)
           2   ->  Staph-k4-ESS                (ESS input file in FASTA format)
           3   ->  4                           k value of k-mer tokenization of reads
           4   ->  Staph-k4-Bijection          (Dk output folder and file divided into fragments)
           5   ->  num_blocks                  es. 10
           6   ->  0/1                         no/yes verbose mode
        */

        From_S_and_Fk_rearranged_To_Dk_approx.execute("/home/user/IdeaProjects/TradeOff/experiments/", "Staph","Staph-k4-ESS", 4, "Staph-k4-Bijection", 10, 1);
    }
}
