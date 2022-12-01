import me.lemire.integercompression.differential.IntegratedIntCompressor;

import java.io.*;
import java.nio.file.FileSystems;

public class Lemire_Fk_BIC_Compr_exp {

    private Lemire_Fk_BIC_Compr_exp () {}

    protected static void execute (String path, String Fk_To_compr_in, String scenario, String approx_Suffix, int debug) throws InterruptedException, IOException {
        /*
            0   ->  /home/user/IdeaProjects/TradeOff/exp-fragments/                 (initial path)
            1   ->  Staph-k4-DSK-verbatim                                           (Fk folder with fragments)
            2   ->  b                                                               (a: DPx scenario; b: Gaps scenario)
            3   ->  "" / "_approx"                                                  approx_suffix
            4   ->  0/1                                                             no/yes verbose mode
         */

        //start time
        long start_time = System.currentTimeMillis();

        String Fk_EXT = ".txt";
        String Fk_compr_EXT = ".txt.bic";
        String nfo_EXT = ".nfo";

        // File System File Separator
        String FS_SEPARATOR = FileSystems.getDefault().getSeparator();

        int num_rows_file = 0;
        long total_rows = 0;
        int num_fragments = 0;
        int k = 0;


        // Reading .nfo information: total_rows, k, num_rows_files, num_fragments
        BufferedReader Header_NFO = null;
        try {
            if (scenario.compareTo("a") == 0) {
                Header_NFO = new BufferedReader(new FileReader(path + Fk_To_compr_in + FS_SEPARATOR + Fk_To_compr_in + approx_Suffix + nfo_EXT));
                if (debug==1) System.out.println("\nBIC COMPR" + "\t" + "Read .nfo file from: " + path + Fk_To_compr_in + FS_SEPARATOR + Fk_To_compr_in + approx_Suffix + nfo_EXT);
            }
            else {
                // delete "-Gk" suffix from Fk_To_compr_in for reading .nfo file !!!.
                Header_NFO = new BufferedReader(new FileReader(path + Fk_To_compr_in.substring(0, Fk_To_compr_in.length() -3) + FS_SEPARATOR + Fk_To_compr_in.substring(0, Fk_To_compr_in.length() -3) + approx_Suffix + nfo_EXT));
                if (debug==1) System.out.println("\nBIC COMPR" + "\t" + "Read .nfo file from: " + path + Fk_To_compr_in.substring(0, Fk_To_compr_in.length() -3) + FS_SEPARATOR + Fk_To_compr_in.substring(0, Fk_To_compr_in.length() -3) + approx_Suffix + nfo_EXT);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Extract info from .nfo file with auxiliary variable for reading lines from file
        String l;

        l = Header_NFO.readLine();
        total_rows = Long.parseLong(l);

        l = Header_NFO.readLine();
        k = Integer.parseInt(l);

        l = Header_NFO.readLine();
        num_rows_file = Integer.parseInt(l);

        l = Header_NFO.readLine();
        num_fragments = Integer.parseInt(l);

        // Closing Header File
        Header_NFO.close();


        //FIXME: Reading Fragment Sequence
        int i = 0;
        while (i < num_fragments) {
            // work with the fragments denoted by _0 ... _(num_fragments - 1), i.e. all except the last one
            int rows_i;

            if (i != num_fragments - 1) rows_i = num_rows_file;
            else {
                rows_i = (int) (total_rows - ( (num_fragments - 1) * num_rows_file));
            }
            //System.out.println("Load Fragment " + i + ". Reading: " + rows_i + " rows.");

            //FIXME: choose  DPx Scenario (a) or Gaps Scenario (b) --> tolto col refactoring dello script: (Manage_Fa_fragments)
            // String Sub_Type_exp = ( (scenario.compareTo("a") == 0) ? "" : "-Gk");


            //FIXME: storage fragment frequencies to int array
            long start_time_storage_frequencies = System.currentTimeMillis();

            // opening in reading of the frequency fragment
            BufferedReader Fk_in = null;
            try {

                if (scenario.compareTo("a") == 0) {
                    Fk_in = new BufferedReader(new FileReader(path + Fk_To_compr_in + FS_SEPARATOR + Fk_To_compr_in + approx_Suffix + "_" + i + Fk_EXT));
                    if (debug==1) System.out.println("\nBIC COMPR - Read frequencies from: " + path + Fk_To_compr_in + FS_SEPARATOR + Fk_To_compr_in + approx_Suffix + "_" + i + Fk_EXT);
                }
                else {
                    // delete "-Gk" suffix from Fk_To_compr_in for manage Fk file !!!.
                    Fk_in = new BufferedReader(new FileReader(path + Fk_To_compr_in.substring(0, Fk_To_compr_in.length() -3) + FS_SEPARATOR + Fk_To_compr_in.substring(0, Fk_To_compr_in.length() -3) + approx_Suffix + "-Gk_" + i + Fk_EXT));
                    if (debug==1) System.out.println("\nBIC COMPR - Read frequencies from: " + path + Fk_To_compr_in.substring(0, Fk_To_compr_in.length() -3) + FS_SEPARATOR + Fk_To_compr_in.substring(0, Fk_To_compr_in.length() -3) + approx_Suffix + "-Gk_" + i + Fk_EXT);
                }

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            // Static allocation of the vector of integers with cardinality equal to rows_i
            int[] data = new int[rows_i];

            int aux_freq;

            for (int j = 0; j < data.length; ++j) {

                l = Fk_in.readLine();
                aux_freq = Integer.parseInt(l);
                data[j] = aux_freq;
            }
            // Fk_in input fragment file closure
            Fk_in.close();

            long end_time_storage_frequencies = System.currentTimeMillis();
            long elapsed_time_storage_frequencies = (end_time_storage_frequencies - start_time_storage_frequencies);
            //System.out.println("Elapsed time for storing frequencies into int array: " + elapsed_time_storage_frequencies + " ms");
            // Memory occupation of the data array
            long dim_data = 16 + ( 4 * data.length );   // 16 represent overhead (12 array header + 4 array length)
            //System.out.println("Memory size occupation of array data: " + dim_data + " bytes");

            //FIXME: compression phase
            long start_time_compression_integers = System.currentTimeMillis();

            //BIC
            IntegratedIntCompressor iic = new IntegratedIntCompressor();
            int[] compressed = iic.compress(data);

            // Occupazione in memoria dell'array compressed
            long dim_compressed = 16 + ( 4 * compressed.length );
            //System.out.println("Memory size occupation of array compressed: " + dim_compressed + " bytes");

            long end_time_compression_integers = System.currentTimeMillis();
            long elapsed_time_compression_integers = (end_time_compression_integers - start_time_compression_integers);
            //System.out.println("Elapsed time for compression of integers: " + elapsed_time_compression_integers + " ms");

            //DEBUG
            //int[] recovered = iic.uncompress(compressed);
            //Test
            //Assert.assertArrayEquals(recovered, data);
            //DEBUG
            //System.out.println("Dimensione [#integers] Integer Array Data Uncompressed: " + data.length);
            //System.out.println("Dimensione [#integers] Integer array Data Compressed (BIC): " + compressed.length);
            //System.out.println("Dimensione [#integers] Integer array Data Decompressed (BIC): " + recovered.length);

            //FIXME: storage to file phase
            long start_time_storage_to_file = System.currentTimeMillis();

            // Write Output Lemire_out file in .txt.bic format
            PrintWriter Lemire_BIC_out = null;
            try {

                if (scenario.compareTo("a") == 0) {
                    Lemire_BIC_out = new PrintWriter(new FileWriter(path + Fk_To_compr_in + FS_SEPARATOR + Fk_To_compr_in + approx_Suffix + "_" + i + Fk_compr_EXT));
                    if (debug==1) System.out.println("\nBIC COMPR" + "\t" + "Write data to: " + path + Fk_To_compr_in + FS_SEPARATOR + Fk_To_compr_in + approx_Suffix + "_" + i + Fk_compr_EXT);
                }
                else {
                    // delete "-Gk" suffix from Fk_To_compr_in path for reading .nfo file !!!.
                    Lemire_BIC_out = new PrintWriter(new FileWriter(path + Fk_To_compr_in.substring(0, Fk_To_compr_in.length() -3) + FS_SEPARATOR + Fk_To_compr_in.substring(0, Fk_To_compr_in.length() -3) + approx_Suffix + "-Gk_" + i + Fk_compr_EXT));
                    if (debug==1) System.out.println("\nBIC COMPR" + "\t" + "Write data to: " + path + Fk_To_compr_in.substring(0, Fk_To_compr_in.length() -3) + FS_SEPARATOR + Fk_To_compr_in.substring(0, Fk_To_compr_in.length() -3) + approx_Suffix + "-Gk_" + i + Fk_compr_EXT);
                }

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            //Storage compressed data on Lemire_BIC_out file
            for (int x = 0; x < compressed.length; x++) {
                Lemire_BIC_out.println(compressed[x]);
            }
            // Closing Lemire_BIC_out output file
            Lemire_BIC_out.close();

            long end_time_storage_to_file = System.currentTimeMillis();
            long elapsed_time_storage_to_file = (end_time_storage_to_file - start_time_storage_to_file);
            //System.out.println("Elapsed time for storage to file: " + elapsed_time_storage_to_file + " ms");

            long end_time = System.currentTimeMillis();
            long elapsed_time=(end_time - start_time);
            //System.out.println("Elapsed time: " + elapsed_time + " ms" + "\n\n");

            // update fragment counter
            i++;
        }
        // deallocate int arrays compressed and data
        //FIXME: free compressed and data arrays from memory
        System.gc();
    }

    public static void main(String[] args) throws IOException, InterruptedException {

         /*
            0   ->  /home/user/IdeaProjects/TradeOff/exp-fragments/                 (initial path)
            1   ->  Staph-k4-DSK-verbatim                                           (Fk folder with fragments)
            2   ->  b                                                               (a: DPx scenario; b: Gaps scenario)
            3   ->  "" / "_approx"                                                  approx_suffix
            4   ->  0/1                                                             no/yes verbose mode
         */

        Lemire_Fk_BIC_Compr_exp.execute("/home/user/IdeaProjects/TradeOff/exp-fragments/", "Staph-k4-DSK-verbatim", "b", "", 1);

    }
}
