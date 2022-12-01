import me.lemire.integercompression.*;

import java.io.*;

public class Lemire_Fk_OptPFOR_Decompr_exp {
    
    private Lemire_Fk_OptPFOR_Decompr_exp () {}
    
    protected static void execute (String path, String Fk_To_decompr_in, String scenario, String approx_Suffix) throws InterruptedException, IOException {
        /*
            0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/     (initial path)
            1   ->  Staph-k4-DSK-verbatim                                           (Fk folder with compressed fragments)
            2   ->  b                                                               (a: DPx scenario; b: Gaps scenario)
            3   ->  "" / "_approx"                                                  approx_suffix
         */

        //start time
        long start_time = System.currentTimeMillis();

        String Fk_compr_EXT = ".txt.opt";
        String Fk_decompr_EXT = ".txt.opt.d";
        String nfo_EXT = ".nfo";

        String FS_SEPARATOR = "/";

        int num_rows_file = 0;
        long total_rows = 0;
        int num_fragments = 0;
        int k = 0;

        // Lettura informazioni .nfo per ricavare le informazioni: total_rows, k, num_rows_file, num_fragments
        BufferedReader Header_NFO = null;
        try {

            if (scenario.compareTo("a") == 0) {
                Header_NFO = new BufferedReader(new FileReader(path + Fk_To_decompr_in + FS_SEPARATOR + Fk_To_decompr_in + approx_Suffix + nfo_EXT));
                System.out.println("OptPFOR DECOMPR 001" + "\t" + "Read .nfo file from: " + path + Fk_To_decompr_in + FS_SEPARATOR + Fk_To_decompr_in + approx_Suffix + nfo_EXT);
            }
            else {
                // delete "-Gk" suffix from Fk_To_compr_in for reading .nfo file !!!.
                Header_NFO = new BufferedReader(new FileReader(path + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + FS_SEPARATOR + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + approx_Suffix + nfo_EXT));
                System.out.println("OptPFOR DECOMPR 001" + "\t" + "Read .nfo file from: " + path + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + FS_SEPARATOR + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + approx_Suffix + nfo_EXT);
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Extract info from .nfo file with auxiliary variable for reading lines from file
        String l;

        l = Header_NFO.readLine();
        total_rows = Long.parseLong(l);
        //System.out.println("total_rows = " + total_rows);

        l = Header_NFO.readLine();
        k = Integer.parseInt(l);
        //System.out.println("k = " + k);

        l = Header_NFO.readLine();
        num_rows_file = Integer.parseInt(l);
        //System.out.println("num_rows_file = " + num_rows_file);

        l = Header_NFO.readLine();
        num_fragments = Integer.parseInt(l);
        //System.out.println("num_fragments = " + num_fragments + "\n");

        // Chiusura File
        Header_NFO.close();


        SkippableIntegerCODEC[] codec = {
                new SkippableComposition(new OptPFD(), new VariableByte())
        };


        //TODO: lettura sequenza frammenti compressi
        int i = 0;
        while (i < num_fragments) {
            // se lavoro con i frammenti denotati da _0 .. _(num_fragments - 1), ovvero tutti tranne l'ultimo
            int rows_i;

            if (i != num_fragments - 1) {
            } else {
                rows_i = (int) (total_rows - ((num_fragments - 1) * num_rows_file));
            }
            //System.out.println("Load Fragment " + i + ". Reading: " + rows_i + " rows.");

            //TODO: choose  DPx Scenario (a) or Gaps Scenario (b) --> tolto col refactoring dello script: (Manage_Fa_fragments)
            // String Sub_Type_exp = ( (scenario.compareTo("a") == 0) ? "" : "-Gk");

            //TODO: storage fragment frequencies to int array
            long start_time_storage_frequencies = System.currentTimeMillis();

            //apertura in lettura frammento delle frequenze
            BufferedReader Fk_compr_in = null;
            try {

                if (scenario.compareTo("a") == 0) {
                    Fk_compr_in = new BufferedReader(new FileReader(path + Fk_To_decompr_in + FS_SEPARATOR + Fk_To_decompr_in + approx_Suffix + "_" + i + Fk_compr_EXT));
                    System.out.println("\nOptPFOR Decompr - Read compressed frequencies from: " + path + Fk_To_decompr_in + FS_SEPARATOR + Fk_To_decompr_in + approx_Suffix + "_" + i + Fk_compr_EXT);
                }
                else {
                    // delete "-Gk" suffix from Fk_To_compr_in for reading .nfo file !!!.
                    Fk_compr_in = new BufferedReader(new FileReader(path + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + FS_SEPARATOR + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + approx_Suffix + "-Gk_" + i + Fk_compr_EXT));
                    System.out.println("\nOptPFOR Decompr - Read compressed frequencies from: " + path + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + FS_SEPARATOR + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + approx_Suffix + "-Gk_" + i + Fk_compr_EXT);
                }

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            // variabile indice di conteggio frequenze incontrate
            int countFreq = 0;

            while ((l = Fk_compr_in.readLine()) != null) {
                // incremento il counter delle frequenze incontrate
                countFreq++;
            }

            // Chiusura dei file di input Fk_in
            Fk_compr_in.close();

            //TODO: leggere gli interi direttamente dal file, perchè il numero non è rows_i, ma varia da frammento a frammento

            if (scenario.compareTo("a") == 0) {
                Fk_compr_in = new BufferedReader(new FileReader(path + Fk_To_decompr_in + FS_SEPARATOR + Fk_To_decompr_in + approx_Suffix + "_" + i + Fk_compr_EXT));
                System.out.println("\nOptPFOR Decompr - Read compressed frequencies from: " + path + Fk_To_decompr_in + FS_SEPARATOR + Fk_To_decompr_in + approx_Suffix + "_" + i + Fk_compr_EXT);
            }
            else {
                // delete "-Gk" suffix from Fk_To_compr_in for reading .nfo file !!!.
                Fk_compr_in = new BufferedReader(new FileReader(path + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + FS_SEPARATOR + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + approx_Suffix + "-Gk_" + i + Fk_compr_EXT));
                System.out.println("\nOptPFOR Decompr - Read compressed frequencies from: " + path + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + FS_SEPARATOR + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + approx_Suffix + "-Gk_" + i + Fk_compr_EXT);
            }

            // Allocazione Statica del vettore di interi con cardinalità appena ricevuta dalla scansione del file di input
            int[] data = new int[countFreq];

            int aux_freq;

            //FIXME
            for (int j = 0; j < data.length; ++j) {

                l = Fk_compr_in.readLine();
                aux_freq = Integer.parseInt(l);
                data[j] = aux_freq;
            }
            // Chiusura file Fk_in input fragment file
            Fk_compr_in.close();

            long end_time_storage_frequencies = System.currentTimeMillis();
            long elapsed_time_storage_frequencies = (end_time_storage_frequencies - start_time_storage_frequencies);
            //System.out.println("Elapsed time for storing frequencies: " + elapsed_time_storage_frequencies + " ms");
            // Occupazione in memoria dell'array data
            long dim_data = 16 + ( 4 * data.length );   // 16 represent overheah (12 array header + 4 array length)
            //System.out.println("Memory size occupation of input compressed array data: " + dim_data + " bytes");

            //TODO: OptPFOR decompression phase
            long start_time_decompression_integers = System.currentTimeMillis();

            IntCompressor iic = new IntCompressor(codec[0]);

            //int[] compressed = iic.compress(data);

            int[] recovered = iic.uncompress(data);

            long end_time_decompression_integers = System.currentTimeMillis();
            long elapsed_time_decompression_integers = (end_time_decompression_integers - start_time_decompression_integers);
            //System.out.println("Elapsed time for decompression of integers: " + elapsed_time_decompression_integers + " ms");

            //Test
            //Assert.assertArrayEquals(recovered, data);

            //TODO: storage to file phase
            long start_time_storage_to_file = System.currentTimeMillis();

            // Write Output Lemire_out file in .txt.bic.d format
            PrintWriter Lemire_OptPFOR_out = null;
            try {

                if (scenario.compareTo("a") == 0) {
                    Lemire_OptPFOR_out = new PrintWriter(new FileWriter(path + Fk_To_decompr_in + FS_SEPARATOR + Fk_To_decompr_in + approx_Suffix + "_" + i + Fk_decompr_EXT));
                    System.out.println("OptPFOR Decompr - Write data to: " + path + Fk_To_decompr_in + FS_SEPARATOR + Fk_To_decompr_in + approx_Suffix + "_" + i + Fk_decompr_EXT);
                }
                else {
                    // delete "-Gk" suffix from Fk_To_compr_in path for reading .nfo file !!!.
                    Lemire_OptPFOR_out = new PrintWriter(new FileWriter(path + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + FS_SEPARATOR + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + approx_Suffix + "-Gk_" + i + Fk_decompr_EXT));
                    System.out.println("OptPFOR Decompr - Write data to: " + path + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + FS_SEPARATOR + Fk_To_decompr_in.substring(0, Fk_To_decompr_in.length() -3) + approx_Suffix + "-Gk_" + i + Fk_decompr_EXT);
                }

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            //Salvataggio recovered su file Lemire_OptPFOR_out
            for (int x = 0; x < recovered.length; x++) {
                Lemire_OptPFOR_out.println(recovered[x]);
            }

            // Chiusura file Lemire_OptPFOR_out output file
            Lemire_OptPFOR_out.close();

            long end_time_storage_to_file = System.currentTimeMillis();
            long elapsed_time_storage_to_file = (end_time_storage_to_file - start_time_storage_to_file);
            //System.out.println("Elapsed time for storage to file: " + elapsed_time_storage_to_file + " ms");


            long end_time = System.currentTimeMillis();
            long elapsed_time=(end_time - start_time);
            //System.out.println("Elapsed time: " + elapsed_time + " ms");

            // update fragment counter
            i++;
        }
        // deallocate int arrays recovered and data
        //TODO: free recovered and data arrays from memory
        System.gc();
    }

    public static void main(String[] args) throws IOException, InterruptedException {

         /*
            0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/     (initial path)
            1   ->  Staph-k4-DSK-verbatim                                           (Fk folder with compressed fragments)
            2   ->  b                                                               (a: DPx scenario; b: Gaps scenario)
            3   ->  "" / "_approx"                                                  approx_suffix
         */

        Lemire_Fk_OptPFOR_Decompr_exp.execute("/home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/", "Staph-k4-DSK-verbatim", "b", "");
    }
}
