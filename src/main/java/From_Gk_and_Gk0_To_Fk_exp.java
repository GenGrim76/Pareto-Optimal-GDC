import java.io.*;
import java.nio.file.FileSystems;

public class From_Gk_and_Gk0_To_Fk_exp {
    
    private From_Gk_and_Gk0_To_Fk_exp () {}
    
    protected static void execute (String path, String Gk_folder_in, int approx, int debug) throws InterruptedException, IOException {

        /*
            0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/    (initial path to read and write files)
            1   ->  Staph-k4-DSK-verbatim                                          (Gk input file divided into fragments)
            2   ->  0/1                                                            exact/approximate solution
            3   ->  0/1                                                            no/yes verbose mode
        */

        String Fk_out_EXT = ".txt";
        String Gk_EXT = ".txt";
        String Gk0_EXT = ".txt";
        String nfo_EXT = ".nfo";

        String FS_SEPARATOR = FileSystems.getDefault().getSeparator();

        // Approximate Solution Check
        String approx_Suffix = null;
        if (approx == 0) {
            approx_Suffix = "";
        }
        else {
            approx_Suffix = "_approx";
        }


        int num_rows_file = 0;
        long total_rows = 0;
        int num_fragments = 0;
        int k = 0;


        // Lettura informazioni .nfo per ricavare le informazioni: total_rows, k, num_rows_file, num_fragments
        BufferedReader Header_NFO = null;
        try {
            Header_NFO = new BufferedReader(new FileReader(path + Gk_folder_in + FS_SEPARATOR + Gk_folder_in + approx_Suffix + nfo_EXT));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Extract info from .nfo file
        String l;

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

        //TODO: lettura first_frequency from Gk0 file
        // Read Gk file info Gk0_in .txt format
        BufferedReader Gk0_in = null;
        try {
            Gk0_in = new BufferedReader(new FileReader(path + Gk_folder_in + FS_SEPARATOR + Gk_folder_in + approx_Suffix + "-Gk0" + Gk0_EXT));
            if (debug==1) System.out.println("Read info from: " + path + Gk_folder_in + FS_SEPARATOR + Gk_folder_in + approx_Suffix + "-Gk0" + Gk0_EXT + "\n");

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        // Scan Gk0 to get the first frequency from the file Gk0_in
        long Gk0_freq = Long.parseLong(Gk0_in.readLine());
        // GK0_in File Closure
        Gk0_in.close();


        //FIXME: reading fragments sequence
        int i = 0;
        while (i < num_fragments) {

            // I work with the fragments denoted by _0 ... _(num_fragments - 1), that is, all but the last one.
            int rows_i;

            if (i != num_fragments - 1) rows_i = num_rows_file;
            else {
                rows_i = (int) (total_rows - ( (num_fragments - 1) * num_rows_file));
            }


            // Read Gk_i file in .txt format
            BufferedReader Gk_in = null;
            try {
                Gk_in = new BufferedReader(new FileReader(path + Gk_folder_in + FS_SEPARATOR + Gk_folder_in + approx_Suffix + "-Gk_" + i + Gk_EXT));
                if (debug==1) System.out.println("Read from: " + path + Gk_folder_in + FS_SEPARATOR + Gk_folder_in + approx_Suffix + "-Gk_" + i + Gk_EXT);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }


            // Write Output Fk_out_i file in .txt format
            PrintWriter Fk_out = null;
            try {
                Fk_out = new PrintWriter(new FileWriter(path + Gk_folder_in + FS_SEPARATOR + Gk_folder_in + approx_Suffix + "-Fk-Post_" + i + Fk_out_EXT));
                if (debug==1) System.out.println("Write to: " + path + Gk_folder_in + FS_SEPARATOR + Gk_folder_in + approx_Suffix + "-Fk-Post_" + i + Fk_out_EXT);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            // Gk scan for reading the frequency gaps from the Gk_in file for reconstruction of the Fk file fragments

            // Gap count index variable encountered for which to sum Gk0
            long countGap = 0;

            long current_i = Gk0_freq;

            for (int j = 0; j < rows_i; j++) {

                l = Gk_in.readLine();

                // update the counter of the frequencies encountered
                countGap++;

                long gap_i = Long.parseLong(l);

                current_i += gap_i;

                Fk_out.println(current_i);
            }

            // I save the last frequency of the fragment for the first gap of the next possible fragment.
            Gk0_freq = current_i;

            // update fragment counter
            i++;

            // Closing the input files Gk_in and the output file Fk_out
            Fk_out.close();
            Gk_in.close();

            // for removing use (i-1) frament index
            int pos = i-1;

            // remove i_th fragment Fk-Post files, after identity check
            File file = new File(path + Gk_folder_in + FS_SEPARATOR + Gk_folder_in + approx_Suffix + "-Fk-Post_" + pos + Fk_out_EXT);
            if (file.exists() && file.isFile()) {
                file.delete();
            } else {
                System.out.println(path + Gk_folder_in + FS_SEPARATOR + Gk_folder_in + approx_Suffix + "-Fk-Post_" + pos + Fk_out_EXT + " not a file or non-existent file !!!.");
            }
        }

        // Cleanup Garbage Collector
        System.gc();
    }

    public static void main(String[] args) throws IOException, InterruptedException {

        /*
            0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/    (initial path to read and write files)
            1   ->  Staph-k4-DSK-verbatim                                          (Gk input file divided into fragments)
            2   ->  0/1                                                            exact/approximate solution
            3   ->  0/1                                                            no/yes verbose mode
        */

        From_Gk_and_Gk0_To_Fk_exp.execute("/home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/", "Staph-k4-DSK-verbatim", 0, 1);
    }
}
