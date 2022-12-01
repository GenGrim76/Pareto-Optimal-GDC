import java.io.*;
import java.nio.file.FileSystems;
import java.util.Scanner;

public class From_Fk_To_Gk_and_Gk0_exp {

    private From_Fk_To_Gk_and_Gk0_exp () {}

    protected static void execute (String path, String Fk_folder_in, int approx, int debug) throws InterruptedException, IOException {

        /*
            0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/    (initial path to read and write files)
            1   ->  Staph-k4-DSK-verbatim                                          (Fk input file divided into fragments)
            2   ->  0/1                                                            exact/approximate solution
            3   ->  0/1                                                            no/yes verbose mode

        */

        String Fk_EXT = ".txt";
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

        // auxiliary variables
        String l;
        long aux_freq = 0, gap_i = 0, frequency_i = 0;

        int num_rows_file = 0;
        long total_rows = 0;
        int num_fragments = 0;
        int k = 0;

        // Reading .nfo information to extract information: total_rows, k, num_rows_file, num_fragments
        BufferedReader Header_NFO = null;
        try {
            Header_NFO = new BufferedReader(new FileReader(path + Fk_folder_in + FS_SEPARATOR + Fk_folder_in + approx_Suffix + nfo_EXT));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        //TODO: Extract info from .nfo file
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


        // FIXME: Opening Fragment file #0
        // Read current fragment Fk_0 file in .txt format
        BufferedReader Fk_in = null;
        try {
            Fk_in = new BufferedReader(new FileReader(path + Fk_folder_in + FS_SEPARATOR + Fk_folder_in + approx_Suffix + "_0" + Fk_EXT));
            if (debug==1) System.out.println("Read: " + path + Fk_folder_in + FS_SEPARATOR + Fk_folder_in + "_0" + Fk_EXT);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Write Output fragment Gk_out file in .txt format
        PrintWriter Gk_out = null;
        try {
            Gk_out = new PrintWriter(new FileWriter(path + Fk_folder_in + FS_SEPARATOR + Fk_folder_in + approx_Suffix + "-Gk" + "_0" + Gk_EXT));
            if (debug==1) System.out.println("Write: " + path + Fk_folder_in + FS_SEPARATOR + Fk_folder_in + approx_Suffix + "-Gk" + "_0" + Gk_EXT);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Opening Gk0 file for saving the first frequency of the fragment Fk_0
        // Write Output Gk0_out file in .txt format
        PrintWriter Gk0_out = null;
        try {
            Gk0_out = new PrintWriter(new FileWriter(path + Fk_folder_in + FS_SEPARATOR + Fk_folder_in + approx_Suffix + "-Gk0" + Gk0_EXT));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        //FIXME: Scan first fragment Fk file to get the first frequency from the file
        long first_frequency = Long.parseLong(Fk_in.readLine());

        //Write first frequency to Gk0_out output file
        Gk0_out.println(first_frequency);
        // close Gk0_out output file
        Gk0_out.close();

        // Set first frequency in Gk to 0, represent first_frequency in Fk.
        Gk_out.println(0);

        int i = 0;
        while (i < num_fragments) {

            // I work with the fragments denoted by _0 ... _(num_fragments - 1), that is, all but the last one.
            int rows_i;
            int in_j;

            if (i != num_fragments - 1) rows_i = num_rows_file;
            else {
                rows_i = (int) (total_rows - ( (num_fragments - 1) * num_rows_file));
            }

            //System.out.println("Load Fragment " + i + ". Reading: " + rows_i + " rows.");

            if (i == 0) {
                // already a line from the Fk fragment I read it
                in_j = 1;

            }
            else {
                // I have not yet read the first line of the fragment
                in_j = 0;
            }

            // auxiliary variable for computing gaps
            aux_freq = first_frequency;

            for (int j = in_j; j < rows_i; j++) {

                l = Fk_in.readLine();

                frequency_i = Long.parseLong(l);

                gap_i = frequency_i - aux_freq;

                Gk_out.println(gap_i);

                aux_freq = frequency_i;
            }

            // I save the last frequency of the fragment for the first gap of the next possible fragment
            first_frequency = aux_freq;

            // Fk and Gk current fragment closure
            Gk_out.close();
            Fk_in.close();

            // update fragment counter
            i++;

            if ( i != (num_fragments) ) {
                // FIXME: Reopen frequency Fk and Gk fragments file
                Fk_in = new BufferedReader(new FileReader(path + Fk_folder_in + FS_SEPARATOR + Fk_folder_in + approx_Suffix + "_" + i + Fk_EXT));
                if (debug==1) System.out.println("Read:  " + path + Fk_folder_in + FS_SEPARATOR + Fk_folder_in + approx_Suffix + "_" + i + Fk_EXT);

                Gk_out = new PrintWriter(new FileWriter(path + Fk_folder_in + FS_SEPARATOR + Fk_folder_in + approx_Suffix + "-Gk" + "_" + i + Gk_EXT));
                if (debug==1) System.out.println("Write:  " + path + Fk_folder_in + FS_SEPARATOR + Fk_folder_in + approx_Suffix + "-Gk" + "_" + i + Gk_EXT);
            }
        }

        // Cleanup Garbage Collector
        System.gc();
    }

    public static void main(String[] args) throws IOException, InterruptedException {

        /*
            0   ->  /home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/    (initial path to read and write files)
            1   ->  Staph-k4-DSK-verbatim                                          (Fk input file divided into fragments)
            2   ->  0/1                                                            exact/approximate solution
            3   ->  0/1                                                            no/yes verbose mode
        */

        From_Fk_To_Gk_and_Gk0_exp.execute("/home/ubuntu-fs-server/IdeaProjects/TradeOff/exp-fragments/", "Staph-k4-DSK-verbatim", 0, 1);

    }
}
