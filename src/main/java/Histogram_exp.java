import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;

import java.awt.*;
import java.io.*;
import java.nio.file.FileSystems;

public class Histogram_exp {

    private static final String Histo_EXT = ".png";

    private Histogram_exp() {
    }

    protected static void execute(String path, String exp_name, int debug) throws InterruptedException, IOException {

        //start time
        long start_time = System.currentTimeMillis();

        String[] bin_interval_new = {
                "1-5",      // 0
                "6-10",     // 1
                "11-20",    // 2
                "21-30",    // 3
                "31-40",    // 4
                "41-50",    // 5
                "51-60",    // 6
                "61-70",    // 7
                "71-80",    // 8
                "81-90",    // 9
                "91-100",   // 10
                "101-200",  // 11
                "201-300",  // 12
                "301-400",  // 13
                "401-500",  // 14
                "501-600",  // 15
                ">600",     // 16
        };

        String nfo_EXT = ".nfo";
        String Fk_EXT = ".txt";
        String FS_SEPARATOR = FileSystems.getDefault().getSeparator();

        int num_rows_file;
        long total_rows;
        int num_fragments;
        int k;

        //FIXME: Extract informations from .nfo file for retrieval: total_rows, k, num_rows_file, num_fragments
        BufferedReader Header_NFO = null;
        try {

            Header_NFO = new BufferedReader(new FileReader(path + exp_name + FS_SEPARATOR + exp_name + nfo_EXT));
            if (debug==1) System.out.println("Histogram creation." + "\t" + "Read .nfo file from: " + path + exp_name + FS_SEPARATOR + exp_name + nfo_EXT);

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        //FIXME: Extract info from .nfo file with auxiliary variable for reading lines from file
        String l;

        assert Header_NFO != null;

        l = Header_NFO.readLine();
        total_rows = Long.parseLong(l);

        l = Header_NFO.readLine();
        k = Integer.parseInt(l);

        l = Header_NFO.readLine();
        num_rows_file = Integer.parseInt(l);

        l = Header_NFO.readLine();
        num_fragments = Integer.parseInt(l);

        // Chiusura File
        Header_NFO.close();


        // Reading fragment sequence for construction of histogram data set

        // fragment counter
        int i = 0;

        // Auxiliary variable for managing fragment frequencies
        int aux_freq;

        //FIXME: new interval
        // Static allocation of the binning vector containing the intervals:
        // 1-5; 6-10; 11-20; 21-30; 31-40; 41-50; 51-60; 61-70; 71-80; 81-90; 91-100;
        // 101-200; 201-300; 301-400; 401-500; 501-600; >600
        // for a total of 17 bins
        double[] bins = new double[17];

        //initialize bins values
        for (int j = 0; j < bins.length; j++) {
            bins[i] = 0;
        }

        while (i < num_fragments) {

            // I work with the fragments denoted by _0 ... _(num_fragments - 1), that is, all but the last one.
            int rows_i;

            //FIXME: Compute how many rows are present in current fragment
            if (i != num_fragments - 1) rows_i = num_rows_file;
            else {
                rows_i = (int) (total_rows - ((num_fragments - 1) * num_rows_file));
            }
            if (debug==1) System.out.println("Load Fragment " + i + ". Reading: " + rows_i + " rows.");

            //FIXME: Opening Fragment of Frequencies
            BufferedReader Fk_fragment_in = null;
            try {
                Fk_fragment_in = new BufferedReader(new FileReader(path + exp_name + FS_SEPARATOR + exp_name + "_" + i + Fk_EXT));
                if (debug==1) System.out.println("Fk current fragment. - Read from: " + path + exp_name + FS_SEPARATOR + exp_name + "_" + i + Fk_EXT);

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            // FIXME: frequency occurrences computing and updating
            for (int j = 0; j < rows_i; ++j) {

                assert Fk_fragment_in != null;

                l = Fk_fragment_in.readLine();
                aux_freq = Integer.parseInt(l);

                if (isBetween(aux_freq, 1, 5)) {
                    bins[0]++;
                } else if (isBetween(aux_freq, 6, 10)) {
                    bins[1]++;
                } else if (isBetween(aux_freq, 11, 20)) {
                    bins[2]++;
                } else if (isBetween(aux_freq, 21, 30)) {
                    bins[3]++;
                } else if (isBetween(aux_freq, 31, 40)) {
                    bins[4]++;
                } else if (isBetween(aux_freq, 41, 50)) {
                    bins[5]++;
                } else if (isBetween(aux_freq, 51, 60)) {
                    bins[6]++;
                } else if (isBetween(aux_freq, 61, 70)) {
                    bins[7]++;
                } else if (isBetween(aux_freq, 71, 80)) {
                    bins[8]++;
                } else if (isBetween(aux_freq, 81, 90)) {
                    bins[9]++;
                } else if (isBetween(aux_freq, 91, 100)) {
                    bins[10]++;
                } else if (isBetween(aux_freq, 101, 200)) {
                    bins[11]++;
                } else if (isBetween(aux_freq, 201, 300)) {
                    bins[12]++;
                } else if (isBetween(aux_freq, 301, 400)) {
                    bins[13]++;
                } else if (isBetween(aux_freq, 401, 500)) {
                    bins[14]++;
                } else if (isBetween(aux_freq, 501, 600)) {
                    bins[15]++;
                } else if (isBetween(aux_freq, 601, Integer.MAX_VALUE)) {
                    bins[16]++;
                }
            }

            //update fragment counter
            i++;

        }


        // FIXME: Write bins info to histo.txt file
        PrintWriter histo_out = null;
        try {
            histo_out = new PrintWriter(new FileWriter(path + exp_name + FS_SEPARATOR + exp_name + "_histo" + Fk_EXT));
            if (debug==1) System.out.println("Write Histo file to: " + path + exp_name + FS_SEPARATOR + exp_name + "_histo" + Fk_EXT);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        assert histo_out != null;

        histo_out.println("interval, occurrencies, percent");
        for (int j = 0; j < bins.length; j++) {
            //"%d%n" + "\t" + "%.0f%n", j+1, bins[j]
            histo_out.println(bin_interval_new[j] + ", " + bins[j] + ", " + bins[j]/total_rows);
        }
        histo_out.close();


        //FIXME: BarChart creation
        //dataset1
        DefaultCategoryDataset dataset1 = new DefaultCategoryDataset();

        for (int j = 0; j < bins.length; j++) {
            dataset1.setValue(bins[j]/total_rows, "occurrences", bin_interval_new[j]);
        }

        String plotTitle = "Frequency Distribution";
        String xaxis = "bins";
        String yaxis = "frequency occurences";
        PlotOrientation orientation = PlotOrientation.VERTICAL;
        boolean show = false;
        boolean toolTips = false;
        boolean urls = false;
        JFreeChart barchart = ChartFactory.createBarChart( plotTitle, xaxis, yaxis, dataset1, orientation, show, toolTips, urls);
        barchart.setBackgroundPaint(Color.white);


        //FIXME: save PNG image to experiment folder
        int width = 1980;
        int heigth = 1080;
        try {
            ChartUtils.saveChartAsPNG(new File(path + exp_name + FS_SEPARATOR + exp_name + "_histo" + Histo_EXT), barchart, width, heigth);
        } catch (IOException e) {
            e.printStackTrace();
        }

        //stop time
        long time_elapsed = System.currentTimeMillis() - start_time;

        //DEBUG
        if (debug==1) System.out.println("time_elapsed = " + time_elapsed);
    }

    public static boolean isBetween (int x, int lower, int upper) {
        return lower <= x && x <= upper;
    }

    public static void main(String[] args) throws IOException, InterruptedException {

        /* Example:

            0   ->  /home/user/IdeaProjects/TradeOff/experiments/StaphAU/k4/                (initial path)
            1   ->  StaphAU-DSK-k4-verbatim                                                 (Fk folder with fragments)
            2   ->  0/1                                                                     no/yes verbose mode

         */

        Histogram_exp.execute("/home/user/IdeaProjects/TradeOff/experiments/AssPlants/k16/", "AssPlants-DSK-k16-verbatim", 1);
    }
}
