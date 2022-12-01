import java.util.Locale;

public class StringUtils {

    //Costruttore private: non sarà possibile istanziare oggetti di questa classe
    private StringUtils() {
    }

    public static String ReverseComplement(String kmer) {

        int len = kmer.length();

        // transform the input kmer into StringBuilder
        char C_kmer[] = kmer.toCharArray();

        // Auxiliary StringBuilder which will contain the reverse Complement of kmer
        char C_RCkmer[] = new char[len];

        //Scan k-mer Left to Right
        for (int i = 0; i < len; i++) {
            switch (C_kmer[i]) {
                case 'A':
                    C_RCkmer[len - i - 1] = 'T';
                    break;
                case 'C':
                    C_RCkmer[len - i - 1] = 'G';
                    break;
                case 'G':
                    C_RCkmer[len - i - 1] = 'C';
                    break;
                case 'T':
                    C_RCkmer[len - i - 1] = 'A';
                    break;
            }
        }
        return String.valueOf(C_RCkmer);
    }

    public static String Canonical(String kmer) {

        String canonical = ReverseComplement(kmer);

        if (kmer.compareTo(canonical) < 0)
            return kmer;
        else
            return canonical;

    }

    public static String Canonical_DSK_RelOrd(String kmer) {

        // https://github.com/GATB/dsk#kmers-and-their-reverse-complements
        // https://qastack.it/software/127639/why-do-some-sorting-methods-sort-by-1-10-2-3

        String canonical = ReverseComplement(kmer);

        //DEBUG
        // System.out.println("Reverse Complement of " + kmer + " is: " + canonical);

        int k = kmer.length();

        //DEBUG
        // System.out.println("Length of kmer is: " + k);

        // control, position by position, for sorting according to order relation: (  A < C < T < G  )
        // consider A with 0, C with 1, T with 2, G with 3.
        // nucleotide base conversions A, C, G, T with 0, 1, 3, 2


        // transform the input kmer and canonical into StringBuilder
        char C_kmer[] = new char[k];
        char C_canonical[] = new char[k];

        C_kmer = kmer.toCharArray();

        //DEBUG
        // System.out.print(C_kmer);
        // System.out.println(" con lunghezza: " + C_kmer.length);

        C_canonical = canonical.toCharArray();

        //DEBUG
        // System.out.print(C_canonical);
        // System.out.println(" witrh length: " + C_canonical.length);

        // Scan k-mer Left to Right mode (Most Significative to Least Significative Nucleotide)
        for (int i = 0; i < k; i++) {
            switch (C_kmer[i]) {
                case 'A':
                    C_kmer[i] = '0';
                    break;
                case 'C':
                    C_kmer[i] = '1';
                    break;
                case 'T':
                    C_kmer[i] = '2';
                    break;
                case 'G':
                    C_kmer[i] = '3';
                    break;
            }
        }

        //DEBUG
        // System.out.print(C_kmer);
        // System.out.println(" with length: " + C_kmer.length);


        //Scan canonical Left to Right (Most Significative to Least Significative Nucleotide)
        for (int i = 0; i < k; i++) {
            switch (C_canonical[i]) {
                case 'A':
                    C_canonical[i] = '0';
                    break;
                case 'C':
                    C_canonical[i] = '1';
                    break;
                case 'T':
                    C_canonical[i] = '2';
                    break;
                case 'G':
                    C_canonical[i] = '3';
                    break;
            }
        }

        //DEBUG
        // System.out.print(C_canonical);
        // System.out.println(" with length: " + C_canonical.length);


        // long value_kmer = Long.parseLong(String.valueOf(C_kmer));
        String str_kmer = String.valueOf(C_kmer);


        // long value_canonical = Long.parseLong(String.valueOf(C_canonical));
        String str_canonical = String.valueOf(C_canonical);


        //DEBUG
        // System.out.println("kmer and its RC: " + kmer + " " + canonical);

        // System.out.println("Transform :    " + value_kmer + " " + value_canonical);
        // System.out.println("Transform :    " + str_kmer + " " + str_canonical);

        // System.out.print("min (rel. DSK): ");


        // if (value_kmer <= value_canonical) {
        if (str_kmer.compareTo(str_canonical) <= 0) {
            // System.out.print(value_kmer);
            // System.out.println("");

            // return String.valueOf(kmer);
            return kmer;
        } else {
            // System.out.print(value_canonical);
            // System.out.println("");

            // return String.valueOf(canonical);
            return canonical;
        }
    }

    public static String Canonical_DSK_RelOrd_speedup(String kmer) {

        //TODO: A, C, T, G DSK relation order
        String reverse_kmer = StringUtils.ReverseComplement(kmer);

        int continua = 1;
        int index = 0;
        int scelta = -1;
        String aux;

        while (index <kmer.length() && (scelta == -1) ) {

            char c1 = kmer.charAt(index);
            char c2 = reverse_kmer.charAt(index);

            //System.out.println("c1 = " + c1);
            //System.out.println("c2 = " + c2);

            if (c1 != c2) {
                if ( (c1 == 'A') && ( (c2 == 'C') || (c2=='T') || (c2=='G') ) ) {
                    //System.out.println("c1=" + c1 + " - index=" + index);
                    scelta = 0;
                    continua = 0;
                }

                if ( (c1 == 'C') && ( (c2 == 'T') || (c2=='G') ) ) {
                    //System.out.println("c1=" + c1 + " - index=" + index);
                    scelta = 0;
                }
                else scelta =1;

                if ((c1 == 'T') && (c2 == 'G')) {
                    //System.out.println("c1=" + c1 + " - index=" + index);
                    scelta = 0;
                }
                else scelta = 1;

                if  (c1 == 'G') {
                    //System.out.println("c1=" + c1 + " - index=" + index);
                    scelta = 1;
                }
                else scelta =0;
            }
            else { //c1==c2
                index++;
            }
        }
        if (scelta == 0) {
            return kmer;
        } else if (scelta == 1) {
            return reverse_kmer;
        } else return kmer; // //scelta == -1, I choose arbitrarily, they represent the same k-mer
    }


}
