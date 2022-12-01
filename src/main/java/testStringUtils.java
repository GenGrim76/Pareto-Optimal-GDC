import java.io.IOException;

public class testStringUtils {

    public static void main(String[] args) throws IOException {

        // kmer = "GGAC";       //DSK behaviour
        // rev_compl_kmer = "GTCC";
        // can_DSK = "GTCC", why in pos.=1 T<G .

        String kmer = "GGAC";

        String canonical = StringUtils.Canonical_DSK_RelOrd_speedup(kmer);

        System.out.println("Canonical(DSK) of " + kmer + " is: " + canonical);
        System.out.println("");
   }
}
