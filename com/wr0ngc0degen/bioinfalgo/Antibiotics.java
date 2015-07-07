package com.wr0ngc0degen.bioinfalgo;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by Alena on 08.11.2014.
 */
public class Antibiotics
{
    public static HashMap<String, String> dnaToAA = new HashMap<String, String>();
    public static HashMap<String, Integer> aaMasses = new HashMap<>();

    static
    {
        aaMasses.put("G", 57);
        aaMasses.put("A", 71);
        aaMasses.put("S", 87);
        aaMasses.put("P", 97);
        aaMasses.put("V", 99);
        aaMasses.put("T", 101);
        aaMasses.put("C", 103);
        aaMasses.put("I", 113);
        aaMasses.put("L", 113);
        aaMasses.put("N", 114);
        aaMasses.put("D", 115);
        aaMasses.put("K", 128);
        aaMasses.put("Q", 128);
        aaMasses.put("E", 129);
        aaMasses.put("M", 131);
        aaMasses.put("H", 137);
        aaMasses.put("F", 147);
        aaMasses.put("R", 156);
        aaMasses.put("Y", 163);
        aaMasses.put("W", 186);

        dnaToAA.put("AAA", "K");
        dnaToAA.put("AAC", "N");
        dnaToAA.put("AAG", "K");
        dnaToAA.put("AAU", "N");
        dnaToAA.put("ACA", "T");
        dnaToAA.put("ACC", "T");
        dnaToAA.put("ACG", "T");
        dnaToAA.put("ACU", "T");
        dnaToAA.put("AGA", "R");
        dnaToAA.put("AGC", "S");
        dnaToAA.put("AGG", "R");
        dnaToAA.put("AGU", "S");
        dnaToAA.put("AUA", "I");
        dnaToAA.put("AUC", "I");
        dnaToAA.put("AUG", "M");
        dnaToAA.put("AUU", "I");
        dnaToAA.put("CAA", "Q");
        dnaToAA.put("CAC", "H");
        dnaToAA.put("CAG", "Q");
        dnaToAA.put("CAU", "H");
        dnaToAA.put("CCA", "P");
        dnaToAA.put("CCC", "P");
        dnaToAA.put("CCG", "P");
        dnaToAA.put("CCU", "P");
        dnaToAA.put("CGA", "R");
        dnaToAA.put("CGC", "R");
        dnaToAA.put("CGG", "R");
        dnaToAA.put("CGU", "R");
        dnaToAA.put("CUA", "L");
        dnaToAA.put("CUC", "L");
        dnaToAA.put("CUG", "L");
        dnaToAA.put("CUU", "L");
        dnaToAA.put("GAA", "E");
        dnaToAA.put("GAC", "D");
        dnaToAA.put("GAG", "E");
        dnaToAA.put("GAU", "D");
        dnaToAA.put("GCA", "A");
        dnaToAA.put("GCC", "A");
        dnaToAA.put("GCG", "A");
        dnaToAA.put("GCU", "A");
        dnaToAA.put("GGA", "G");
        dnaToAA.put("GGC", "G");
        dnaToAA.put("GGG", "G");
        dnaToAA.put("GGU", "G");
        dnaToAA.put("GUA", "V");
        dnaToAA.put("GUC", "V");
        dnaToAA.put("GUG", "V");
        dnaToAA.put("GUU", "V");
        dnaToAA.put("UAA", " ");
        dnaToAA.put("UAC", "Y");
        dnaToAA.put("UAG", " ");
        dnaToAA.put("UAU", "Y");
        dnaToAA.put("UCA", "S");
        dnaToAA.put("UCC", "S");
        dnaToAA.put("UCG", "S");
        dnaToAA.put("UCU", "S");
        dnaToAA.put("UGA", " ");
        dnaToAA.put("UGC", "C");
        dnaToAA.put("UGG", "W");
        dnaToAA.put("UGU", "C");
        dnaToAA.put("UUA", "L");
        dnaToAA.put("UUC", "F");
        dnaToAA.put("UUG", "L");
        dnaToAA.put("UUU", "F");
    }

    public static void main(String[] args)
    {
        theSpectralConvolutionSavesTheDayStep4("188 1008 950 1065 588 877 319 507 218 978 87 1081 608 202 521 0 675 764 232 289 101 689 358 246 101 925 838 1040 822 156 776 188 374 115 606 721 806 707 822 420 287 1109 533 289 432 1008 1196 907 131 590 620 402 1010 333 907 707 317 776 576 794 1065 778 420 953 489 1095 879 620 475 475 489 243 1095 964 721 271 576 418 131 101 87 691 1109 186 505 374 390 1095 663 863 909 994");
    }

    //The Spectral Convolution Saves the Day | Step 4
    private static void theSpectralConvolutionSavesTheDayStep4(String spectrum)
    {
        ArrayList<Integer> convolution = spectralConvolutionProblem(spectrum);
        for (Integer integer : convolution)
        {
            System.out.print(integer + " ");
        }
    }

    //From 20 to More than 100 Amino Acids | Step 1
    private static void from20ToMoreThan100AminoAcids(String spectrum, int N)
    {
        //the idea is thes same as in Mass Spectrometry Meets Golf | Step 9
        //modifying aaMasses to include all the aminoacids with weights from 57 to 200
        aaMasses.clear();
        for (int i = 57; i < 201; i++)
        {
            char c = (char) i;
            aaMasses.put(String.valueOf(c), i);
        }
        Set<String> peptides = leaderboardCyclopeptideSequencingSet(spectrum, N);
        for (String peptide : peptides)
        {
            System.out.println(peptide);
        }
    }

    //Mass Spectrometry Meets Golf | Step 9
    private static void massSpectrometryMeetsGolfStep9(String spectrum, int N)
    {
        //to make it work we should modify aaMasses first.
        //by removing duplicates
        //in order to expand by mass rather than by letter
        // that approach allows to store more peptides in the leaderboard
        aaMasses.remove("L");
        aaMasses.remove("Q");
        Set<String> peptides = leaderboardCyclopeptideSequencingSet(spectrum, N);
        for (String peptide : peptides)
        {
            System.out.println(peptide);
        }
    }

    public static ArrayList<Integer> spectralConvolutionProblem(String spec)
    {
        ArrayList<Integer> result = new ArrayList<>();
        int[] spectrum = parseSpectrum(spec);
        Arrays.sort(spectrum);
        for (int i = 1; i < spectrum.length; i++)
        {
            int firstElement = spectrum[i];
            for (int j = i - 1; j >= 0; j--)
            {
                int secondElement = spectrum[j];
                if (firstElement - secondElement != 0)
                {
                    result.add(firstElement - secondElement);
                }
            }
        }
        return result;
    }

    public static Set<String> parsePeptides(String peptides)
    {
        Set<String> result = new HashSet<String>();
        result.addAll(Arrays.asList(peptides.split(" ")));
        return result;
    }

    public static void trim(Set<String> leaderBoard, String spectrum, int N/*, boolean last*/)
    {
        if (leaderBoard.isEmpty() || leaderBoard.size() < N)
        {
            return;
        }
        ArrayList<String> toRemove = new ArrayList<String>();
        ArrayList<Integer> scores = new ArrayList<Integer>();
        Map<String, Integer> peptideScore = new HashMap<String, Integer>();
        for (String peptide : leaderBoard)
        {
            int score;
            score = linearScore(peptide, spectrum);
            peptideScore.put(peptide, score);
            scores.add(score);
        }

        Collections.sort(scores);
        Integer max = 0;
        for (int i = 0; i < N; i++)
        {
            max = Collections.max(scores);
            scores.remove(max);
        }
        for (Map.Entry<String, Integer> stringIntegerEntry : peptideScore.entrySet())
        {
            if (stringIntegerEntry.getValue() < max)
            {
                toRemove.add(stringIntegerEntry.getKey());
            }
        }
        leaderBoard.removeAll(toRemove);
    }

    public static String leaderboardCyclopeptideSequencing(String s, int N)
    {
        Set<String> leaderboard = new HashSet<String>();
        leaderboard.add("");
        String leaderPeptide = "";
        int[] spectrum = parseSpectrum(s);
        int parentMass = parentMass(spectrum);
        while (!leaderboard.isEmpty())
        {
            expand(leaderboard);
            Set<String> toRemove = new HashSet<String>();
            for (String peptide : leaderboard)
            {
                if (massOfFragment(peptide) == parentMass)
                {
                    if (score(peptide, s) > score(leaderPeptide, s))
                    {
                        leaderPeptide = peptide;
                    }
                } else if (massOfFragment(peptide) > parentMass)
                {
                    toRemove.add(peptide);
                }
            }
            leaderboard.removeAll(toRemove);
            trim(leaderboard, s, N);
        }
        return outputPeptide(leaderPeptide);
    }

    public static Set<String> leaderboardCyclopeptideSequencingSet(String s, int N)
    {
        Set<String> result = new HashSet<>();
        Set<String> leaderboard = new HashSet<String>();
        HashMap<String, Integer> leadersScore = new HashMap<>();
        leaderboard.add("");
        String leaderPeptide = "";
        int[] spectrum = parseSpectrum(s);
        int parentMass = parentMass(spectrum);
        while (!leaderboard.isEmpty())
        {
            expand(leaderboard);
            Set<String> toRemove = new HashSet<>();
            for (String peptide : leaderboard)
            {
                if (massOfFragment(peptide) == parentMass)
                {
                    if (score(peptide, s) >= score(leaderPeptide, s))
                    {
                        leaderPeptide = peptide;
                        leadersScore.put(peptide, score(peptide, s));
                    }
                } else if (massOfFragment(peptide) > parentMass)
                {
                    toRemove.add(peptide);
                }
            }
            leaderboard.removeAll(toRemove);
            trim(leaderboard, s, N);
        }

        int score = score(leaderPeptide, s);
        for (Map.Entry<String, Integer> stringIntegerEntry : leadersScore.entrySet())
        {
            if (stringIntegerEntry.getValue().equals(score))
            {
                result.add(outputPeptide(stringIntegerEntry.getKey()));
            }
        }
        return result;
    }

    public static int score(String peptide, String spectrumString)
    {
        int score = 0;
        int[] spectrum = parseSpectrum(spectrumString);
        ArrayList<Integer> spec = listFromArray(spectrum);
        ArrayList<Integer> linearSpectrum = cyclospectrum(peptide);
        for (Integer integer : linearSpectrum)
        {
            if (spec.remove(integer))
            {
                score++;
            }
        }
        return score;
    }

    public static int linearScore(String peptide, String spectrumString)
    {
        int score = 0;
        int[] spectrum = parseSpectrum(spectrumString);
        ArrayList<Integer> spec = listFromArray(spectrum);
        ArrayList<Integer> linearSpectrum = linearSpectrum(peptide);
        for (Integer integer : linearSpectrum)
        {
            if (spec.remove(integer))
            {
                score++;
            }
        }
        return score;
    }

    public static Set<String> cyclopeptideSequencing(String s)
    {
        Set<String> result = new HashSet<String>();
        int[] spectrum = parseSpectrum(s);
        Set<String> peptides = new HashSet<String>();
        int parentMass = parentMass(spectrum);
        peptides.add("");
        while (!peptides.isEmpty())
        {
            expand(peptides);
            Set<String> peptidesToRemove = new HashSet<String>();
            for (String peptide : peptides)
            {
                if (massOfFragment(peptide) == parentMass)
                {
                    int[] cyclospecInts = getCyclospecInts(peptide);
                    Arrays.sort(cyclospecInts);
                    if (Arrays.equals(cyclospecInts, spectrum))
                    {
                        result.add(outputPeptide(peptide));
                    } else
                    {
                        peptidesToRemove.add(peptide);
                    }
                } else
                {
                    if (!isConsistentWithSpectrum(peptide, spectrum))
                    {
                        peptidesToRemove.add(peptide);
                    }
                }
            }
            peptides.removeAll(peptidesToRemove);
        }
        return result;
    }

    public static boolean isConsistentWithSpectrum(String peptide, int[] spectrum)
    {
        ArrayList<Integer> c = linearSpectrum(peptide);
        ArrayList<Integer> s = listFromArray(spectrum);
        for (Integer integer : c)
        {
            if (!s.remove(integer))
            {
                return false;
            }
        }
        return true;
    }

    private static ArrayList<Integer> listFromArray(int[] spectrum)
    {
        ArrayList<Integer> s = new ArrayList<Integer>(spectrum.length);
        for (int aSpectrum : spectrum)
        {
            s.add(aSpectrum);
        }
        return s;
    }

    private static int[] getCyclospecInts(String peptide)
    {
        ArrayList<Integer> cyclospec = cyclospectrum(peptide);
        int[] a = new int[cyclospec.size()];
        for (int i = 0; i < a.length; i++)
        {
            a[i] = cyclospec.get(i);
        }
        return a;
    }

    public static String outputPeptide(String peptide)
    {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < peptide.length(); i++)
        {
            char aa = peptide.charAt(i);
            sb.append(aaMasses.get(aa + "")).append("-");
        }
        sb.deleteCharAt(sb.lastIndexOf("-"));
        return sb.toString();
    }

    public static int parentMass(int[] spectrum)
    {
        //        Arrays.sort(spectrum);
        return spectrum[spectrum.length - 1];
    }

    public static void expand(Set<String> peptides)
    {
        Set<String> expanded = new HashSet<String>();
        Set<String> remove = new HashSet<String>();
        for (String peptide : peptides)
        {
            for (String s : aaMasses.keySet())
            {
                expanded.add(peptide + s);
                remove.add(peptide);
            }
        }
        peptides.addAll(expanded);
        peptides.removeAll(remove);
    }

    public static int[] parseSpectrum(String spectrum)
    {
        String[] nums = spectrum.split(" ");
        int[] result = new int[nums.length];
        for (int i = 0; i < nums.length; i++)
        {
            result[i] = Integer.valueOf(nums[i]);
        }
        return result;
    }


    public static int countLinearPeptides(int n)
    {
        return ((1 + n) * n / 2) + 1;
    }

    public static ArrayList<Integer> cyclospectrum(String peptide)
    {
        ArrayList<Integer> result = new ArrayList<Integer>();
        String p = peptide + peptide;
        for (int j = peptide.length() - 1; j > 0; j--)
        {
            for (int i = 0; i < peptide.length(); i++)
            {
                String subfragment = p.substring(i, i + j);
                result.add(massOfFragment(subfragment));
            }
        }
        result.add(0);
        result.add(massOfFragment(peptide));
        return result;
    }

    public static ArrayList<Integer> linearSpectrum(String peptide)
    {
        ArrayList<Integer> result = new ArrayList<Integer>();
        for (int j = peptide.length() - 1; j > 0; j--)
        {
            for (int i = 0; i < peptide.length() - j + 1; i++)
            {
                String subfragment = peptide.substring(i, i + j);
                result.add(massOfFragment(subfragment));
            }
        }
        result.add(0);
        result.add(massOfFragment(peptide));
        return result;
    }

    public static int massOfFragment(String fragment)
    {
        int result = 0;
        char[] aas = fragment.toCharArray();
        for (char aa : aas)
        {
            result += aaMasses.get(aa + "");
        }
        return result;
    }

    public static String translate(String rna)
    {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < rna.length(); i += 3)
        {
            String triplet = rna.substring(i, i + 3);
            String aminoAcid = dnaToAA.get(triplet);
            if (!aminoAcid.equals(" "))
            {
                sb.append(aminoAcid);
            } else
            {
                return sb.toString();
            }
        }
        return sb.toString();
    }

    public static ArrayList<String> peptideEncoding(String genome, String peptide)
    {
        ArrayList<String> result = new ArrayList<String>();
        Set<String> allRNASeq = getAllRNASeqForPeptide(peptide);
        for (String rna : allRNASeq)
        {
            String dna = rnaToDNA(rna);
            String reverseDNA = PatternCount.reverseCompliment(dna);
            Pattern p = Pattern.compile(dna);
            Matcher m = p.matcher(genome);
            while (m.find())
            {
                result.add(dna);
            }
            p = Pattern.compile(reverseDNA);
            m = p.matcher(genome);
            while (m.find())
            {
                result.add(reverseDNA);
            }
        }
        return result;
    }


    public static void addNextAAToRNA(Set<String> rnas, char aa)
    {
        Set<String> allRNAForAA = getAllRNAForAA(aa);
        Set<String> rnasToDelete = new HashSet<String>(rnas);
        if (rnas.size() == 0)
        {
            for (String s : allRNAForAA)
            {
                rnas.add(s);
            }

        } else
        {
            for (String s : allRNAForAA)
            {
                for (String rna : rnasToDelete)
                {
                    rnas.add(rna + s);
                }
            }
            rnas.removeAll(rnasToDelete);
        }
    }

    public static String rnaToDNA(String rna)
    {
        return rna.replaceAll("U", "T");
    }

    public static Set<String> getAllRNASeqForPeptide(String peptide)
    {
        Set<String> result = new HashSet<String>();
        char[] aas = peptide.toCharArray();
        for (char aa : aas)
        {
            addNextAAToRNA(result, aa);
        }
        return result;
    }

    public static Set<String> getAllRNAForAA(char aa)
    {
        Set<String> result = new HashSet<String>();
        for (Map.Entry<String, String> stringStringEntry : dnaToAA.entrySet())
        {
            if (stringStringEntry.getValue().charAt(0) == aa)
            {
                result.add(stringStringEntry.getKey());
            }
        }

        return result;
    }
}
