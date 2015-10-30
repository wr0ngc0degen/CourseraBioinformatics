package com.wr0ngc0degen.bioinfalgo;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by Alena on 30.11.2014.
 */
public class MotifFinding
{
    public static char[] alphabet = new char[]{'A', 'C', 'G', 'T'};

    public static void main(String[] args) throws FileNotFoundException
    {
        //        motifEnumeration("dataset_156_7.txt");
        //        medianString("dataset_158_9.txt");
        //        profileMostProbableKmer("dataset_159_3.txt");
        //        greedyMotifSearch("dataset_159_5.txt");
    }

    /*
    CODE CHALLENGE: Implement GREEDYMOTIFSEARCH.
    Input: Integers k and t, followed by a collection of strings Dna.
    Output: A collection of strings BestMotifs resulting from applying GREEDYMOTIFSEARCH(Dna,k,t).
    If at any step you find more than one Profile-most probable k-mer in a given string, use the
    one occurring first.
     */
    //Greedy Motif Search | Step 5
    private static void greedyMotifSearch(String fileName) throws FileNotFoundException
    {
        Scanner scanner = new Scanner(new File(fileName));
        String[] kt = scanner.nextLine().split(" ");
        int k = Integer.parseInt(kt[0]);
        int t = Integer.parseInt(kt[1]);
        List<String> dnas = new ArrayList<>();
        while (scanner.hasNextLine())
        {
            dnas.add(scanner.nextLine());
        }
        List<String> bestMotifs = greedyMotifSearch(dnas, k, t);
        bestMotifs.stream().forEach(System.out::println);
    }

    private static List<String> greedyMotifSearch(List<String> dnas, int k, int numberOfSequences)
    {
        List<String> bestMotifs = dnas.stream().map(s -> s.substring(0, k)).collect(Collectors.toList());
        String firstDNA = dnas.get(0);

        for (int i = 0; i < firstDNA.length() - k + 1; i++)
        {
            List<String> motifs = new ArrayList<>(numberOfSequences);
            String kmer = firstDNA.substring(i, i + k);
            motifs.add(kmer);
            for (int j = 1; j < numberOfSequences; j++)
            {
                String currentDNA = dnas.get(j);
                double[][] profile = calculateProfile(motifs);
                String kmerToAdd = profileMostProbableKmer(currentDNA, k, profile);
                motifs.add(kmerToAdd);
            }
            if (distanceFromPatternToDNAs(findConsensus(bestMotifs), bestMotifs) > distanceFromPatternToDNAs(findConsensus(motifs), motifs))
            {
                bestMotifs = motifs;
            }
        }
        return bestMotifs;
    }

    private static String findConsensus(List<String> motifs)
    {
        double[][] profile = calculateProfile(motifs);
        char[] consensus = new char[profile[0].length];

        //transpose the profile to ease calculation
        double[][] transposed = new double[profile[0].length][profile.length];
        for (int i = 0; i < profile.length; i++)
        {
            double[] doubles = profile[i];
            for (int j = 0; j < doubles.length; j++)
            {
                transposed[j][i] = profile[i][j];
            }
        }
        for (int i = 0; i < transposed.length; i++)
        {
            double[] doubles = transposed[i];
            double max = Arrays.stream(doubles).max().getAsDouble();
            for (int j = 0; j < doubles.length; j++)
            {
                double aDouble = doubles[j];
                if (aDouble == max)
                {
                    consensus[i] = PatternCount.intToChar(j);
                }
            }
        }
        return new String(consensus);
    }

    private static double[][] calculateProfile(List<String> motifs)
    {
        String firstMotif = motifs.get(0);
        double[][] profile = new double[4][firstMotif.length()];

        int numberOfSequences = motifs.size();

        //filling the matrix with count of every letter
        for (int i = 0; i < firstMotif.length(); i++)
        {
            for (String motif : motifs)
            {
                char charAtI = motif.charAt(i);
                profile[PatternCount.charToInt(charAtI)][i] += 1;
            }
        }
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < firstMotif.length(); j++)
            {
                profile[i][j] = profile[i][j] / numberOfSequences;
            }
        }
        return profile;
    }

    /*
    CODE CHALLENGE: Solve the Profile-most Probable k-mer Problem.
    Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.
     Input: A string Text, an integer k, and a 4 ? k matrix Profile.
     Output: A Profile-most probable k-mer in Text.
     */
    //Greedy Motif Search | Step 3
    private static void profileMostProbableKmer(String fileName) throws FileNotFoundException
    {
        Scanner scanner = new Scanner(new File(fileName));
        String text = scanner.nextLine();
        int k = Integer.parseInt(scanner.nextLine());
        double[][] matrix = new double[4][text.length()];
        for (int i = 0; i < 4; i++)
        {
            String line = scanner.nextLine();
            matrix[i] = Arrays.stream(line.split(" ")).mapToDouble(Double::parseDouble).toArray();
        }
        System.out.println(profileMostProbableKmer(text, k, matrix));
    }

    private static String profileMostProbableKmer(String text, int k, double[][] matrix)
    {
        List<String> kmers = new ArrayList<>();
        for (int i = 0; i < text.length() - k + 1; i++)
        {
            String kmer = text.substring(i, i + k);
            kmers.add(kmer);
        }
        return kmers.stream().max((o1, o2) -> ((Double) getScore(o1, matrix)).compareTo(getScore(o2, matrix))).get();
    }

    private static double getScore(String motif, double[][] matrix)
    {
        double result = 1;
        char[] chars = motif.toCharArray();
        for (int i = 0; i < chars.length; i++)
        {
            result *= getProbForChar(matrix, i, chars[i]);
        }
        return result;
    }

    private static double getProbForChar(double[][] matrix, int pos, char ch)
    {
        return matrix[PatternCount.charToInt(ch)][pos];
    }

    /*
    CODE CHALLENGE: Implement MEDIANSTRING.
     Input: An integer k, followed by a collection of strings Dna.
     Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all k-mers Pattern. (If there are
     multiple such strings Pattern, then you may return any one.)
     */
    //From Motif Finding to Finding a Median String | Step 9
    private static void medianString(String fileName) throws FileNotFoundException
    {
        Scanner scanner = new Scanner(new File(fileName));
        String line = scanner.nextLine();
        int k = Integer.parseInt(line);
        List<String> dnas = new ArrayList<>();
        while (scanner.hasNextLine())
        {
            dnas.add(scanner.nextLine());
        }
        System.out.println(medianString(k, dnas));
    }

    private static String medianString(int k, List<String> dnas)
    {
        Set<String> allKmers = PatternCount.generateAllKmers(k);
        int distance = Integer.MAX_VALUE;
        String median = "";
        for (String pattern : allKmers)
        {
            int distanceFromPatternToDNAs = distanceFromPatternToDNAs(pattern, dnas);
            if (distance > distanceFromPatternToDNAs)
            {
                distance = distanceFromPatternToDNAs;
                median = pattern;
            }
        }
        return median;
    }

    private static int distanceFromPatternToDNAs(String pattern, List<String> dnas)
    {
        return dnas.stream().mapToInt(value -> minDistancePatternFromDNA(pattern, value)).sum();
    }

    private static int minDistancePatternFromDNA(String pattern, String dna)
    {
        int k = pattern.length();
        int numKmersInDNA = dna.length() - k + 1;
        List<String> dnaKmers = new ArrayList<>(numKmersInDNA);
        for (int i = 0; i < numKmersInDNA; i++)
        {
            dnaKmers.add(dna.substring(i, i + k));
        }
        return dnaKmers.stream().mapToInt(value -> PatternCount.hammingDistance(value, pattern)).min().getAsInt();
    }

    /*
    CODE CHALLENGE: Implement MOTIFENUMERATION.
     Input: Integers k and d, followed by a collection of strings Dna.
     Output: All (k, d)-motifs in Dna
     */
    //Motif Finding Is More Difficult Than You Think | Step 7
    private static void motifEnumeration(String fileName) throws FileNotFoundException
    {
        Scanner scanner = new Scanner(new File(fileName));
        String line = scanner.nextLine();
        String[] kd = line.split(" ");
        int k = Integer.parseInt(kd[0]);
        int d = Integer.parseInt(kd[1]);
        Set<String> dnas = new HashSet<>();
        while (scanner.hasNextLine())
        {
            dnas.add(scanner.nextLine());
        }

        Set<String> patterns = motifEnumeration(dnas, k, d);
        List<String> patternList = new ArrayList<>(patterns);
        Collections.sort(patternList);
        for (String s : patternList)
        {
            System.out.print(s + " ");
        }
    }

    private static Set<String> motifEnumeration(Set<String> dnas, int k, int d)
    {
        Set<String> patterns = new HashSet<>();
        Set<String> allKMersFromDNAs = getAllKmers(dnas, k);

        for (String kmer : allKMersFromDNAs)
        {
            Set<String> kmerNeighbors = neighbors(kmer, d);
            for (String kmerNeighbor : kmerNeighbors)
            {
                Set<String> newNeighbors = neighbors(kmerNeighbor, d);
                if (checkAllDNAContainsKmerFamily(dnas, newNeighbors))
                {
                    patterns.add(kmerNeighbor);
                }
            }
        }
        return patterns;
    }

    private static boolean checkAllDNAContainsKmerFamily(Set<String> dnas, Set<String> neighbors)
    {
        for (String dna : dnas)
        {
            if (!checkIfDNAContainsKmerFamily(dna, neighbors))
            {
                return false;
            }
        }
        return true;
    }

    private static boolean checkIfDNAContainsKmerFamily(String dna, Set<String> neighbors)
    {
        for (String neighbor : neighbors)
        {
            if (dna.contains(neighbor))
            {
                return true;
            }
        }
        return false;
    }

    private static Set<String> getAllKmers(Set<String> dnas, int k)
    {
        Set<String> allKMersFromDNAs = new HashSet<>();
        for (String dna : dnas)
        {
            for (int i = 0; i < dna.length() - k + 1; i++)
            {
                String kmer = dna.substring(i, i + k);
                allKMersFromDNAs.add(kmer);
            }
        }
        return allKMersFromDNAs;
    }

    public static Set<String> immediateNeighbors(String pattern)
    {
        Set<String> neighborhood = new HashSet<>();
        neighborhood.add(pattern);
        int length = pattern.length();
        for (int i = 0; i < length; i++)
        {
            char symbol = pattern.charAt(i);

            for (char c : alphabet)
            {
                if (c != symbol)
                {
                    String neighbor = pattern.substring(0, i) + c + pattern.substring(i + 1, length);
                    neighborhood.add(neighbor);
                }
            }
        }
        return neighborhood;
    }

    public static Set<String> neighbors(String pattern, int d)
    {
        Set<String> neighborhood = new HashSet<>();
        if (d == 0)
        {
            neighborhood.add(pattern);
            return neighborhood;
        }
        if (pattern.length() == 1)
        {
            for (char c : alphabet)
            {
                neighborhood.add(c + "");
            }
            return neighborhood;
        }
        String suffix = pattern.substring(1, pattern.length());
        Set<String> suffixNeighbors = neighbors(suffix, d);
        for (String suffixNeighbor : suffixNeighbors)
        {
            if (PatternCount.hammingDistance(suffix, suffixNeighbor) < d)
            {
                for (char c : alphabet)
                {
                    neighborhood.add(c + suffixNeighbor);
                }
            } else
            {
                neighborhood.add(pattern.charAt(0) + suffixNeighbor);
            }
        }
        return neighborhood;
    }

}
