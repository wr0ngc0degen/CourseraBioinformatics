package com.wr0ngc0degen.bioinfalgo;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created by Alena on 30.11.2014.
 */
public class MotifFinding
{
    public static char[] alphabet = new char[]{'A', 'C', 'G', 'T'};

    public static void main(String[] args) throws FileNotFoundException
    {
        motifEnumeration("dataset_156_7.txt");
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
        Set<String> patterns = new HashSet<String>();
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
        Set<String> allKMersFromDNAs = new HashSet<String>();
        for (String dna : dnas)
        {
            for (int i = 0; i < dna.length() - k; i++)
            {
                String kmer = dna.substring(i, i + k);
                allKMersFromDNAs.add(kmer);
            }
        }
        return allKMersFromDNAs;
    }

    public static Set<String> immediateNeighbors(String pattern)
    {
        Set<String> neighborhood = new HashSet<String>();
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
        Set<String> neighborhood = new HashSet<String>();
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
