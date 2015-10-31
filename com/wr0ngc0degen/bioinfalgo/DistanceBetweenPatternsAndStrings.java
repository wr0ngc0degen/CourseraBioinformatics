package com.wr0ngc0degen.bioinfalgo;

/**
 * Created by Alena on 10/31/2015.
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.util.stream.Collectors;

class DistanceBetweenPatternsAndStrings
{
    public static void main(String[] args)
    {
        Scanner scanner = new Scanner(System.in);
        String pattern = scanner.nextLine();
        String dnasString = scanner.nextLine();
        List<String> dnas = Arrays.stream(dnasString.split(" ")).collect(Collectors.toList());
        System.out.println(dnas.stream().mapToInt(value -> minDistancePatternFromDNA(pattern, value)).sum());
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
        return dnaKmers.stream().mapToInt(value -> hammingDistance(value, pattern)).min().getAsInt();
    }

    public static int hammingDistance(String s1, String s2)
    {
        int distance = 0;
        char[] charS1 = s1.toCharArray();
        for (int i = 0; i < charS1.length; i++)
        {
            char c = charS1[i];
            if (c != s2.charAt(i))
            {
                distance++;
            }
        }
        return distance;
    }
}