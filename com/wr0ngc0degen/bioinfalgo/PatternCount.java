package com.wr0ngc0degen.bioinfalgo;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Created by Alena on 31.10.2014.
 */
public class PatternCount
{
    public static void main(String[] args)
    {
/*
        Scanner scanner = readFromFile("E-coli.txt");
        String genome = "";
        if (scanner.hasNextLine())
        {
            genome = scanner.nextLine();
        }
        Set<String> fw = betterClumpFinding(genome, 9, 500, 3);
        for (String s : fw)
        {
            System.out.print(s + " ");
        }
        System.out.println(fw.size());
*/

        System.out.println(patternToNumber2("CTTCTCACGTACAACAAAATC"));
        System.out.println(patternToNumber("CTTCTCACGTACAACAAAATC"));

        ///////////////////////////
/*
        Set<String> result = new HashSet<String>();
        Set<String> kmers = new HashSet<String>();
        for (int i = 0; i < (int) Math.pow(4, 4); i++)
        {
            kmers.add(numberToPattern(i, 4));
        }
        for (String kmer : kmers)
        {
            if(hammingDistance("TGCA", kmer) <= 3)
            {
                result.add(kmer);
            }
        }
        System.out.println(result.size());

        System.out.println("reverseCompliment");
        System.out.println(reverseCompliment("GCTAGCT"));

        System.out.println("hammingDistance");
        System.out.println(hammingDistance("CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA", "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG"));

        System.out.println("minSkew");
        int[] computeSkew = computeSkew("CATTCCAGTACTTCATGATGGCGTGAAGA");
        Set<Integer> minSkew = minSkew(computeSkew);
        System.out.println(minSkew);
        for (int i = 0; i < computeSkew.length; i++)
        {
            System.out.print(computeSkew[i]);

        }
        System.out.println();

        System.out.println("count");
        System.out.println(count("TACGCATTACAAAGCACA", "AA", 1));
*/
    }

    private static void quizOne()
    {
        //pattern count
        System.out.println("patternCount = \n" + patternCount("ACTGTACGATGATGTGTGTCAAAG", "TGT"));
        //frequent words
        Set<String> frequentWords = frequentWords("CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA", 3);
        System.out.println("frequentWords");
        for (String frequentWord : frequentWords)
        {
            System.out.print(frequentWord + " ");
        }
        //reverse complement
        System.out.println("\nreverse complement");
        System.out.println(reverseCompliment("CCAGATC"));
    }

    public static Set<String> frequentWordsWithMismatchAndReverse(String text, int k, int dist)
    {
        Set<String> result = new HashSet<String>();
        Set<String> kmers = new HashSet<String>();
        for (int i = 0; i < (int) Math.pow(4, k); i++)
        {
            kmers.add(numberToPattern(i, k));
        }

        for (int i = 0; i < text.length() - k + 1; i++)
        {
            String kmer = text.substring(i, i + k);
            kmers.add(kmer);
        }
        HashMap<String, Integer> kmersFreq = new HashMap<String, Integer>();
        for (String kmer : kmers)
        {
            for (int i = 0; i < text.length() - k + 1; i++)
            {
                String sub = text.substring(i, i + k);
                if (hammingDistance(kmer, sub) <= dist)
                {
                    if (kmersFreq.get(kmer) != null)
                    {
                        int f = kmersFreq.get(kmer);
                        kmersFreq.put(kmer, f + 1);
                    } else
                    {
                        kmersFreq.put(kmer, 1);
                    }
                }
                if (hammingDistance(reverseCompliment(kmer), sub) <= dist)
                {
                    if (kmersFreq.get(kmer) != null)
                    {
                        int f = kmersFreq.get(kmer);
                        kmersFreq.put(kmer, f + 1);
                    } else
                    {
                        kmersFreq.put(kmer, 1);
                    }
                }
            }
        }
        int max = Collections.max(kmersFreq.values());
        for (Map.Entry<String, Integer> stringIntegerEntry : kmersFreq.entrySet())
        {
            if (stringIntegerEntry.getValue().equals(max))
            {
                result.add(stringIntegerEntry.getKey());
            }
        }
        return result;
    }

    public static Set<String> frequentWordsWithMismatch(String text, int k, int dist)
    {
        Set<String> result = new HashSet<String>();
        Set<String> kmers = new HashSet<String>();
        for (int i = 0; i < (int) Math.pow(4, k); i++)
        {
            kmers.add(numberToPattern(i, k));
        }

        for (int i = 0; i < text.length() - k + 1; i++)
        {
            String kmer = text.substring(i, i + k);
            kmers.add(kmer);
        }
        HashMap<String, Integer> kmersFreq = new HashMap<String, Integer>();
        for (String kmer : kmers)
        {
            for (int i = 0; i < text.length() - k + 1; i++)
            {
                String sub = text.substring(i, i + k);
                if (hammingDistance(kmer, sub) <= dist)
                {
                    if (kmersFreq.get(kmer) != null)
                    {
                        int f = kmersFreq.get(kmer);
                        kmersFreq.put(kmer, f + 1);
                    } else
                    {
                        kmersFreq.put(kmer, 1);
                    }
                }
            }
        }
        int max = Collections.max(kmersFreq.values());
        for (Map.Entry<String, Integer> stringIntegerEntry : kmersFreq.entrySet())
        {
            if (stringIntegerEntry.getValue().equals(max))
            {
                result.add(stringIntegerEntry.getKey());
            }
        }
        return result;
    }

    public static int count(String genome, String pattern, int d)
    {
        return approximatePatternMatch(genome, pattern, d).size();
    }

    public static List<Integer> approximatePatternMatch(String text, String pattern, int maxDist)
    {
        List<Integer> result = new ArrayList<Integer>();
        for (int i = 0; i < text.length() - pattern.length() + 1; i++)
        {
            String piece = text.substring(i, i + pattern.length());
            if (hammingDistance(pattern, piece) <= maxDist)
            {
                result.add(i);
            }
        }
        return result;
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

    public static Set<Integer> minSkew(int[] skews)
    {
        Set<Integer> result = new HashSet<Integer>();
        int[] skewsCopy = Arrays.copyOf(skews, skews.length);
        Arrays.sort(skewsCopy);
        int min = skewsCopy[0];
        for (int i = 0; i < skews.length; i++)
        {
            int skew = skews[i];
            if (skew == min)
            {
                result.add(i);
            }
        }
        return result;
    }

    public static int[] computeSkew(String genome)
    {
        int[] skews = new int[genome.length() + 1];
        char[] genomeChars = genome.toCharArray();
        skews[0] = 0;
        for (int i = 1; i < skews.length; i++)
        {
            if (genomeChars[i - 1] == 'C')
            {
                skews[i] = skews[i - 1] - 1;
            } else if (genomeChars[i - 1] == 'G')
            {
                skews[i] = skews[i - 1] + 1;
            } else
            {
                skews[i] = skews[i - 1];
            }

        }
        return skews;
    }

    public static Set<String> fasterFrequentWords(String text, int k)
    {
        Set<String> frequentPatterns = new HashSet<String>();
        int[] freqArray = frequencyArray(text, k);
        int max = 0;
        for (int i : freqArray)
        {
            if (i > max)
            {
                max = i;
            }
        }
        for (int i = 0; i < freqArray.length; i++)
        {
            if (freqArray[i] == max)
            {
                String pattern = numberToPattern(i, k);
                frequentPatterns.add(pattern);
            }

        }
        return frequentPatterns;
    }


    public static Set<String> findFrequentWordsBySorting(String text, int k)
    {
        Set<String> frequentPatterns = new HashSet<String>();
        int[] index = new int[text.length() - k + 1];
        int[] count = new int[index.length];
        for (int i = 0; i < text.length() - k + 1; i++)
        {
            String pattern = text.substring(i, i + k);
            index[i] = patternToNumber(pattern);
            count[i] = 1;
        }

        Arrays.sort(index);
        for (int i = 1; i < count.length; i++)
        {
            if (index[i] == index[i - 1])
            {
                count[i] = count[i - 1] + 1;
            }
        }
        int max = 0;
        for (int i : count)
        {
            if (i > max)
            {
                max = i;
            }
        }
        for (int i = 0; i < count.length; i++)
        {
            if (count[i] == max)
            {
                String pattern = numberToPattern(index[i], k);
                frequentPatterns.add(pattern);
            }
        }
        return frequentPatterns;
    }

    public static int[] frequencyArray(String text, int k)
    {
        int[] freqArray = new int[(int) Math.pow(4, k)];
        for (int i = 0; i < freqArray.length; i++)
        {
            freqArray[i] = 0;
        }
        for (int i = 0; i < text.length() - k + 1; i++)
        {
            String pattern = text.substring(i, i + k);
            int j = patternToNumber(pattern);
            freqArray[j]++;
        }
        return freqArray;
    }

    public static int patternToNumber(String pattern)
    {
        char[] patternChars = pattern.toCharArray();
        int result = 0;
        for (int i = patternChars.length - 1, j = 0; i >= 0; i--, j++)
        {
            char patternChar = patternChars[i];
            result += Math.pow(4, j) * charToInt(patternChar);
        }

        return result;
    }

    public static int patternToNumber2(String pattern)
    {
        int length = pattern.length();
        if (length == 0)
        {
            return 0;
        } else
        {
            char symbol = pattern.charAt(length - 1);
            pattern = pattern.substring(0, length - 1);
            return 4 * patternToNumber2(pattern) + charToInt(symbol);
        }
    }

    public static int charToInt(char c)
    {
        switch (c)
        {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
            default:
                return 0;
        }
    }

    public static char intToChar(int i)
    {
        switch (i)
        {
            case 0:
                return 'A';
            case 1:
                return 'C';
            case 2:
                return 'G';
            case 3:
                return 'T';
            default:
                return 'A';
        }
    }

    public static String numberToPattern(int number, int size)
    {
        StringBuffer sb = new StringBuffer();
        for (int i = size - 1; i >= 0; i--)
        {
            int letter = (int) (number / Math.pow(4, i));
            sb.append(intToChar(letter));
            number = (int) (number % Math.pow(4, i));
        }
        return sb.toString();
    }

    public static int patternCount(String text, String pattern)
    {
        int count = 0;
        for (int i = 0; i < text.length() - pattern.length() + 1; i++)
        {
            if (text.substring(i, i + pattern.length()).equals(pattern))
            {
                count++;
            }

        }
        return count;
    }

    public static Scanner readFromFile(String fileName)
    {
        Scanner scanner = null;
        try
        {

            File file = new File(fileName);
            if (file.exists())
            {
                scanner = new Scanner(file, "UTF-8");
            }
        } catch (IOException ioe)
        {
            System.err.println("Could not open " + fileName);
        }

        return scanner;
    }

    public static Set<String> frequentWords(String text, int k)
    {
        Set<String> frequentPatterns = new HashSet<String>();
        ArrayList<Integer> count = new ArrayList<Integer>(text.length() - k + 1);
        String pattern;
        for (int i = 0; i < text.length() - k; i++)
        {
            pattern = text.substring(i, i + k);
            count.add(i, patternCount(text, pattern));
        }
        Integer max = Collections.max(count);
        for (int i = 0; i < text.length() - k; i++)
        {
            if (count.get(i).equals(max))
            {
                frequentPatterns.add(text.substring(i, i + k));
            }
        }
        return frequentPatterns;
    }

    public static String reverseCompliment(String text)
    {
        StringBuilder sb = new StringBuilder();
        for (int i = text.toCharArray().length - 1; i >= 0; i--)
        {
            char c = text.toCharArray()[i];
            if (c == 'A')
            {
                sb.append('T');
            }
            if (c == 'T')
            {
                sb.append('A');
            }
            if (c == 'G')
            {
                sb.append('C');
            }
            if (c == 'C')
            {
                sb.append('G');
            }
        }
        return sb.toString();
    }

    public static List<Integer> occurrencesOfPattern(String pattern, String text)
    {
        List<Integer> occurrences = new ArrayList<Integer>();
        for (int i = 0; i < text.length() - pattern.length() + 1; i++)
        {
            if (text.substring(i, i + pattern.length()).equals(pattern))
            {
                occurrences.add(i);
            }
        }
        return occurrences;
    }

    public static Set<String> findClumps(String genome, int k, int L, int t)
    {
        Set<String> result = new HashSet<String>();
        for (int i = 0; i < genome.length() - L + 1; i++)
        {
            String piece = genome.substring(i, L + i);
            Set<String> allPatterns = findFrequentWordsBySorting(piece, k);
            for (String pattern : allPatterns)
            {
                if (patternCount(piece, pattern) >= t)
                {
                    result.add(pattern);
                }
            }
        }
        return result;
    }

    public static Set<String> betterClumpFinding(String genome, int k, int L, int t)
    {
        Set<String> frequentPatterns = new HashSet<String>();
        int[] clumps = new int[(int) Math.pow(4, k)];
        for (int i = 0; i < clumps.length; i++)
        {
            clumps[i] = 0;
        }
        String text = genome.substring(0, L);
        int[] freqArray = frequencyArray(text, k);
        for (int i = 0; i < freqArray.length; i++)
        {
            if (freqArray[i] >= t)
            {
                clumps[i] = 1;
            }
        }
        for (int i = 1; i < genome.length() - L; i++)
        {
            String firstPattern = genome.substring(i - 1, i - 1 + k);
            int j = patternToNumber(firstPattern);
            freqArray[j]--;
            String lastPattern = genome.substring(i + L - k, i + L);
            j = patternToNumber(lastPattern);
            freqArray[j]++;
            if (freqArray[j] >= t)
            {
                clumps[j] = 1;
            }
        }
        for (int i = 0; i < clumps.length; i++)
        {
            if (clumps[i] == 1)
            {
                String pattern = numberToPattern(i, k);
                frequentPatterns.add(pattern);
            }

        }
        return frequentPatterns;
    }
}
