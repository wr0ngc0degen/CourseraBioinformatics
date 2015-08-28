import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.LinkedList;
import java.util.Scanner;

/**
 * Created by Alena on 26.08.2015.
 */
public class FragileRegions
{
    static StringBuilder builder = new StringBuilder();
    public static void main(String[] args) throws Exception
    {
        greedySorting("dataset_286_3.txt");
        PrintWriter out = new PrintWriter("result.txt");
        out.print(builder);
        out.flush();
    }

    //A Greedy Algorithm for Sorting by Reversals | Step 3
    private static void greedySorting(String fileName) throws FileNotFoundException
    {
        Scanner scanner = new Scanner(new File(fileName));
        String input = scanner.nextLine();
        LinkedList<Integer> permutation = convertStringToPermutation(input);
        greedySorting(permutation);
    }

    private static int greedySorting(LinkedList<Integer> permutation)
    {
        int approxReversalDistance = 0;
        for (int k = 0; k < permutation.size(); k++)
        {
            //element k is not sorted
            if (Math.abs(permutation.get(k)) != k + 1)
            {
                int endPos = permutation.indexOf(k + 1) == -1 ? permutation.indexOf(-(k + 1)) : permutation.indexOf(k + 1);
                kSortingReversal(permutation, k, endPos);
                builder.append(convertPermutationToString(permutation)).append("\n");
                approxReversalDistance++;
            }
            if (permutation.get(k) == -(k + 1))
            {
                permutation.set(k, k + 1);
                builder.append(convertPermutationToString(permutation)).append("\n");
                approxReversalDistance++;
            }
        }
        return approxReversalDistance;
    }

    private static void kSortingReversal(LinkedList<Integer> permutation, int startPos, int endPos)
    {
        for (int i = startPos, j = endPos; i <= Math.floorDiv((endPos - startPos) , 2) + startPos; i++, j--)
        {
            exchange(permutation, i, j);
        }
    }

    private static void exchange(LinkedList<Integer> permutation, int posFirst, int posLast)
    {
        int first = permutation.get(posFirst);
        permutation.set(posFirst, -permutation.get(posLast));
        permutation.set(posLast, -first);
    }

    private static LinkedList<Integer> convertStringToPermutation(String input)
    {
        LinkedList<Integer> permutation = new LinkedList<>();
        String[] elements = input.replace("(", "").replace(")", "").split(" ");
        for (String element : elements)
        {
            permutation.add(Integer.parseInt(element));
        }
        return permutation;
    }

    private static String convertPermutationToString(LinkedList<Integer> permutation)
    {
        StringBuilder sb = new StringBuilder("(");
        for (Integer integer : permutation)
        {
            if (integer > 0)
            {
                sb.append("+");
            }
            sb.append(integer).append(" ");
        }
        return sb.deleteCharAt(sb.length() - 1).append(")").toString();
    }
}
