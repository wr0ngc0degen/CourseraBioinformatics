import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class SARS
{
    public static void main(String[] args) throws Exception
    {
        limbLengthProblem("dataset_10329_11.txt");
        //        distanceBetweenLeaves("dataset_10328_11.txt");
    }

    //Toward An Algorithm for Distance-Based Phylogeny Construction | Step 11
    /*
     CODE CHALLENGE: Solve the Limb Length Problem.
     Input: An integer n, followed by an integer j between 0 and n - 1, followed by a space-separated
     additive distance matrix D (whose elements are integers).
     Output: The limb length of the leaf in Tree(D) corresponding to row j of this distance
     matrix (use 0-based indexing).
     */
    private static void limbLengthProblem(String fileName) throws FileNotFoundException
    {
        Scanner scanner = new Scanner(new File(fileName));
        String line = scanner.nextLine();
        int n = Integer.parseInt(line);
        line = scanner.nextLine();
        int j = Integer.parseInt(line);
        int[][] matrix = new int[n][n];
        for (int i = 0; i < n; i++)
        {
            line = scanner.nextLine();
            matrix[i] = Arrays.stream(line.split(" ")).mapToInt(Integer::parseInt).toArray();
        }
        System.out.println(limbLenghtProblem(n, j, matrix));
    }

    private static int limbLenghtProblem(int n, int j, int[][] matrix)
    {
        List<Integer> distances = new ArrayList<>();
        for (int i = 0; i < matrix.length; i++)
        {
            if (i == j)
            {
                continue;
            }
            for (int k = i + 1; k < matrix.length; k++)
            {
                if(k == j)
                {
                    continue;
                }
                distances.add(limbLenght(i, j, k, matrix));
            }
        }
        return distances.stream().min(Integer::compare).get();
    }

    private static int limbLenght(int i, int j, int k, int[][] matrix)
    {
        return (matrix[i][j] + matrix[j][k] - matrix[i][k]) / 2;
    }

    //Transforming Distance Matrices into Evolutionary Trees | Step 11
    private static void distanceBetweenLeaves(String fileName) throws Exception
    {
        Scanner scanner = new Scanner(new File(fileName));
        String line = scanner.nextLine();
        Integer leavesNumber = Integer.parseInt(line);

        EdgeWeightedDigraph edgeWeightedDigraph = Utilities.constructEdgeWeightedDigraphFromString(fileName, 1);
        Digraph digraph = Utilities.constructDigraphFromString(fileName, 1);
        for (int leaveFrom = 0; leaveFrom < leavesNumber; leaveFrom++)
        {
            BreadthFirstDirectedPaths paths = new BreadthFirstDirectedPaths(digraph, leaveFrom);
            for (int leaveTo = 0; leaveTo < leavesNumber; leaveTo++)
            {
                Iterable<DirectedEdge> edges = edgeWeightedDigraph.edges();
                int distanceBetweenLeaves = 0;
                Integer vertexFrom = leaveFrom;
                for (Integer vertexTo : paths.pathTo(leaveTo))
                {
                    Stream<DirectedEdge> edgesStream = StreamSupport.stream(edges.spliterator(), false);
                    if (!vertexFrom.equals(vertexTo))
                    {
                        Integer fromInner = vertexFrom;
                        Integer toInner = vertexTo;

                        distanceBetweenLeaves += edgesStream.filter((s) -> s.from() == fromInner && s.to() == toInner).findFirst().get().weight();
                        vertexFrom = vertexTo;
                    }
                }
                System.out.print(distanceBetweenLeaves + " ");
            }
            System.out.println("");
        }

    }
}
