import java.io.File;
import java.util.Scanner;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class SARS
{
    public static void main(String[] args) throws Exception
    {
        distanceBetweenLeaves("dataset_10328_11.txt");
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
