import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * Created by Alena on 10/15/2015.
 */
public class Utilities
{
    /**
     * graph is encoded like this:
     * <p>
     * 0->4:11
     * 1->4:2
     * 2->5:6
     * 3->5:7
     * 4->0:11
     * 4->1:2
     * 4->5:4
     * 5->4:4
     * 5->3:7
     * 5->2:6
     *
     * @param fileName       name of file containing the graph
     * @param numLinesToSkip number of lines before the graph structure starts (probably with some other parameters)
     * @return EdgeWeightedDigraph
     * @throws FileNotFoundException
     */
    public static EdgeWeightedDigraph constructEdgeWeightedDigraphFromString(String fileName, int numLinesToSkip) throws FileNotFoundException
    {
        Scanner scanner = new Scanner(new File(fileName));

        //skipping lines
        for (int i = 0; i < numLinesToSkip; i++)
        {
            scanner.nextLine();
        }

        //creating list of Directed Edges first
        ArrayList<DirectedEdge> edges = new ArrayList<>();
        int numVertices = 0;
        while (scanner.hasNextLine())
        {
            String line = scanner.nextLine();
            //splitting from and to vertices. "0->4:11" will result in "0" and "4:11"
            String[] fromTo = line.split("->");
            //parsing from vertex
            int from = Integer.parseInt(fromTo[0]);
            //splitting to and distance. "4:11" will result in "4" and "11"
            String[] toDist = fromTo[1].split(":");
            int to = Integer.parseInt(toDist[0]);
            int dist = Integer.parseInt(toDist[1]);
            //constructing a DirectedEdge object which consists of from, to and distance
            DirectedEdge directedEdge = new DirectedEdge(from, to, dist);
            int max = Math.max(from, to);
            numVertices = Math.max(max, numVertices);
            edges.add(directedEdge);
        }

        //now adding ther vertices to a graph
        EdgeWeightedDigraph digraph = new EdgeWeightedDigraph(numVertices + 1);
        for (DirectedEdge edge : edges)
        {
            digraph.addEdge(edge);
        }

        return digraph;
    }

    public static Digraph constructDigraphFromString(String fileName, int numLinesToSkip) throws FileNotFoundException
    {
        Scanner scanner = new Scanner(new File(fileName));

        //skipping lines
        for (int i = 0; i < numLinesToSkip; i++)
        {
            scanner.nextLine();
        }

        //creating list of Directed Edges first
        ArrayList<Edge> edges = new ArrayList<>();
        int numVertices = 0;
        while (scanner.hasNextLine())
        {
            String line = scanner.nextLine();
            //splitting from and to vertices. "0->4:11" will result in "0" and "4:11"
            String[] fromTo = line.split("->");
            //parsing from vertex
            int from = Integer.parseInt(fromTo[0]);
            //splitting to and distance. "4:11" will result in "4" and "11"
            String[] toDist = fromTo[1].split(":");
            int to = Integer.parseInt(toDist[0]);
            Edge edge = new Edge(from, to);
            int max = Math.max(from, to);
            numVertices = Math.max(max, numVertices);
            edges.add(edge);
        }

        //now adding ther vertices to a graph
        Digraph digraph = new Digraph(numVertices + 1);
        for (Edge edge : edges)
        {
            digraph.addEdge(edge.from, edge.to);
        }

        return digraph;

    }
}
