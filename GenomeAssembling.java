import java.io.*;
import java.util.*;

/**
 * Created by Alena on 01.12.2014.
 */
public class GenomeAssembling
{
    public static ArrayList<String> nodesMapping;

    public static void main(String[] args) throws Exception
    {
        //        deBruijnFromKmers("dataset_200_7.txt");
        //                        eulerianCycle("dataset_203_2.txt");
        //        eulerianPath("dataset_203_5.txt");

        stringReconstructionProblem("test.txt");
    }

    public static void stringReconstructionProblem(String fileName) throws IOException
    {
        ArrayList<String> list = new ArrayList<String>();
        Scanner scanner = new Scanner(new File(fileName));
        int fragmentLenght = Integer.parseInt(scanner.nextLine());
        while (scanner.hasNextLine())
        {
            list.add(scanner.nextLine());
        }
        StringBuilder sb = deBruijnFromKmers(list);
        BufferedWriter writer = null;
        String pathname = "temp_de_bruijn.txt";
        File file = new File(pathname);
        try
        {
            writer = new BufferedWriter(new FileWriter(file));
            writer.append(sb);
        } finally
        {
            if (writer != null)
            {
                writer.close();
            }
        }
        StringBuilder finalString = new StringBuilder();
        ArrayList<Integer> eulerianCycle = eulerianPath(pathname);
        for (int i = 0; i < eulerianCycle.size() - 1; i++)
        {
            Integer integer = eulerianCycle.get(i);
            String node = nodesMapping.get(integer);
            if (i == 0)
            {
                finalString.append(node);
            } else
            {
                finalString.append(node.substring(fragmentLenght - 2));
            }
        }
        System.out.println(finalString);
    }

    private static void deBruijnFromKmers(String fileName) throws FileNotFoundException
    {
        ArrayList<String> list = new ArrayList<String>();
        Scanner scanner = new Scanner(new File(fileName));
        while (scanner.hasNextLine())
        {
            list.add(scanner.nextLine());
        }
        deBruijnFromKmers(list);
    }

    public static void overlap(ArrayList<String> patterns)
    {
        Digraph digraph = new Digraph(patterns.size());
        for (String pattern : patterns)
        {
            for (String s : patterns)
            {
                if (suffix(pattern).equals(prefix(s)))
                {
                    digraph.addEdge(patterns.indexOf(pattern), patterns.indexOf(s));
                }
            }
        }

        for (String pattern : patterns)
        {
            Iterable<Integer> adj = digraph.adj(patterns.indexOf(pattern));
            for (Integer integer : adj)
            {
                System.out.println(pattern + " -> " + patterns.get(integer));
            }
        }
    }

    public static ArrayList<Integer> eulerianPath(String fileName) throws FileNotFoundException
    {
        Digraph digraph = readDigraphFromFileWithNames(fileName);
        int[] vertices = findStartAndStopVertexForEulerianPath(digraph);
        Integer startVertex = vertices[0];
        Integer endVertex = vertices[1];
        digraph.addEdge(endVertex, startVertex);
        ArrayList<Integer> cycle = eulerianCycle(digraph, startVertex);
        cycle = startCycleFromNewVertex(startVertex, cycle);
        printEulerianCycleWithNames(cycle);
        return cycle;
    }

    public static void eulerianCycle(String fileName) throws FileNotFoundException
    {
        Digraph digraph = readDigraphFromFileWithNums(fileName);
        printEulerianCycleWithNums(eulerianCycle(digraph, 0));
    }

    private static Digraph readDigraphFromFileWithNames(String fileName) throws FileNotFoundException
    {
        List<SymbolEdge> edges = new ArrayList<>();
        Set<String> vertices = new HashSet<>();
        Scanner scanner = new Scanner(new File(fileName));
        while (scanner.hasNextLine())
        {
            String line = scanner.nextLine();
            String[] fromTo = line.split(" -> ");
            String from = fromTo[0];
            vertices.add(from);
            String[] toList = fromTo[1].split(",");
            for (String to : toList)
            {
                vertices.add(to);
                edges.add(new SymbolEdge(from, to));
            }
        }

        nodesMapping = new ArrayList<>(vertices.size());
        for (String vertice : vertices)
        {
            nodesMapping.add(vertice);
        }

        Digraph digraph = new Digraph(vertices.size());
        for (SymbolEdge edge : edges)
        {
            digraph.addEdge(nodesMapping.indexOf(edge.getFrom()), nodesMapping.indexOf(edge.getTo()));
        }
        return digraph;
    }

    private static Digraph readDigraphFromFileWithNums(String fileName) throws FileNotFoundException
    {
        List<Edge> edges = new ArrayList<Edge>();
        Scanner scanner = new Scanner(new File(fileName));
        int max = 0;
        while (scanner.hasNextLine())
        {
            String line = scanner.nextLine();
            String[] fromTo = line.split(" -> ");
            Integer from = Integer.parseInt(fromTo[0]);
            max = Math.max(max, from);
            String[] toList = fromTo[1].split(",");
            for (String toString : toList)
            {
                Integer to = Integer.parseInt(toString);
                edges.add(new Edge(from, to));
            }
        }

        Digraph digraph = new Digraph(max + 1);
        for (Edge edge : edges)
        {
            digraph.addEdge(edge.from, edge.to);
        }
        return digraph;
    }

    public static ArrayList<Integer> eulerianCycle(Digraph digraph, int startVertex)
    {
        List<Edge>[] unusedEdges = new ArrayList[digraph.V()];

        for (int i = 0; i < unusedEdges.length; i++)
        {
            unusedEdges[i] = new ArrayList<>();
            Iterable<Integer> adj = digraph.adj(i);
            for (Integer integer : adj)
            {
                unusedEdges[i].add(new Edge(i, integer));
            }
        }

        ArrayList<Integer> nodesWithUnusedEdges = new ArrayList<Integer>();
        for (int i = 0; i < unusedEdges.length; i++)
        {
            nodesWithUnusedEdges.add(i);
        }

        ArrayList<Integer> currentCycle = new ArrayList<Integer>();
        currentCycle.add(startVertex);
        currentCycle = createCycle(startVertex, unusedEdges, currentCycle, true);

        for (int i = 0; i < unusedEdges.length; i++)
        {
            List<Edge> unusedEdge = unusedEdges[i];
            if (unusedEdge.isEmpty())
            {
                nodesWithUnusedEdges.remove(new Integer(i));
            }
        }

        while (nodesWithUnusedEdges.size() > 0)
        {
            int nodeWithUnusedEdges = getNodesWithUnusedEdgesWithinTheCycle(currentCycle, nodesWithUnusedEdges);
            currentCycle = createCycle(nodeWithUnusedEdges, unusedEdges, currentCycle, true);
            for (int i = 0; i < unusedEdges.length; i++)
            {
                List<Edge> unusedEdge = unusedEdges[i];
                if (unusedEdge.isEmpty())
                {
                    nodesWithUnusedEdges.remove(new Integer(i));
                }
            }
        }
        return currentCycle;
    }

    private static Integer getNodesWithUnusedEdgesWithinTheCycle(ArrayList<Integer> cycle, ArrayList<Integer> nodesWithUnusedEdges)
    {
        for (Integer nodeWithUnusedEdges : nodesWithUnusedEdges)
        {
            if (cycle.contains(nodeWithUnusedEdges))
            {
                return nodeWithUnusedEdges;
            }
        }
        return null;
    }

    private static int[] findStartAndStopVertexForEulerianPath(Digraph digraph)
    {
        int[] result = new int[2];
        int[] out = new int[digraph.V()];
        int[] in = new int[digraph.V()];
        fillOutVertexNumArray(digraph, out);
        Digraph reversedDigraph = digraph.reverse();
        fillOutVertexNumArray(reversedDigraph, in);
        for (int i = 0; i < in.length; i++)
        {
            int moreOut = out[i] - in[i];
            if (moreOut > 0)
            {
                result[0] = i;
            } else if (moreOut < 0)
            {
                result[1] = i;
            }
        }
        return result;
    }

    private static void fillOutVertexNumArray(Digraph digraph, int[] out)
    {
        for (int i = 0; i < digraph.V(); i++)
        {
            int numOut = 0;
            Iterator<Integer> iterator = digraph.adj(i).iterator();
            while (iterator.hasNext())
            {
                numOut++;
                iterator.next();
            }
            out[i] = numOut;
        }
    }

    private static ArrayList<Integer> createCycle(int startVertex, List<Edge>[] unusedEdges, ArrayList<Integer> cycle, boolean external)
    {
        if (external && cycle.size() > 1)
        {
            cycle = startCycleFromNewVertex(startVertex, cycle);
        }
        List<Edge> unusedEdgesForStartVertex = unusedEdges[startVertex];
        if (unusedEdgesForStartVertex.isEmpty())
        {
            return cycle;
        }
        Edge unusedEdge = unusedEdgesForStartVertex.get(0);
        unusedEdgesForStartVertex.remove(unusedEdge);
        cycle.add(unusedEdge.to);
        return createCycle(unusedEdge.to, unusedEdges, cycle, false);
    }

    private static ArrayList<Integer> startCycleFromNewVertex(int startVertex, ArrayList<Integer> cycle)
    {
        ArrayList<Integer> newCycle = new ArrayList<>();
        int startPos = cycle.indexOf(new Integer(startVertex));
        int size = cycle.size();
        for (int i = 0; i < size; i++)
        {
            int pos = i + startPos;
            if (pos >= size)
            {
                pos = i + startPos - size;
            }
            if (!(pos == size - 1))
            {
                newCycle.add(cycle.get(pos));
            }
        }
        newCycle.add(startVertex);
        cycle = newCycle;
        return cycle;
    }

    private static void printEulerianCycleWithNames(ArrayList<Integer> cycle)
    {
        StringBuilder sb = new StringBuilder();
        for (Integer integer : cycle)
        {
            sb.append(nodesMapping.get(integer)).append("->");
        }
        sb.delete(sb.length() - 2, sb.length());
        System.out.println(sb.toString());
    }

    private static void printEulerianCycleWithNums(ArrayList<Integer> cycle)
    {
        StringBuilder sb = new StringBuilder();
        for (Integer integer : cycle)
        {
            sb.append(integer).append("->");
        }
        sb.delete(sb.length() - 2, sb.length());
        System.out.println(sb.toString());
    }

    public static StringBuilder deBruijnFromKmers(ArrayList<String> stringComposition)
    {
        SortedSet<String> compositionSet = new TreeSet<>();
        for (String s : stringComposition)
        {
            compositionSet.add(prefix(s));
            compositionSet.add(suffix(s));
        }
        HashMap<String, Integer> compositionMap = new HashMap<>();
        HashMap<Integer, String> inverseCompositionMap = new HashMap<>();

        int num = 0;
        for (String s : compositionSet)
        {
            compositionMap.put(s, num);
            inverseCompositionMap.put(num, s);
            num++;
        }

        Digraph d = new Digraph(compositionSet.size());
        for (String s : stringComposition)
        {
            String kmer = prefix(s);
            String nextKmer = suffix(s);
            d.addEdge(compositionMap.get(kmer), compositionMap.get(nextKmer));
        }

        StringBuilder result = new StringBuilder();

        for (String s : compositionSet)
        {
            Iterable<Integer> adj = d.adj(compositionMap.get(s));
            StringBuilder sb = new StringBuilder();
            for (Integer integer : adj)
            {
                sb.append(inverseCompositionMap.get(integer)).append(",");
            }
            if (sb.length() > 0)
            {
                sb.deleteCharAt(sb.length() - 1);
                result.append(s).append(" -> ");
                result.append(sb).append("\n");
            }
        }
        return result;
    }

    public static void deBruijn(String text, int k)
    {

        ArrayList<String> stringComposition = stringComposition(text, k - 1);
        SortedSet<String> compositionSet = new TreeSet<String>(stringComposition);
        HashMap<String, Integer> compositionMap = new HashMap<String, Integer>();
        HashMap<Integer, String> inverseCompositionMap = new HashMap<Integer, String>();

        int num = 0;
        for (String s : compositionSet)
        {
            compositionMap.put(s, num);
            inverseCompositionMap.put(num, s);
            num++;
        }

        Digraph d = new Digraph(compositionSet.size());
        for (int i = 0; i < text.length() - k + 1; i++)
        {
            String kmer = text.substring(i, i + k - 1);
            String nextKmer = text.substring(i + 1, i + k);
            d.addEdge(compositionMap.get(kmer), compositionMap.get(nextKmer));
        }

        for (String s : compositionSet)
        {
            Iterable<Integer> adj = d.adj(compositionMap.get(s));
            StringBuilder sb = new StringBuilder();
            for (Integer integer : adj)
            {
                sb.append(inverseCompositionMap.get(integer)).append(",");
            }
            if (sb.length() > 0)
            {
                sb.deleteCharAt(sb.length() - 1);
                System.out.print("\n" + s + " -> ");
                System.out.print(sb);
            }
        }
    }


    public static String suffix(String s)
    {
        return s.substring(1, s.length());
    }

    public static String prefix(String s)
    {
        return s.substring(0, s.length() - 1);
    }

    public static ArrayList<String> stringComposition(String text, int k)
    {
        ArrayList<String> composition = new ArrayList<String>();
        for (int i = 0; i < text.length() - k + 1; i++)
        {
            composition.add(text.substring(i, i + k));
        }
        Collections.sort(composition);
        return composition;
    }

    public static String reconstructFromPath(ArrayList<String> path)
    {
        StringBuilder sb = new StringBuilder();
        sb.append(path.get(0));
        int k = path.get(0).length();
        for (int i = 1; i < path.size(); i++)
        {
            String s = path.get(i);
            sb.append(s.substring(k - 1, k));
        }
        return sb.toString();
    }
}
