import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

/**
 * Created by Alena on 21.12.2014.
 */
public class SequenceAlignment
{
    public static final int DOWN = 1;
    public static final int RIGHT = 2;
    public static final int DIAGONAL = 3;

    public static void main(String[] args) throws Exception
    {
        Scanner scanner = new Scanner(new File("dataset_245_7.txt"));
        ArrayList<DirectedEdge> edges = new ArrayList<DirectedEdge>();
        int numVertices = 0;
        int source = Integer.parseInt(scanner.nextLine().trim());
        int sink = Integer.parseInt(scanner.nextLine().trim());
        while (scanner.hasNextLine())
        {
            String line = scanner.nextLine();
            String[] fromTo = line.split("->");
            int from = Integer.parseInt(fromTo[0]);
            String[] toDist = fromTo[1].split(":");
            int to = Integer.parseInt(toDist[0]);
            int dist = Integer.parseInt(toDist[1]);
            DirectedEdge directedEdge = new DirectedEdge(from, to, dist);
            int max = Math.max(from, to);
            numVertices = Math.max(max, numVertices);
            edges.add(directedEdge);
        }
        EdgeWeightedDigraph digraph = new EdgeWeightedDigraph(numVertices + 1);
        for (DirectedEdge edge : edges)
        {
            digraph.addEdge(edge);
        }

        AcyclicLP acyclicLP = new AcyclicLP(digraph, source);
        System.out.println((int) acyclicLP.distTo(sink));
        Iterable<DirectedEdge> path = acyclicLP.pathTo(sink);
        for (DirectedEdge directedEdge : path)
        {
            System.out.print(directedEdge.from() + "->");
        }
        System.out.print(sink);
    }

    private static void calcLCS()
    {
        String v = "A";
        String w = "A";
        int[][] backtrack = backTrackLCS(v, w);
        outputLCS(backtrack, v, v.length() - 1, w.length() - 1);
    }

    public static void outputLCS(int[][] backtrack, String v, int i, int j)
    {
        if (i < 0 || j < 0)
        {
            return;
        }
        if (backtrack[i][j] == DOWN)
        {
            outputLCS(backtrack, v, i - 1, j);
        } else if (backtrack[i][j] == RIGHT)
        {
            outputLCS(backtrack, v, i, j - 1);
        } else if (backtrack[i][j] == DIAGONAL)
        {
            outputLCS(backtrack, v, i - 1, j - 1);
            System.out.print(v.charAt(i));
        }
    }

    public static int[][] backTrackLCS(String v, String w)
    {
        int[][] backtrack = new int[v.length()][w.length()];
        int s[][] = new int[v.length() + 1][w.length() + 1];
        for (int i = 1; i <= v.length(); i++)
        {
            for (int j = 1; j <= w.length(); j++)
            {
                int diagonal = s[i - 1][j - 1];
                if (v.charAt(i - 1) == w.charAt(j - 1))
                {
                    diagonal++;
                }
                s[i][j] = Math.max(Math.max(s[i - 1][j], s[i][j - 1]), diagonal);
                if (s[i][j] == s[i - 1][j])
                {
                    backtrack[i - 1][j - 1] = DOWN;
                } else if (s[i][j] == s[i][j - 1])
                {
                    backtrack[i - 1][j - 1] = RIGHT;
                } else if (s[i][j] == s[i - 1][j - 1] + 1)
                {
                    backtrack[i - 1][j - 1] = DIAGONAL;
                }
            }
        }
        return backtrack;
    }

    private static void calcTopologicalOrdering() throws FileNotFoundException
    {
        Scanner scanner = new Scanner(new File("dataset_254_2.txt"));
        ArrayList<Pair> pairs = new ArrayList<Pair>();
        int numVertices = 0;
        while (scanner.hasNextLine())
        {
            String line = scanner.nextLine();
            String[] fromTo = line.split(" -> ");
            int from = Integer.parseInt(fromTo[0]);
            String[] tos = fromTo[1].split(",");
            for (String toString : tos)
            {
                int to = Integer.parseInt(toString);
                Pair pair = new Pair(from, to);
                int max = Math.max(from, to);
                numVertices = Math.max(max, numVertices);
                pairs.add(pair);
            }
        }
        Digraph digraph = new Digraph(numVertices + 1);
        for (Pair pair : pairs)
        {
            digraph.addEdge(pair.from, pair.to);
        }
        Iterable<Integer> topo = topologicalOrdering(digraph);
        for (Integer integer : topo)
        {
            System.out.print(integer + ", ");
        }
    }

    static class Pair
    {
        int from;
        int to;

        public Pair(int from, int to)
        {
            this.from = from;
            this.to = to;
        }
    }

    private static void calcManhattanTourist() throws FileNotFoundException
    {
        Scanner scanner = new Scanner(new File("dataset_261_9.txt"));
        String[] nm = scanner.nextLine().split(" ");
        int n = Integer.parseInt(nm[0]);
        int m = Integer.parseInt(nm[1]);
        int[][] down = new int[n][m + 1];
        int[][] right = new int[n + 1][m];
        boolean isDown = true;
        int i = 0;
        while (scanner.hasNextLine())
        {
            String line = scanner.nextLine();
            if (line.equals("-"))
            {
                isDown = false;
                i = 0;
                continue;
            }
            String[] nums = line.split(" ");
            for (int j = 0; j < nums.length; j++)
            {
                int num = Integer.parseInt(nums[j]);
                if (isDown)
                {
                    down[i][j] = num;
                } else
                {
                    right[i][j] = num;
                }
            }
            i++;
        }
        System.out.println(manhattanTourist(n, m, down, right));
    }

    private static void calcDPChange()
    {
        System.out.println(DPChange(19489, new int[]{24, 18, 17, 16, 5, 3, 1}));
    }

    public static Iterable<Integer> topologicalOrdering(Digraph digraph)
    {
        Topological topological = new Topological(digraph);
        if (topological.hasOrder())
        {
            return topological.order();
        } else
        {
            System.out.println("the input graph is not a DAG");
            return null;
        }
    }

    private static Set<Integer> setOfAllNodesWithNoIncomingEdges(Digraph digraph)
    {
        Set<Integer> nodes = new HashSet<Integer>();
        for (int i = 0; i < digraph.V(); i++)
        {
            nodes.add(i);
        }
        for (int i = 0; i < digraph.V(); i++)
        {
            for (Integer wittIncomingEdge : digraph.adj(i))
            {
                nodes.remove(wittIncomingEdge);
            }
        }
        return nodes;
    }

    public static int DPChange(int money, int[] coins)
    {
        ArrayList<Integer> minNumCoins = new ArrayList<Integer>();
        minNumCoins.add(0, 0);
        for (int m = 1; m <= money; m++)
        {
            minNumCoins.add(m, Integer.MAX_VALUE);
            for (int i = 0; i < coins.length; i++)
            {
                int coin = coins[i];
                if (m >= coin)
                {
                    if (minNumCoins.get(m - coin) + 1 < minNumCoins.get(m))
                    {
                        minNumCoins.set(m, minNumCoins.get(m - coin) + 1);
                    }
                }
            }
        }
        return minNumCoins.get(money);
    }

    public static int manhattanTourist(int n, int m, int[][] down, int[][] right)
    {
        int[][] result = new int[n + 1][m + 1];
        result[0][0] = 0;
        for (int i = 1; i <= n; i++)
        {
            result[i][0] = result[i - 1][0] + down[i - 1][0];
        }
        for (int i = 1; i <= m; i++)
        {
            result[0][i] = result[0][i - 1] + right[0][i - 1];
        }

        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= m; j++)
            {
                int fromUp = result[i - 1][j] + down[i - 1][j];
                int fromLeft = result[i][j - 1] + right[i][j - 1];
                result[i][j] = Math.max(fromUp, fromLeft);
            }
        }
        return result[n][m];
    }
}