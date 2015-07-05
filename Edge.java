
class Edge
{
    int from;
    int to;

    Edge(int from, int to)
    {
        this.from = from;
        this.to = to;
    }

    @Override
    public boolean equals(Object obj)
    {
        if(!(obj instanceof Edge) )
        {
            return false;
        }
        Edge edge = (Edge) obj;
        return from == edge.from && to == edge.to;
    }
}
