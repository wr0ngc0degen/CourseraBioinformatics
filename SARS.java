import java.io.File;
import java.util.Scanner;

public class SARS
{
    public static void main(String[] args) throws Exception
    {
        distanceBetweenLeaves("test.txt");
    }

    //Transforming Distance Matrices into Evolutionary Trees | Step 11
    private static void distanceBetweenLeaves(String fileName) throws Exception
    {
        Scanner scanner = new Scanner(new File(fileName));
        while (scanner.hasNext())
        {
            System.out.println(scanner.nextLine());
        }
    }
}
