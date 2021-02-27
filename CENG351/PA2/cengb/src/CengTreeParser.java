import java.io.*;
import java.util.ArrayList;
import java.util.Locale;
import java.util.Scanner;

public class CengTreeParser
{
    public static ArrayList<CengVideo> parseVideosFromFile(String filename)
    {
        ArrayList<CengVideo> videoList = new ArrayList<CengVideo>();


        // You need to parse the input file in order to use GUI tables.
        // TODO: Parse the input file, and convert them into CengVideos
        File f = new File(filename);
        Scanner s = null;
        try {
            s = new Scanner(f);
            s.useDelimiter("\\Z");
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        while (true){
            assert s != null;
            if (!s.hasNextLine()) break;
            String[] line_arr = s.nextLine().split("\\|");
            Integer key = Integer.parseInt(line_arr[0]);
            String videoTitle = line_arr[1];
            String channelName = line_arr[2];
            String category = line_arr[3];
            videoList.add(new CengVideo(key, videoTitle, channelName, category));
        }
        return videoList;
    }

    public static void startParsingCommandLine() throws IOException
    {
        // TODO: Start listening and parsing command line -System.in-.
        // There are 4 commands:
        // 1) quit : End the app, gracefully. Print nothing, call nothing, just break off your command line loop.
        // 2) add : Parse and create the video, and call CengVideoRunner.addVideo(newlyCreatedVideo).
        // 3) search : Parse the key, and call CengVideoRunner.searchVideo(parsedKey).
        // 4) print : Print the whole tree, call CengVideoRunner.printTree().
        // Commands (quit, add, search, print) are case-insensitive.
        BufferedReader r = new BufferedReader(new InputStreamReader(System.in));
        String newline;
        label:
        while (true){
            newline = r.readLine();
            String[] line_arr = newline.split("\\|");
            String comm = line_arr[0].toLowerCase();
            switch (comm) {
                case "add":
                    Integer key = Integer.parseInt(line_arr[1]);
                    String videoTitle = line_arr[2];
                    String channelName = line_arr[3];
                    String category = line_arr[4];
                    CengVideoRunner.addVideo(new CengVideo(key, videoTitle, channelName, category));
                    break;
                case "search":
                    Integer parsedKey = Integer.parseInt(line_arr[1]);
                    CengVideoRunner.searchVideo(parsedKey);
                    break;
                case "print":
                    CengVideoRunner.printTree();
                    break;
                case "quit":
                    break label;
            }
        }
    }
}
