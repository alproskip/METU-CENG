package ceng.ceng351.cengtubedb;

import java.sql.*;
import java.sql.PreparedStatement;
import java.util.ArrayList;

public class CengTubeDB implements ICengTubeDB {
    
    private Connection conn = null;
    
    @Override
    public void initialize() {
        try {
            this.conn = DriverManager.getConnection("jdbc:mysql://144.122.71.168:3306/db223716","e223716","f16ed713");
        } catch (SQLException ex) {
            System.out.println("SQLException: "+ex.getMessage());
        }
    }

    @Override
    public int createTables() {
        int createdTables = 0;

        String UserTable = "CREATE TABLE IF NOT EXISTS user(" +
                "userID INT NOT NULL, " +
                "userName VARCHAR(30), " +
                "email VARCHAR(30), " +
                "password VARCHAR(30), " +
                "status VARCHAR(15), " +
                "PRIMARY KEY (userID));";

        String VideoTable = "CREATE TABLE IF NOT EXISTS video(" +
                "videoID INT NOT NULL, "+
                "userID INT, " +
                "videoTitle VARCHAR(60), " +
                "likeCount INT, " +
                "dislikeCount INT, " +
                "datePublished DATE, " +
                "FOREIGN KEY (userID) REFERENCES user(userID) ON DELETE CASCADE ON UPDATE CASCADE, " +
                "PRIMARY KEY (videoID));";

        String CommentTable = "CREATE TABLE IF NOT EXISTS comment(" +
                "commentID INT NOT NULL, " +
                "userID INT," +
                "videoID INT, " +
                "commentText VARCHAR(1000), " +
                "dateCommented DATE, " +
                "FOREIGN KEY (userID) REFERENCES video(videoID) ON DELETE SET NULL ON UPDATE CASCADE, " +
                "FOREIGN KEY (videoID) REFERENCES video(videoID) ON DELETE CASCADE ON UPDATE CASCADE, " +
                "PRIMARY KEY (commentID));";

        String WatchTable = "CREATE TABLE IF NOT EXISTS watch(" +
                "userID INT NOT NULL, " +
                "videoID INT NOT NULL, " +
                "dateWatched DATE, " +
                "FOREIGN KEY (userID) REFERENCES video(videoID) ON DELETE CASCADE ON UPDATE CASCADE, " +
                "FOREIGN KEY (videoID) REFERENCES video(videoID) ON DELETE CASCADE ON UPDATE CASCADE, " +
                "PRIMARY KEY (userID, videoID));";

        try { // CREATE USER TABLE
            Statement s = conn.createStatement();
            s.executeUpdate(UserTable);
            createdTables++;
            s.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: "+ex.getMessage());
        }

        try { // CREATE VIDEO TABLE
            Statement s = conn.createStatement();
            s.executeUpdate(VideoTable);
            createdTables++;
            s.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: "+ex.getMessage());
        }

        try { // CREATE COMMENT TABLE
            Statement s = conn.createStatement();
            s.executeUpdate(CommentTable);
            createdTables++;
            s.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: "+ex.getMessage());
        }

        try { // WATCH USER TABLE
            Statement s = conn.createStatement();
            s.executeUpdate(WatchTable);
            createdTables++;
            s.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: "+ex.getMessage());
        }

        return createdTables;
    }

    @Override
    public int dropTables() {
        int droppedTables = 0;

        String dropUser = "DROP TABLE IF EXISTS user;";
        String dropVideo = "DROP TABLE IF EXISTS video;";
        String dropComment = "DROP TABLE IF EXISTS comment;";
        String dropWatch = "DROP TABLE IF EXISTS watch;";

        try { // DROP COMMENT TABLE
            Statement s = conn.createStatement();
            s.executeUpdate(dropComment);
            droppedTables++;
            s.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: "+ex.getMessage());
        }

        try { // DROP WATCH TABLE
            Statement s = conn.createStatement();
            s.executeUpdate(dropWatch);
            droppedTables++;
            s.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: "+ex.getMessage());
        }

        try { // DROP VIDEO TABLE
            Statement s = conn.createStatement();
            s.executeUpdate(dropVideo);
            droppedTables++;
            s.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: "+ex.getMessage());
        }

        try { // DROP USER TABLE
            Statement s = conn.createStatement();
            s.executeUpdate(dropUser);
            droppedTables++;
            s.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: "+ex.getMessage());
        }

        return droppedTables;
    }

    @Override
    public int insertUser(User[] users) {
        int inserted = 0;

        for (User user : users) {
            try {
                PreparedStatement ps = conn.prepareStatement("INSERT INTO user VALUES(?,?,?,?,?);");
                ps.setInt(1,user.getUserID());
                ps.setString(2,user.getUserName());
                ps.setString(3,user.getEmail());
                ps.setString(4, user.getPassword());
                ps.setString(5, user.getStatus());
                ps.executeUpdate();
                ps.close();
                inserted++;
            } catch (SQLException ex) {
                System.out.println("SQLException: " + ex.getMessage());
            }
        }
        return inserted;
    }

    @Override
    public int insertVideo(Video[] videos) {
        int inserted = 0;

        for (Video video : videos) {
            try { // CREATE USER TABLE
                PreparedStatement ps = conn.prepareStatement("INSERT INTO video VALUES(?,?,?,?,?,?);");
                ps.setInt(1, video.getVideoID());
                ps.setInt(2, video.getUserID());
                ps.setString(3,video.getVideoTitle());
                ps.setInt(4, video.getLikeCount());
                ps.setInt(5, video.getDislikeCount());
                ps.setString(6, video.getDatePublished());
                ps.executeUpdate();
                ps.close();
                inserted++;
            } catch (SQLException ex) {
                System.out.println("SQLException: " + ex.getMessage());
            }
        }
        return inserted;
    }

    @Override
    public int insertComment(Comment[] comments) {
        int inserted = 0;

        for (Comment comment : comments) {
            try { // CREATE USER TABLE
                PreparedStatement ps = conn.prepareStatement("INSERT INTO comment VALUES(?,?,?,?,?);");
                ps.setInt(1, comment.getCommentID());
                ps.setInt(2, comment.getUserID());
                ps.setInt(3, comment.getVideoID());
                ps.setString(4, comment.getCommentText());
                ps.setString(5, comment.getDateCommented());
                ps.executeUpdate();
                ps.close();
                inserted++;
            } catch (SQLException ex) {
                System.out.println("SQLException: " + ex.getMessage());
            }
        }
        return inserted;
    }

    @Override
    public int insertWatch(Watch[] watchEntries) {
        int inserted = 0;

        for (Watch watch : watchEntries) {
            try { // CREATE USER TABLE
                PreparedStatement ps = conn.prepareStatement("INSERT INTO watch VALUES(?,?,?);");
                ps.setInt(1, watch.getUserID());
                ps.setInt(2, watch.getVideoID());
                ps.setString(3, watch.getDateWatched());
                ps.executeUpdate();
                ps.close();
                inserted++;
            } catch (SQLException ex) {
                System.out.println("SQLException: " + ex.getMessage());
            }
        }
        return inserted;
    }

    @Override
    public QueryResult.VideoTitleLikeCountDislikeCountResult[] question3() {

        ArrayList<QueryResult.VideoTitleLikeCountDislikeCountResult> resultList = new ArrayList<>();

        String query = "SELECT videoTitle, likeCount, dislikeCount FROM video WHERE likeCount > dislikeCount ORDER BY video.videoTitle;";
        try {
            Statement st = conn.createStatement();
            ResultSet rSet = st.executeQuery(query);

            while(rSet.next()){
                String videoTitle = rSet.getString("videoTitle");
                int likeCount = rSet.getInt("likeCount");
                int dislikeCount = rSet.getInt("dislikeCount");

                QueryResult.VideoTitleLikeCountDislikeCountResult temp = new QueryResult.VideoTitleLikeCountDislikeCountResult(videoTitle, likeCount, dislikeCount);
                resultList.add(temp);
            }
            st.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: " + ex.getMessage());
        }

        QueryResult.VideoTitleLikeCountDislikeCountResult[] result = new QueryResult.VideoTitleLikeCountDislikeCountResult[resultList.size()];
        return resultList.toArray(result);
    }

    @Override
    public QueryResult.VideoTitleUserNameCommentTextResult[] question4(Integer userID) {
        ArrayList<QueryResult.VideoTitleUserNameCommentTextResult> resultList = new ArrayList<>();
        String uID = userID.toString();
        String query = "SELECT videoTitle, userName, commentText "+
                "FROM user U, video V, comment C "+
                "WHERE U.userID = C.userID AND V.videoID = C.videoID AND U.userID = "+ uID + " ORDER BY V.videoTitle;";
        try {
            Statement st = conn.createStatement();
            ResultSet rSet = st.executeQuery(query);

            while(rSet.next()){
                String videoTitle = rSet.getString("videoTitle");
                String userName = rSet.getString("userName");
                String commentText = rSet.getString("commentText");

                QueryResult.VideoTitleUserNameCommentTextResult temp = new QueryResult.VideoTitleUserNameCommentTextResult(videoTitle, userName, commentText);
                resultList.add(temp);
            }
            st.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: " + ex.getMessage());
        }

        QueryResult.VideoTitleUserNameCommentTextResult[] result = new QueryResult.VideoTitleUserNameCommentTextResult[resultList.size()];

        return resultList.toArray(result);
    }

    @Override
    public QueryResult.VideoTitleUserNameDatePublishedResult[] question5(Integer userID) {
        ArrayList<QueryResult.VideoTitleUserNameDatePublishedResult> resultList = new ArrayList<>();
        String uID = userID.toString();
        String query = "SELECT V.videoTitle, U.userName, V.datePublished "+
                "FROM user U, video V "+
                "WHERE U.userID=V.userID AND " +
                "datePublished = (SELECT min(datePublished) FROM user U1, video V1 "+
                "WHERE U1.userID=V1.userID AND V1.videoTitle NOT LIKE '%VLOG%' AND U1.userID = " + uID + ") ORDER BY V.videoTitle;";
        try {
            Statement st = conn.createStatement();
            ResultSet rSet = st.executeQuery(query);

            while(rSet.next()){
                String videoTitle = rSet.getString("videoTitle");
                String userName = rSet.getString("userName");
                String datePublished = rSet.getString("datePublished");

                QueryResult.VideoTitleUserNameDatePublishedResult temp = new QueryResult.VideoTitleUserNameDatePublishedResult(videoTitle, userName, datePublished);
                resultList.add(temp);
            }
            st.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: " + ex.getMessage());
        }

        QueryResult.VideoTitleUserNameDatePublishedResult[] result = new QueryResult.VideoTitleUserNameDatePublishedResult[resultList.size()];

        return resultList.toArray(result);
    }

    @Override
    public QueryResult.VideoTitleUserNameNumOfWatchResult[] question6(String dateStart, String dateEnd) {

        ArrayList<QueryResult.VideoTitleUserNameNumOfWatchResult> resultList = new ArrayList<>();

        String query = "SELECT V.videoTitle, U.userName, COUNT(*) FROM user U, video V, watch W "+
                "WHERE U.userID = V.userID AND V.videoID=W.videoID AND "+
                "W.dateWatched<='" + dateEnd + "' AND W.dateWatched>='" + dateStart + "' GROUP BY W.videoID ORDER BY COUNT(*) DESC LIMIT 3;";
        try {
            Statement st = conn.createStatement();
            ResultSet rSet = st.executeQuery(query);

            while(rSet.next()){
                String videoTitle = rSet.getString("videoTitle");
                String userName = rSet.getString("userName");
                Integer numofWatch = rSet.getInt("COUNT(*)");

                QueryResult.VideoTitleUserNameNumOfWatchResult temp = new QueryResult.VideoTitleUserNameNumOfWatchResult(videoTitle, userName, numofWatch);
                resultList.add(temp);
            }
            st.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: " + ex.getMessage());
        }

        QueryResult.VideoTitleUserNameNumOfWatchResult[] result = new QueryResult.VideoTitleUserNameNumOfWatchResult[resultList.size()];

        return resultList.toArray(result);
    }

    @Override
    public QueryResult.UserIDUserNameNumOfVideosWatchedResult[] question7() {

        ArrayList<QueryResult.UserIDUserNameNumOfVideosWatchedResult> resultList = new ArrayList<>();

        String query = "SELECT U.userID, U.userName, COUNT(*) "+
                "FROM user U, video V WHERE EXISTS "+
                "(SELECT W.videoID FROM watch W WHERE W.videoID=V.videoID AND W.userID=U.userID) AND NOT EXISTS "+
                "(SELECT W.videoID FROM watch W WHERE W.videoID=V.videoID AND W.userID<>U.userID) "+
                "GROUP BY U.userID ORDER BY U.userID;";

        try {
            Statement st = conn.createStatement();
            ResultSet rs = st.executeQuery(query);

            while (rs.next()){
                Integer uID = rs.getInt("userID");
                String name = rs.getString("userName");
                Integer count = rs.getInt("COUNT(*)");

                QueryResult.UserIDUserNameNumOfVideosWatchedResult temp = new QueryResult.UserIDUserNameNumOfVideosWatchedResult(uID, name, count);
                resultList.add(temp);
            }
            st.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: "+ ex.getMessage());
        }


        QueryResult.UserIDUserNameNumOfVideosWatchedResult[] result = new QueryResult.UserIDUserNameNumOfVideosWatchedResult[resultList.size()];
        return resultList.toArray(result);
    }

    @Override
    public QueryResult.UserIDUserNameEmailResult[] question8() {
        ArrayList<QueryResult.UserIDUserNameEmailResult> resultList = new ArrayList<>();

        String query = "SELECT U.userID, U.userName, U.email "+
                "FROM user U "+
                "WHERE EXISTS "+
                "(SELECT V1.videoID FROM video V1 WHERE U.userID=V1.userID) "+
                "AND NOT EXISTS "+
                "(SELECT V.videoID FROM video V WHERE V.userID=U.userID "+
                "AND NOT EXISTS "+
                "(SELECT V.videoID FROM watch W, comment C "+
                "WHERE U.userID=V.userID "+
                "AND W.videoID = V.videoID AND W.userID=V.userID AND C.userID=U.userID AND C.videoID=V.videoID)) "+
                "ORDER BY U.userID;";

        try {
            Statement st = conn.createStatement();
            ResultSet rs = st.executeQuery(query);

            while (rs.next()){
                Integer uID = rs.getInt("userID");
                String name = rs.getString("userName");
                String email = rs.getString("email");

                QueryResult.UserIDUserNameEmailResult temp = new QueryResult.UserIDUserNameEmailResult(uID, name, email);
                resultList.add(temp);
            }
            st.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: "+ ex.getMessage());
        }

        QueryResult.UserIDUserNameEmailResult[] result = new QueryResult.UserIDUserNameEmailResult[resultList.size()];
        return resultList.toArray(result);
    }

    @Override
    public QueryResult.UserIDUserNameEmailResult[] question9() {

        ArrayList<QueryResult.UserIDUserNameEmailResult> resultList = new ArrayList<>();

        String query = "SELECT U.userID, U.userName, U.email "+
                "FROM user U WHERE EXISTS "+
                "(SELECT W.userID FROM watch W WHERE W.userID=U.userID) " +
                "AND NOT EXISTS "+
                "(SELECT C.userID FROM comment C WHERE C.userID=U.userID) "+
                "ORDER BY U.userID;";

        try {
            Statement st = conn.createStatement();
            ResultSet rs = st.executeQuery(query);

            while (rs.next()){
                Integer uID = rs.getInt("userID");
                String name = rs.getString("userName");
                String email = rs.getString("email");

                QueryResult.UserIDUserNameEmailResult temp = new QueryResult.UserIDUserNameEmailResult(uID, name, email);
                resultList.add(temp);
            }
            st.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: "+ ex.getMessage());
        }

        QueryResult.UserIDUserNameEmailResult[] result = new QueryResult.UserIDUserNameEmailResult[resultList.size()];
        return resultList.toArray(result);

    }

    @Override
    public int question10(int givenViewCount) {
        int result = 0;
        String givenCount = String.valueOf(givenViewCount);
        String query = "UPDATE user U "+
                "SET U.status='verified' "+
                "WHERE EXISTS (SELECT U.userName, COUNT(*) "+
                    "FROM video V, watch W WHERE W.videoID = V.videoID AND U.userID = V.userID "+
                    "GROUP BY U.userID HAVING COUNT(*) > "+ givenCount + " order by count(*) desc);";
        try {
            Statement st = conn.createStatement();
            result = st.executeUpdate(query);
            st.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: " + ex.getMessage());
        }

        return result;
    }

    @Override
    public int question11(Integer videoID, String newTitle) {
        int result = 0;

        String query = "UPDATE video SET videoTitle=\""+ newTitle +"\" WHERE videoID="+ videoID.toString() + ";";
        try {
            Statement st = conn.createStatement();
            result = st.executeUpdate(query);
            st.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: " + ex.getMessage());
        }

        return result;
    }

    @Override
    public int question12(String videoTitle) {
        int result = 0;
        String query = "DELETE FROM video WHERE videoTitle =\"" + videoTitle + "\";";
        String query2 = "SELECT COUNT(*) FROM video;";
        try {
            Statement st = conn.createStatement();
            st.executeUpdate(query);
            //ResultSet rs = st.executeQuery(query2);
            //result = rs.getInt("COUNT(*)");
            //ResultSetMetaData rsmd = rs.getMetaData();
            //String namex = rsmd.getColumnName(1);
            //System.out.println(namex);
            st.close();
        } catch (SQLException ex) {
            System.out.println("SQLException: " + ex.getMessage());
        }

        try {
            Statement st = conn.createStatement();
            ResultSet rs = st.executeQuery(query2);
            rs.next();
            result = rs.getInt("COUNT(*)");
            st.close();
        } catch (SQLException ex) {
            System.out.println(ex.getMessage());
        }

        return result;
    }
}
