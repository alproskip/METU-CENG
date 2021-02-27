import java.util.ArrayList;

public class CengTreeNodeLeaf extends CengTreeNode
{
    private ArrayList<CengVideo> videos;
    // TODO: Any extra attributes

    public CengTreeNodeLeaf(CengTreeNode parent)
    {
        super(parent);

        // TODO: Extra initializations
        type = CengNodeType.Leaf;
        videos = new ArrayList<CengVideo>();
    }

    public CengTreeNodeLeaf(CengTreeNode parent, ArrayList<CengVideo> videos) {
        super(parent);
        type = CengNodeType.Leaf;
        this.videos = videos;
    }

    // GUI Methods - Do not modify
    public int videoCount()
    {
        return videos.size();
    }
    public Integer videoKeyAtIndex(Integer index)
    {
        if(index >= this.videoCount()) {
            return -1;
        } else {
            CengVideo video = this.videos.get(index);

            return video.getKey();
        }
    }
    // Extra Functions

    public void addNewVideo(Integer idx,CengVideo video) {
        videos.add(idx, video);
    }

    public CengVideo getVideo(Integer index){
        return videos.get(index);
    }

    public Integer InsertIdx(Integer key) {
        int idx = 0;
        int size = videoCount();
        if (size == 0) {return 0;}
        for (int i = 0; i<size; i++) {
            if (videos.get(i).getKey() < key){
                idx = i+1;
            }
            else {break;}
        }
        return idx;
    }

    public ArrayList<CengVideo> getVideos(Integer start, Integer end) {
        ArrayList<CengVideo> result = new ArrayList<>(videos.subList(start, end));
        return result;
    }

    public void deleteInRange(Integer start, Integer end) {
        videos.subList(start,end).clear();
    }
}
