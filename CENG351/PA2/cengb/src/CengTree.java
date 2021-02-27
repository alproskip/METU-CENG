import java.util.ArrayList;

public class CengTree
{
    public CengTreeNode root;
    // Any extra attributes...

    public CengTree(Integer order)
    {
        CengTreeNode.order = order;
        // TODO: Initialize the class
        root = new CengTreeNodeLeaf(null);
    }

    public void addVideo(CengVideo video)
    {
        // TODO: Insert Video to Tree
        CengTreeNodeLeaf Leaf = findLeaf(video.getKey()); //Find correct leaf
        Integer InsertIndex = Leaf.InsertIdx(video.getKey());
        Leaf.addNewVideo(InsertIndex, video); //Insert to leaf

        int videoCount = Leaf.videoCount();
        if (Leaf.getParent() == null) { // if parent is null then its root
            if (videoCount > 2*CengTreeNode.order) { //overflow happened at root | new level
                CengTreeNodeInternal newRoot = new CengTreeNodeInternal(null);

                Integer midIdx = videoCount/2;
                Integer InternalNodeKey = Leaf.videoKeyAtIndex(midIdx);
                newRoot.addKey(InternalNodeKey);

                ArrayList<CengVideo> l2videos = Leaf.getVideos(midIdx, videoCount);
                CengTreeNodeLeaf Leaf2 = new CengTreeNodeLeaf(newRoot, l2videos);

                Leaf.deleteInRange(midIdx, videoCount);
                Leaf.setParent(newRoot);

                newRoot.addChild(0,Leaf);
                newRoot.addChild(1,Leaf2);
                root = newRoot;
            }
        }
        else { //leaf's parent is not null
            if (videoCount > 2*CengTreeNode.order){ //Overflow in Leaf Node
                CengTreeNodeInternal parentNode = (CengTreeNodeInternal) Leaf.getParent();
                Integer midIdx = videoCount/2;
                Integer InternalNodeKey = Leaf.videoKeyAtIndex(midIdx);
                parentNode.addKey(InternalNodeKey);

                ArrayList<CengVideo> l2videos = Leaf.getVideos(midIdx, videoCount);
                CengTreeNodeLeaf Leaf2 = new CengTreeNodeLeaf(parentNode, l2videos);
                Leaf.deleteInRange(midIdx, videoCount);
                parentNode.addChild(parentNode.getChildInsertIndex(InternalNodeKey), Leaf2);

                while (true){
                    Integer keyCount = parentNode.keyCount();
                    if (keyCount <= 2*CengTreeNode.order) break;
                    int midIntIdx = keyCount/2;
                    Integer middleKey = parentNode.keyAtIndex(midIntIdx);

                    CengTreeNodeInternal parent_of_parent = (CengTreeNodeInternal) parentNode.getParent();

                    if (parent_of_parent==null){ //Current Internal Node is Root
                        CengTreeNodeInternal newParent = new CengTreeNodeInternal(null);
                        newParent.addKey(middleKey);

                        //Construct new Internal Node
                        ArrayList<Integer> curr_keys = parentNode.getAllKeys();
                        ArrayList<CengTreeNode> curr_children = parentNode.getAllChildren();

                        ArrayList<Integer> keys2 = new ArrayList<>(curr_keys.subList(midIntIdx+1,keyCount));
                        ArrayList<CengTreeNode> children2 = new ArrayList<>(curr_children.subList(midIntIdx+1,keyCount+1));

                        CengTreeNodeInternal Internal2 = new CengTreeNodeInternal(newParent, keys2, children2);
                        ArrayList<CengTreeNode> temp = Internal2.getAllChildren();
                        for (CengTreeNode cengTreeNode : temp) {
                            cengTreeNode.setParent(Internal2);
                        }
                        parentNode.deleteInRange(midIntIdx, keyCount);
                        parentNode.setParent(newParent);

                        // Add 2 Internal nodes to the parent
                        newParent.addChild(0,parentNode);
                        newParent.addChild(1,Internal2);
                        root = newParent;
                        break;
                    }
                    else { //Parent also has a parent
                        //will add to parent_of_parent
                        ArrayList<Integer> curr_keys = parentNode.getAllKeys();
                        ArrayList<CengTreeNode> curr_children = parentNode.getAllChildren();
                        ArrayList<Integer> keys2 = new ArrayList<>(curr_keys.subList(midIntIdx+1,keyCount));
                        ArrayList<CengTreeNode> children2 = new ArrayList<>(curr_children.subList(midIntIdx+1,keyCount+1));

                        CengTreeNodeInternal Internal2 = new CengTreeNodeInternal(parent_of_parent, keys2, children2);
                        ArrayList<CengTreeNode> temp = Internal2.getAllChildren();
                        for (CengTreeNode cengTreeNode : temp) {
                            cengTreeNode.setParent(Internal2);
                        }
                        parentNode.deleteInRange(midIntIdx, keyCount);
                        parent_of_parent.addKey(middleKey);
                        parent_of_parent.addChild(parent_of_parent.getChildInsertIndex(middleKey), Internal2);
                        parentNode = parent_of_parent;
                    }
                }
            }
        }
    }

    public ArrayList<CengTreeNode> searchVideo(Integer key)
    {
        // TODO: Search within whole Tree, return visited nodes.
        // Return null if not found.
        CengTreeNode curr = root;
        ArrayList<CengTreeNode> visited = new ArrayList<>();
        int numOfTab = 0;
        while (curr.getType()==CengNodeType.Internal){
            visited.add(curr);
            String tabs = "\t".repeat(numOfTab*2);
            System.out.println(tabs+"<index>");
            int size = ((CengTreeNodeInternal) curr).keyCount();
            int idx = 0;
            for (int i=0; i<size; i++) {
                Integer ckey = ((CengTreeNodeInternal) curr).keyAtIndex(i);
                System.out.println(tabs+ckey);
                if (ckey <= key) {
                    idx = i+1;
                }
            }
            System.out.println(tabs+"</index>");
            numOfTab++;
            curr = ((CengTreeNodeInternal) curr).getAllChildren().get(idx);
        }
        int size = ((CengTreeNodeLeaf) curr).videoCount();
        for (int i=0; i<size; i++){
            if (((CengTreeNodeLeaf) curr).videoKeyAtIndex(i).equals(key)){
                visited.add(curr);
                String record = "<record>"+((CengTreeNodeLeaf) curr).getVideo(i).fullName()+"</record>";
                String tabs = "\t".repeat(numOfTab*2);
                System.out.println(tabs+record);
                return visited;
            }
        }
        System.out.println("Could not find <"+key+">.");
        return null;
    }

    public void printTree()
    {
        // TODO: Print the whole tree to console
        CengTreeNode curr = root;
        printTreeHelper(curr, 0);
    }

    // Any extra functions...
    public CengTreeNodeLeaf findLeaf(Integer key){
        CengTreeNode curr = root;
        while(curr.getType() == CengNodeType.Internal) {
            int idx = 0;
            for (int i = 0; i<((CengTreeNodeInternal) curr).keyCount(); i++){
                if (((CengTreeNodeInternal) curr).keyAtIndex(i) < key) idx = i + 1;
                else break;
            }
            curr = ((CengTreeNodeInternal) curr).getAllChildren().get(idx);
        }
        return (CengTreeNodeLeaf) curr;
    }

    public void printTreeHelper(CengTreeNode curr, Integer numOfTabs){
        String tabs = "\t".repeat(numOfTabs*2);
        if (curr.getType()==CengNodeType.Internal){
            int size = ((CengTreeNodeInternal) curr).keyCount();
            System.out.println(tabs+"<index>");
            for (int i=0; i < size; i++){
                Integer ckey = ((CengTreeNodeInternal) curr).keyAtIndex(i);
                System.out.println(tabs+ckey);
            }
            System.out.println(tabs+"</index>");
            for (int i = 0; i < size+1; i++) {
                printTreeHelper(((CengTreeNodeInternal) curr).getAllChildren().get(i), numOfTabs + 1);
            }
        }
        else{
            int size = ((CengTreeNodeLeaf) curr).videoCount();
            System.out.println(tabs+"<data>");
            for (int i=0; i<size; i++){
                String record = "<record>"+((CengTreeNodeLeaf) curr).getVideo(i).fullName()+"</record>";
                System.out.println(tabs+record);
            }
            System.out.println(tabs+"</data>");
        }
    }
}
