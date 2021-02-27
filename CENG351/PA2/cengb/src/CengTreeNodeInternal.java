import java.util.ArrayList;

public class CengTreeNodeInternal extends CengTreeNode
{
    private ArrayList<Integer> keys;
    private ArrayList<CengTreeNode> children;

    public CengTreeNodeInternal(CengTreeNode parent)
    {
        super(parent);

        // TODO: Extra initializations, if necessary.
        this.type = CengNodeType.Internal;
        keys = new ArrayList<>();
        children = new ArrayList<>();
    }
    public CengTreeNodeInternal(CengTreeNode parent, ArrayList<Integer> newKeys, ArrayList<CengTreeNode> newChildren){
        super(parent);
        type = CengNodeType.Internal;
        keys = newKeys;
        children = newChildren;
    }

    // GUI Methods - Do not modify
    public ArrayList<CengTreeNode> getAllChildren()
    {
        return this.children;
    }
    public Integer keyCount()
    {
        return this.keys.size();
    }
    public Integer keyAtIndex(Integer index)
    {
        if(index >= this.keyCount() || index < 0)
        {
            return -1;
        }
        else
        {
            return this.keys.get(index);
        }
    }
    // Extra Functions
    public ArrayList<Integer> getAllKeys(){
        return this.keys;
    }

    public void addKey(Integer key){
        int idx = 0;
        for (int i=0; i<keys.size();i++) {
            if (key > keyAtIndex(i)) idx = i + 1;
        }
        keys.add(idx, key);
    }

    public void addChild(Integer idx, CengTreeNode child){
        children.add(idx, child);
    }

    public Integer getChildInsertIndex(Integer key){
        int idx = 0;
        for (int i=0; i<keys.size();i++){
            if (key <= keyAtIndex(i)) {
                continue;
            }
            idx = i+1;
        }
        return idx+1;
    }

    public void deleteInRange(Integer start, Integer end) {
        keys.subList(start, end).clear();
        children.subList(start+1, end+1).clear();
    }
}
