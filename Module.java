public class Module implements Comparable<Module>{
    private int size;
    private int start;
    public Module(int size, int start){
        this.size = size;
        this.start = start;
    }

    public int getSize() {
        return size;
    }

    public int getStart() {
        return start;
    }

    @Override
    public int compareTo(Module o) {
        if (this.size < o.getSize()){
            return -1;
        }
        else if (this.size == o.getSize()){
            return 0;
        }
        else{
            return 1;
        }
    }
}
