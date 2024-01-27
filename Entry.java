public class Entry implements Comparable<Entry>{
    private Vertex key;
    private double value;

    public Entry(Vertex key, double value){
        this.key = key;
        this.value = value;
    }
    public double getValue() {
        return value;
    }

    public Vertex getKey() {
        return key;
    }

    @Override
    public int compareTo(Entry o) {
        if (this.value < o.getValue()){
            return -1;
        }
        else if(this.value == o.getValue()){
            return 0;
        }
        else{
            return 1;
        }
    }
}
