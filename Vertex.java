import java.io.PrintWriter;

public class Vertex implements Comparable<Vertex> { //so I can compare the vertices in a priority queue
    public String drugBankID;
    public String genericName;
    public String SMILES;
    public String url;
    public String drugGroups;
    public double score;
    public int value;
    public boolean wasVisited;
    public double dist;
    public Vertex path;
    private int module = -1;
    public Vertex(String drugBankID, String genericName, String SMILES, String url, String drugGroups, double score, boolean wasVisited, double dist){
        this.drugBankID = drugBankID;
        this.genericName = genericName;
        this.SMILES = SMILES;
        this.url = url;
        this.drugGroups = drugGroups;
        this.score = score;
        this.value = Integer.parseInt(drugBankID.split("DB")[1]);//this line will take just the number part of the drugID to compare
        this.wasVisited = wasVisited;
        this.dist = dist;
        this.path = null;
    }
    public void displayDrug(){
        System.out.print(drugBankID + ": ");
        System.out.print(wasVisited + ": ");
        System.out.print(dist + ": ");
        System.out.print(path + ": ");
        System.out.print(genericName + " ");
        System.out.print(SMILES + " ");
        System.out.print(url + " ");
        System.out.print(drugGroups + " | ");
        System.out.println(score);
    }
    public void displayDrug(PrintWriter pw){
        pw.print(drugBankID + ": ");
        pw.print(dist + ": ");
        pw.println();
    }
    public int getValue(){
        return this.value;
    }
    public int getModule(){return this.module;}
    public void setModule(int module) {this.module = module;}

    @Override
    public int compareTo(Vertex o) {
        if (this.dist < o.dist){
            return -1;
        }
        else if (this.dist == o.dist){
            return 0;
        }
        else{
            return 1;
        }
    }
}