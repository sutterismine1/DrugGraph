//When I first wrote this, only God and me understood how it works. Now only God does.
import java.awt.*;
import java.io.*;
import java.util.*;

public class DrugGraph {
	private static final double threshold = 0.7;
	public int moduleLen = 0;
	private int[] tempmodules;	//helper array so BFS can add a module to it without worrying about the amount of modules (initialized to capacity)
	public Module[] modules;
	private int[] tempBFTstarts; //same reason as tempmodules, needed a temporary array to have a fixed size before I find out the real size
	private int capacity;
	private Vertex[] vertices;
	private Vertex[] connectedVertices;
	private double[][] W; //weighted matrix created from (1-similarityMat[i][j])
	private int[][] A; //adjacency matrix based on results of the weighted matrix
	File f = new File("MSTPrimResult.tab");
	FileOutputStream output = new FileOutputStream(f);
	PrintWriter pw = new PrintWriter(output);
	public DrugGraph() throws FileNotFoundException {
		this(3000);
	}
	public DrugGraph(int capacity) throws FileNotFoundException {
		this.capacity = capacity;
		this.vertices = new Vertex[capacity];
		this.W = new double[capacity][capacity];
		this.A = new int[capacity][capacity];
		this.tempmodules = new int[capacity];
		this.tempBFTstarts = new int[capacity];
	}
	//Helper method to read data from a given matrix (in this case the source is hard coded into the method but this can be modified)
	//Creates an adjacency matrix from a list of drugs and a similarity matrix
	//
	//Param: Threshold - defaults to 0.7 and must be between 0 and 1
	// is used to create an adjancency matrix from given similarity matrix where (1-sim[i][j] <= threshold)
	//
	//Result is stored in instance variables (vertices, W, A)
	public void readData(double threshold) throws FileNotFoundException {
		File f = new File("dockedApproved.tab");
		String drugBankID;
		String genericName;
		String SMILES;
		String url;
		String drugGroups;
		double score;
		double[][] sim = new double[vertices.length][vertices.length];
		Scanner in = new Scanner(f);
		in.useDelimiter("\t|\n");
		for (int j = 0; j<6; j++){     //ignores the header at the top of input
			in.next();
		}
		int i = 0;
		while (in.hasNext()) {
			genericName = in.next();
			SMILES = in.next();
			drugBankID = in.next();
			url = in.next();
			drugGroups = in.next();
			score = Double.parseDouble(in.next());
			vertices[i] = new Vertex(drugBankID, genericName, SMILES, url, drugGroups, score, false, Double.POSITIVE_INFINITY);
			i++;
		}
		f = new File("sim_mat.tab");
		in = new Scanner(f);
		in.useDelimiter("\t|\n");
		for (int j = 0; j < vertices.length; j++) {
			for (int k = 0; k < vertices.length; k++) {
				sim[j][k] = Double.parseDouble(in.next());
			}
		}
		for (int j = 0; j < sim.length; j++) {
			for (int k = 0; k < sim.length; k++) {
				if (j!=k && 1-sim[j][k]<=threshold){
					W[j][k] = 1-sim[j][k];
				}
				else{
					W[j][k] = Double.POSITIVE_INFINITY;
				}
			}
		}
		for (int j = 0; j < sim.length; j++) {
			for (int k = 0; k < sim.length; k++) {
				if(W[j][k] != Double.POSITIVE_INFINITY){
					A[j][k] = 1;
				}
				else{
					A[j][k] = 0;
				}
			}
		}
	}
	//helper method to find adjacency array for one vertex
	private int find(Vertex[] a, int id){
		for (int i = 0; i < a.length; i++) {
			if (a[i].value == id){
				return i;
			}
		}
		return -1;
	}
	//breadth-first search starting at index i
	public void BFS(int i){
		int total = 0;
		Queue queue = new LinkedList();
		Vertex v = vertices[i];
		v.wasVisited = true;
		queue.add(v);
		while (!queue.isEmpty()){
			Vertex cur = (Vertex) queue.remove();
			if (cur.getModule() == -1) {
				cur.setModule(i);
			}
			total++;
			int curindex = find(vertices, cur.value);
			for (int j = 0; j < vertices.length; j++) {
				if (A[curindex][j] == 1 && vertices[j].wasVisited == false){	//for neighbours of u
					vertices[j].wasVisited = true;
					queue.add(vertices[j]);
				}
			}
		}
		for (Vertex vertex:
			 vertices) {
			vertex.wasVisited = false;
		}
		tempmodules[moduleLen] = (total);
		tempBFTstarts[moduleLen] = i;
		/*v.wasVisited = true;
		Q.enqueue(v);
		while !Q.isEmpty() // loop1
		{
			v = Q.dequeue();
		while there is an unvisited vertex u adjacent to v // loop2
		{
			u.wasVisited = true;
			Q.enqueue(u);
		}*/
	}
	//sets the module id of all vertices along a path to the given id (starts from i with BFT)
	private void setModuleId(int i, int id){
		int total = 0;
		Queue queue = new LinkedList();
		Vertex v = vertices[i];
		v.wasVisited = true;
		queue.add(v);
		while (!queue.isEmpty()){
			Vertex cur = (Vertex) queue.remove();
			cur.setModule(id);
			total++;
			int curindex = find(vertices, cur.value);
			for (int j = 0; j < vertices.length; j++) {
				if (A[curindex][j] == 1 && vertices[j].wasVisited == false){	//for neighbours of u
					vertices[j].wasVisited = true;
					queue.add(vertices[j]);
				}
			}
		}
		for (Vertex vertex:
				vertices) {
			vertex.wasVisited = false;
		}
	}
	//returns the number of disconnected components in the drug graph
	public int findModules(){
		for (int i = 0; i < vertices.length; i++) {
			if (vertices[i].getModule() == -1){
				BFS(i);
				moduleLen++;
			}
		}
		modules = new Module[moduleLen];
		for (int i = 0; i < moduleLen; i++) {
			modules[i] = new Module(tempmodules[i], tempBFTstarts[i]);
		}
		int total = 0;
		for (Vertex vertex:
			 vertices) {
			if (vertex.getModule() == 2){
				total++;
			}
		}
		Arrays.sort(modules);
		Module sortedModules[] = new Module[moduleLen];
		//Arrays.sort provides a sorted list in increasing order so we need to reverse the output array
		int j = 0;
		for (int i = moduleLen-1; i >= 0; i--) {
			sortedModules[j] = modules[i];
			setModuleId(modules[i].getStart(), j);
			System.out.println(modules[i].getSize());
			j++;
		}
		this.modules = sortedModules;
		return moduleLen;
	}
	public Vertex[] keepAModule(int moduleID){
		//loop through all the modules and return an array with only the vertices with proper module ID
		Vertex arr[] = new Vertex[modules[moduleID].getSize()];
		int i = 0;
		for (Vertex vertex:
			 vertices) {
			if (vertex.getModule() == moduleID){
				arr[i] = vertex;
				i++;
			}
		}
		return arr;
	}
	private void uw(int from, int to){
		int starti = (find(connectedVertices, from));
		Queue queue = new LinkedList();
		Vertex v = connectedVertices[starti];
		v.dist = 0;
		v.wasVisited = true;
		queue.add(v);
		while (!queue.isEmpty()){
			Vertex cur = (Vertex) queue.remove();
			cur.wasVisited = true;
			int curindex = find(vertices, cur.value);
			for (int j = 0; j < vertices.length; j++) {
				if (A[curindex][j] == 1){	//for neighbours of u
					if (vertices[j].dist > cur.dist + 1){	//if distance can be improved by going through v to w
						vertices[j].dist = cur.dist + 1;
						vertices[j].path = cur;
						queue.add(vertices[j]);
					}
				}
			}
		}
		//reset values of wasVisited
		for (Vertex vertex:
			 connectedVertices) {
			vertex.wasVisited = false;
		}
		int endi = (find(connectedVertices, to));
		//now that we have the dists and paths set iterate backwards from the end to the start to see the path
		String temparray[] = new String[(int)connectedVertices[endi].dist+1];
		Vertex curOutputVertex = connectedVertices[endi];
		for (int i = (int)curOutputVertex.dist; i >= 0; i--) {
			temparray[i] = curOutputVertex.drugBankID;
			curOutputVertex = curOutputVertex.path;
		}
		//now we have an array of the path, iterate through the array, printing the results
		for (int i = 0; i < temparray.length; i++) {
			String output = (i != temparray.length-1)? temparray[i] + " - ":temparray[i]; //only prints line if there is a next drug
			System.out.print(output);
		}
		System.out.println();
	}
	public void w(int from, int to){
		int starti = (find(connectedVertices, from));
		Queue queue = new PriorityQueue();
		Vertex v = connectedVertices[starti];
		v.dist = 0;
		queue.add(v);
		while (!queue.isEmpty()){
			Vertex cur = (Vertex) queue.remove();
			if (!cur.wasVisited) {
				cur.wasVisited = true;
				int curindex = find(vertices, cur.value);
				for (int j = 0; j < vertices.length; j++) {
					if (W[curindex][j] != Double.POSITIVE_INFINITY) {    //for neighbours of u
						if (vertices[j].dist > cur.dist + W[curindex][j]) {    //if distance can be improved by going through v to w
							vertices[j].dist = cur.dist + W[curindex][j];
							vertices[j].path = cur;
							queue.add(vertices[j]);
						}
					}
				}
			}
		}
		//reset values of wasVisited
		for (Vertex vertex:
				connectedVertices) {
			vertex.wasVisited = false;
		}
		int endi = (find(connectedVertices, to));
		int pathLength = 0;
		//we need to figure out how many vertices it takes to get from S to F
		Vertex path = connectedVertices[endi]; //start at the end and go back until it can't go back anymore
		while (path != null){
			path = path.path; //this may be the worst line of code I've ever seen
			pathLength++;
		}
		//make an array that has the path from S to F
		String temparray[] = new String[pathLength];
		Vertex curVertex = connectedVertices[endi];
		for (int i = pathLength-1; i >= 0; i--) {
			temparray[i] = curVertex.drugBankID;
			curVertex = curVertex.path;
		}
		//now we have an array of the path, iterate through the array, printing the results
		for (int i = 0; i < temparray.length; i++) {
			String output = (i != temparray.length-1)? temparray[i] + " - ":temparray[i]; //only prints line if there is a next drug
			System.out.print(output);
		}
		System.out.println();
	}
	public void findShortestPath(String fromDrug, String toDrug, String method){
		if (method == "unweighted"){
			uw(Integer.parseInt(fromDrug.split("DB")[1]), Integer.parseInt(toDrug.split("DB")[1]));
		}
		else{
			w(Integer.parseInt(fromDrug.split("DB")[1]), Integer.parseInt(toDrug.split("DB")[1]));
		}
		for (Vertex vertex:
			 connectedVertices) {
			vertex.path = null;
			vertex.dist = Double.POSITIVE_INFINITY;
			vertex.wasVisited = false;
		}
	}
	public double MSTPrim(){
		int i = 0;
		Vertex lastRemoved = null;
		Queue queue = new PriorityQueue();
		Vertex v = connectedVertices[i];
		v.dist = 0;
		queue.add(new Entry(v, v.dist));
		while (!queue.isEmpty()){
			Vertex cur = ((Entry) queue.remove()).getKey();
			if (!cur.wasVisited) {
				cur.wasVisited = true;
				int curindex = find(vertices, cur.value);
				for (int j = 0; j < vertices.length; j++) {
					if (W[curindex][j] != Double.POSITIVE_INFINITY) {    //for neighbours of u
						if ((vertices[j].dist > W[curindex][j]) && (!vertices[j].wasVisited)) {    //if distance can be improved by going through v to w
							vertices[j].dist = W[curindex][j];
							vertices[j].path = cur;
							queue.add(new Entry(vertices[j], vertices[j].dist));
						}
					}
				}
			}
		}
		double total = 0;
		for (Vertex ver: connectedVertices) {
			ver.displayDrug(pw);
			total += ver.dist;
		}
		return total;
	}
	public static void main(String[] args) throws FileNotFoundException {
		DrugGraph d = new DrugGraph(1932);
		d.readData(threshold);
		System.out.println(d.findModules() + " modules found.");
		d.connectedVertices = d.keepAModule(0);
		d.findShortestPath("DB01050", "DB00316", "unweighted");
		d.findShortestPath("DB01050", "DB00316", "weighted");
		System.out.print(Math.round(d.MSTPrim() * 100) / 100.0); //prints the total weight of the MST rounded to 2 decimal places
	}
}
