package rayTracerDistAccelShdPhtnMap;
/////
///		Final ray tracer from cs7490 - supports distribution RT, acceleration structure(BVH), perlin and worley-noise-based textures, photon mapping
/////

import processing.core.*;

import java.awt.event.KeyEvent;
import java.util.*;

//type of myGeomBase object
enum objType {
	None(0),
	BBox(1),
	RenderedBBox(2),
	Instance(3),
	PointLight(4),
	SpotLight(5),
	DiskLight(6),
	Triangle(7),
	Quad(8),
	Plane(9),
	Sphere(10),
	Cylinder(11),
	Hollow_Cylinder(12),
	Torus(13),
	AccelFlatList(14),
	AccelBVH(15);	
	private int value; 
	private static Map<Integer, objType> map = new HashMap<Integer, objType>(); 
	static { for (objType enumV : objType.values()) { map.put(enumV.value, enumV);}}
	private objType(int _val){value = _val;} 
	public int getVal(){return value;}
	public static objType getVal(int idx){return map.get(idx);}
	public static int getNumVals(){return map.size();}						//get # of values in enum
}

public class DistRayTracer extends PApplet {

	public final int sceneCols = 300;
	public final int sceneRows = 300;
	
	//used for determining refinement array
	public final int[] pow2 = new int[]{1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768};
	//used for worley txtrs
	public final int[][] nghbrHdCells = new int[][]{
		{ 0,  0,  0},{ 0,  0,  1},{ 0,  0, -1},{ 0,  1,  0},{ 0,  1,  1},{ 0,  1, -1},{ 0, -1,  0},{ 0, -1,  1},{ 0, -1, -1},
		{ 1,  0,  0},{ 1,  0,  1},{ 1,  0, -1},{ 1,  1,  0},{ 1,  1,  1},{ 1,  1, -1},{ 1, -1,  0},{ 1, -1,  1},{ 1, -1, -1},
		{-1,  0,  0},{-1,  0,  1},{-1,  0, -1},{-1,  1,  0},{-1,  1,  1},{-1,  1, -1},{-1, -1,  0},{-1, -1,  1},{-1, -1, -1}		
	};
	public final double log10_2 = Math.log10(2.0);//div to get base 2 log
	
	public final int hashPrime1 = 1572869;
	public final int hashPrime2 = 6291469;
	
	//map holding all loaded scene descriptions - scene should describe all scene and rendering-specific variables and quantities
	public TreeMap<String, myScene> loadedScenes;

	//max # of prims per leaf of accel structure 
	public final int maxPrimsPerLeaf = 5;
	
	public final float sqrt3 = PApplet.sqrt(3.0f),
			sqrt66 = PApplet.sqrt(6.0f)/6.0f, 
			sqrt612 = .5f*sqrt66;

	//name to save the file with, folder to put picture file in, a global variable for holding current active file name.
	public String gCurrentFile;
	public String currentDir;	//current directory to search for cli file
	//file reader/interpreter
	public myRTFileReader rdr; 	
	
	public boolean[] flags;
	//interface flags	
	public final int altKeyPressed  	= 0;			//alt pressed
	public final int cntlKeyPressed  	= 1;			//cntrl pressed
	public final int numFlags = 2;
	
	public final double epsVal = .0000001;

	//photon kd tree stuff
	public double photon_radius = 2;
	
	public void settings(){	size(sceneCols,sceneRows, P3D);	}	
	public void setup() {
		colorMode(RGB, 1.0f);
		background(0, 0, 0);
		initProg();
		initialize_table();
		gCurrentFile = "";
		currentDir = "";		//TODO set different directories
	}
	
	public void initProg(){	rdr = new myRTFileReader(this);	loadedScenes = new TreeMap<String, myScene>();	initBoolFlags();}//initProg method
	public void draw() {if(!gCurrentFile.equals("")){loadedScenes.get(gCurrentFile).draw();}}
	public void initBoolFlags(){
		flags = new boolean[numFlags];
		for (int i = 0; i < numFlags; ++i) { flags[i] = false;}	
	}		
	
	//print out info at the end of rendering
	public void DispEnd(){
		System.out.println("");
		System.out.println("Image rendered : " +gCurrentFile + " Current directory for cli's : " + getDirName());
		System.out.println("");
	}
	
	//key press IO -ugh.  change to UI
	public void keyPressed() {
		switch(key) {
			case '`' : {if(!gCurrentFile.equals("")){ loadedScenes.get(gCurrentFile).flipNormal();}   break;}//flip poly norms
			case '1':  {gCurrentFile = currentDir + "t01.cli"; break;}		//photon map cli's
			case '2':  {gCurrentFile = currentDir + "t02.cli"; break;}
			case '3':  {gCurrentFile = currentDir + "t03.cli"; break;}
			case '4':  {gCurrentFile = currentDir + "t04.cli"; break;}
			case '5':  {gCurrentFile = currentDir + "t05.cli"; break;}
			case '6':  {gCurrentFile = currentDir + "t06.cli"; break;}
			case '7':  {gCurrentFile = currentDir + "t07.cli"; break;}
			case '8':  {gCurrentFile = currentDir + "t08.cli"; break;}
			case '9':  {gCurrentFile = currentDir + "t09.cli"; break;}
			case '0':  {gCurrentFile = currentDir + "t10.cli"; break;}
			case '-':  {gCurrentFile = currentDir + "t11.cli"; break;}
			
		    case '!':  {gCurrentFile = currentDir + "p4_st01.cli"; break;}		//worley textures on spheres
		    case '@':  {gCurrentFile = currentDir + "p4_st02.cli"; break;}
		    case '#':  {gCurrentFile = currentDir + "p4_st03.cli"; break;}
		    case '$':  {gCurrentFile = currentDir + "p4_st04.cli"; break;}
		    case '%':  {gCurrentFile = currentDir + "p4_st05.cli"; break;}
		    case '^':  {gCurrentFile = currentDir + "p4_st06.cli"; break;}
		    case '&':  {gCurrentFile = currentDir + "p4_st07.cli"; break;}	    
		    case '*':  {gCurrentFile = currentDir + "p4_st08.cli"; break;}		    
		    case 'l':  {gCurrentFile = currentDir + "p4_st09.cli"; break;}		    
		    
		    case ';':  {gCurrentFile = currentDir + "p4_t09.cli"; break;}		    
		    
			case 'O':  {gCurrentFile = currentDir + "p4_t05.cli"; break;}		//wood bunny texture
			case 'P':  {gCurrentFile = currentDir + "p4_t06.cli"; break;}		//marble bunny texture
			case '(':  {gCurrentFile = currentDir + "p4_t07.cli"; break;}		//worley circles bunny texture
			case ')':  {gCurrentFile = currentDir + "p4_t08.cli"; break;}		//crackle bunny texture		

			case '_':  {gCurrentFile = currentDir + "p3_t10.cli"; break;}
			case '=':  {gCurrentFile = currentDir + "plnts3ColsBunnies.cli";  break;}
		    case '+':  {gCurrentFile = currentDir + "p3_t02_sierp.cli"; break;}		//TODO INVESTIGATE THIS
		    
		    case 'A':  {gCurrentFile = currentDir + "p2_t01.cli"; break;}	//these are images from project 2
		    case 'S':  {gCurrentFile = currentDir + "p2_t02.cli"; break;}
		    case 'D':  {gCurrentFile = currentDir + "p2_t03.cli"; break;}
		    case 'F':  {gCurrentFile = currentDir + "p2_t04.cli"; break;}
		    case 'G':  {gCurrentFile = currentDir + "p2_t05.cli"; break;}
		    case 'H':  {gCurrentFile = currentDir + "p2_t06.cli"; break;}
		    case 'J':  {gCurrentFile = currentDir + "p2_t07.cli"; break;}
		    case 'K':  {gCurrentFile = currentDir + "p2_t08.cli"; break;}
		    case 'L':  {gCurrentFile = currentDir + "p2_t09.cli"; break;}
		    case ':':  {gCurrentFile = currentDir + "old_t07c.cli"; break;}
		    
			case 'a':  {gCurrentFile = currentDir + "earthAA1.cli";   break;}
			case 's':  {gCurrentFile = currentDir + "earthAA2.cli";   break;}
			case 'd':  {gCurrentFile = currentDir + "earthAA3.cli";   break;}
			case 'f':  {gCurrentFile = currentDir + "c2clear.cli";   break;}
			case 'g':  {gCurrentFile = currentDir + "c3shinyBall.cli";   break;}
			case 'h':  {gCurrentFile = currentDir + "c4InSphere.cli";   break;}
			case 'j':  {gCurrentFile = currentDir + "c6.cli";   break;}
			case 'k':  {gCurrentFile = currentDir + "c6Fish.cli";   break;}	
			
		    case 'Q':  {gCurrentFile = currentDir + "c2torus.cli"; break;}			    
		    case 'W':  {gCurrentFile = currentDir + "old_t02.cli"; break;}//this is the most recent block of images for project 1b
		    case 'E':  {gCurrentFile = currentDir + "old_t03.cli"; break;}
		    case 'R':  {gCurrentFile = currentDir + "old_t04.cli"; break;}
		    case 'T':  {gCurrentFile = currentDir + "old_t05.cli"; break;}
		    case 'Y':  {gCurrentFile = currentDir + "old_t06.cli"; break;}
		    case 'U':  {gCurrentFile = currentDir + "old_t07.cli"; break;}
		    case 'I':  {gCurrentFile = currentDir + "old_t08.cli"; break;}
			case '{':  {gCurrentFile = currentDir + "old_t09.cli"; break;}
			case '}':  {gCurrentFile = currentDir + "old_t10.cli"; break;}
			
			case 'q':  {gCurrentFile = currentDir + "planets.cli"; break; }
			case 'w':  {gCurrentFile = currentDir + "planets2.cli"; break;}
			case 'e':  {gCurrentFile = currentDir + "planets3.cli"; break;}
			case 'r':  {gCurrentFile = currentDir + "planets3columns.cli";   break;}
			case 't' : {gCurrentFile = currentDir + "trTrans.cli";   break;}
			case 'y':  {gCurrentFile = currentDir + "planets3Ortho.cli"; break;}			
			case 'u':  {gCurrentFile = currentDir + "c1.cli";  break;}
			case 'i':  {gCurrentFile = currentDir + "c2.cli";  break;}
			case 'o':  {gCurrentFile = currentDir + "c3.cli";  break;}
			case 'p':  {gCurrentFile = currentDir + "c4.cli";  break;}
			case '[':  {gCurrentFile = currentDir + "c5.cli";  break;}
			case ']':  {gCurrentFile = currentDir + "c0.cli";  break;}
		    
			case 'Z':  {gCurrentFile = currentDir + "p3_t01.cli"; break;}
			case 'X':  {gCurrentFile = currentDir + "p3_t02.cli"; break;}
			case 'C':  {gCurrentFile = currentDir + "p3_t03.cli"; break;}
			case 'V':  {gCurrentFile = currentDir + "p3_t04.cli"; break;}
			case 'B':  {gCurrentFile = currentDir + "p3_t05.cli"; break;}
			case 'N':  {gCurrentFile = currentDir + "p3_t06.cli"; break;}
			case 'M':  {gCurrentFile = currentDir + "p3_t07.cli"; break;}
			case '<':  {gCurrentFile = currentDir + "p4_t06_2.cli"; break;}
			case '>':  {gCurrentFile = currentDir + "p4_t09.cli"; break;}
			case '?':  {gCurrentFile = currentDir + "p3_t11_sierp.cli"; break;}		//my bunny scene		
			case 'z':  {gCurrentFile = currentDir + "cylinder1.cli";  break;}
			case 'x':  {gCurrentFile = currentDir + "tr0.cli";   break;}
			case 'c':  {gCurrentFile = currentDir + "c0Square.cli";  break;}
			case 'v':  {gCurrentFile = currentDir + "c1octo.cli";  break;}
			case 'b':  {gCurrentFile = currentDir + "old_t0rotate.cli";  break;}
		    case 'n':  {gCurrentFile = currentDir + "old_t03a.cli"; break;}	//this block contains the first set of images given with assignment 1b 
		    case 'm':  {gCurrentFile = currentDir + "old_t04a.cli"; break;}
		    case ',':  {gCurrentFile = currentDir + "old_t05a.cli"; break;}
		    case '.':  {gCurrentFile = currentDir + "old_t06a.cli"; break;}
		    case '/':  {gCurrentFile = currentDir + "old_t07a.cli"; break;}			
			default : {return;}
		}//switch
		if(!gCurrentFile.equals("")){
			rdr.readRTFile(gCurrentFile, null);//pass null as scene so that we don't add to an existing scene
		}
		if((!flags[altKeyPressed])&&(key==CODED)){setFlags(altKeyPressed,(keyCode  == KeyEvent.VK_ALT));}
		if((!flags[cntlKeyPressed])&&(key==CODED)){setFlags(cntlKeyPressed,(keyCode  == KeyEvent.VK_CONTROL));}
	}	
	
	public void keyReleased(){
		if((flags[altKeyPressed])&&(key==CODED)){ if(keyCode == KeyEvent.VK_ALT){endAltKey();}}
		if((flags[cntlKeyPressed])&&(key==CODED)){ if(keyCode == KeyEvent.VK_CONTROL){endCntlKey();}}
	}		
	public void endAltKey(){clearFlags(new int []{altKeyPressed});	}
	public void endCntlKey(){clearFlags(new int []{cntlKeyPressed});}
	
	public void setFlags(int idx, boolean val ){
		flags[idx] = val;
		switch (idx){
			case altKeyPressed 		: { break;}//anything special for altKeyPressed 	
			case cntlKeyPressed 		: { break;}//anything special for altKeyPressed 	
		}
	}
	public void clearFlags(int[] idxs){		for(int idx : idxs){flags[idx]=false;}	}			
	
	public String getDirName(){	if(currentDir.equals("")){	return "data/";}	return "data/"+currentDir;}
	public void setCurDir(int input){
		switch (input){
		case 0 :{currentDir = "";break;}
		case 1 :{currentDir = "old/";break;}
		case 2 :{currentDir = "project1_1/";break;}
		case 3 :{currentDir = "project1_2/";break;}
		case 4 :{currentDir = "project2_1/";break;}
		case 5 :{currentDir = "project2_2/";break;}
		case 6 :{currentDir = "project3/";break;}
		case 7 :{currentDir = "project4/";break;}
		case 8 :{currentDir = "project5/";break;}
		default :{currentDir = "";break;}
		}
	}//setCurDir
	
	public double noise_3d(myVector pt){return noise_3d((float)pt.x, (float)pt.y, (float)pt.z);}
	//from code given by greg for project 4, 3d perlin noise
	public float noise_3d(float x, float y, float z) {		
		// make sure we've initilized table
		if (init_flag == false) {	  initialize_table();	  init_flag = true;	}		
		// Find unit grid cell containing point
		int X = fastfloor(x),Y = fastfloor(y), Z = fastfloor(z);		
		// Get relative xyz coordinates of point within that cell
		x = x - X;	y = y - Y;	z = z - Z;		
		// Wrap the integer cells at 255 (smaller integer period can be introduced here)
		X = X & 255;	Y = Y & 255;	Z = Z & 255;		
		// Calculate a set of eight hashed gradient indices
		int gi000 = perm[X+perm[Y+perm[Z]]] % 12;
		int gi001 = perm[X+perm[Y+perm[Z+1]]] % 12;
		int gi010 = perm[X+perm[Y+1+perm[Z]]] % 12;
		int gi011 = perm[X+perm[Y+1+perm[Z+1]]] % 12;
		int gi100 = perm[X+1+perm[Y+perm[Z]]] % 12;
		int gi101 = perm[X+1+perm[Y+perm[Z+1]]] % 12;
		int gi110 = perm[X+1+perm[Y+1+perm[Z]]] % 12;
		int gi111 = perm[X+1+perm[Y+1+perm[Z+1]]] % 12;
		
		// The gradients of each corner are now:
		// gXXX = grad3[giXXX];
		
		// Calculate noise contributions from each of the eight corners
		float n000= dot(grad3[gi000], x, y, z);
		float n100= dot(grad3[gi100], x-1, y, z);
		float n010= dot(grad3[gi010], x, y-1, z);
		float n110= dot(grad3[gi110], x-1, y-1, z);
		float n001= dot(grad3[gi001], x, y, z-1);
		float n101= dot(grad3[gi101], x-1, y, z-1);
		float n011= dot(grad3[gi011], x, y-1, z-1);
		float n111= dot(grad3[gi111], x-1, y-1, z-1);
		
		// Compute the fade curve value for each of x, y, z
		float u = fade(x), v = fade(y), w = fade(z);
		
//		// Interpolate along x the contributions from each of the corners
//		float nx00 = mix(n000, n100, u);
//		float nx01 = mix(n001, n101, u);
//		float nx10 = mix(n010, n110, u);
//		float nx11 = mix(n011, n111, u);
//		
//		// Interpolate the four results along y
//		float nxy0 = mix(nx00, nx10, v);
//		float nxy1 = mix(nx01, nx11, v);
//		
//		// Interpolate the two last results along z
//	
//		return mix(nxy0, nxy1, w);
		return mix(mix(mix(n000, n100, u), mix(n010, n110, u), v), mix(mix(n001, n101, u), mix(n011, n111, u), v), w);
	
	}//noise_3d

	public boolean init_flag = false;
	public int grad3[][] = {{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1},{0,1,1},{0,-1,1},{0,1,-1},{0,-1,-1}};
	public int p[] = {151,160,137,91,90,15,
				131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
				190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
				88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
				77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
				102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
				135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
				5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
				223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
				129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
				251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
				49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
				138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180};

	// To remove the need for index wrapping, double the permutation table length
	public int perm[] = new int[512];
	public void initialize_table() { for(int i=0; i<512; i++) perm[i]=p[i & 255];}
	// This method is a *lot* faster than using (int)Math.floor(x)
	public int fastfloor(float x) { return x>0 ? (int)x : (int)x-1;}
	public int fastfloor(double x) { return x>0 ? (int)x : (int)x-1;}
	public float dot(int g[], float x, float y, float z) { return g[0]*x + g[1]*y + g[2]*z;}
	public float mix(float a, float b, float t) { return (1-t)*a + t*b;}
	public float fade(float t) { return t*t*t*(t*(t*6-15)+10);}
	//end given code, 3d perlin noise
	
	//utility functions
	//find signed area of enclosed poly, with given points and normal N
	public double calcArea(myVector[] _pts, myVector N){
	    double res = 0;	    
	    if (_pts.length < 3) return 0; 
	    myVector absN = absVec(N);			//find abs val of each coord - want to project to plane normal to biggest coord to find area, and then scale by biggest coord	
	    if((absN.x > absN.y) && (absN.x > absN.z)){//x is max coord
	        for (int i=1, j=2, k=0; i<_pts.length; ++i, ++j, ++k){     res += (_pts[i].y * (_pts[(j%_pts.length)].z - _pts[k].z));}
	        res += (_pts[0].y * (_pts[1].z - _pts[_pts.length-1].z));	     
	        return (res / (2.0f * N.x));    		    	
	    } else if ((absN.y > absN.x) && (absN.x > absN.z)){//y is max coord
	    	for (int i=1, j=2, k=0; i<_pts.length; ++i, ++j, ++k){     res += (_pts[i].z * (_pts[(j%_pts.length)].x - _pts[k].x));}
	        res += (_pts[0].z * (_pts[1].x - _pts[_pts.length-1].x));
	        return (res / (2.0f * N.y));
	    } else {//z is max coord
	    	for (int i=1, j=2, k=0; i<_pts.length; ++i, ++j, ++k){     res += (_pts[i].x * (_pts[(j%_pts.length)].y - _pts[k].y));}
	        res += (_pts[0].x * (_pts[1].y - _pts[_pts.length-1].y));
	        //return (res * (1 / (2.0f * N.z)));   	
	        return (res / (2.0f * N.z));   	
	    }
    }//area	

	//rotate v1 around axis unit vector u, by give angle thet
	public myVector rotVecAroundAxis(myVector v1, myVector u, double thet){		
		double cThet = Math.cos(thet), sThet = Math.sin(thet), oneMC = 1-cThet,
				ux2 = u.x*u.x, uy2 = u.y*u.y, uz2 = u.z*u.z,
				uxy = u.x * u.y, uxz = u.x * u.z, uyz = u.y*u.z,
				uzS = u.z*sThet, uyS = u.y*sThet, uxS = u.x*sThet,
				uxzC1 = uxz *oneMC, uxyC1 = uxy*oneMC, uyzC1 = uyz*oneMC;
		//build rot matrix in vector def
		myVector res = new myVector(
				(ux2*oneMC+cThet) * v1.x + (uxyC1-uzS) 		* v1.y + (uxzC1+uyS) *v1.z,
				(uxyC1+uzS) 	  * v1.x + (uy2*oneMC+cThet)* v1.y + (uyzC1-uxS) *v1.z,
				(uxzC1-uyS) 	  * v1.x + (uyzC1+uxS)		* v1.y + (uz2*oneMC + cThet) * v1.z);
		
		return res;		
	}//rotVecAroundAxis
	
	
	//expand passed bbox to hold passed point - point is in box coords
	public void expandBoxPt(myBBox tarBox, myVector newPt) {
		tarBox.minVals.x = (tarBox.minVals.x < newPt.x) ?tarBox.minVals.x : newPt.x; 
		tarBox.minVals.y = (tarBox.minVals.y < newPt.y) ?tarBox.minVals.y : newPt.y; 
		tarBox.minVals.z = (tarBox.minVals.z < newPt.z) ?tarBox.minVals.z : newPt.z; 
		tarBox.maxVals.x = (tarBox.maxVals.x > newPt.x) ?tarBox.maxVals.x : newPt.x; 
		tarBox.maxVals.y = (tarBox.maxVals.y > newPt.y) ?tarBox.maxVals.y : newPt.y; 
		tarBox.maxVals.z = (tarBox.maxVals.z > newPt.z) ?tarBox.maxVals.z : newPt.z; 
		tarBox.calcMinMaxCtrVals(tarBox.minVals,tarBox.maxVals);
	}

	//expand bbox to encompass passed box
	public void expandBoxByBox(myBBox tarBox, myBBox srcBox, gtMatrix fwdTrans) {
		expandBoxPt(tarBox, getTransformedPt(srcBox.minVals,fwdTrans));
		expandBoxPt(tarBox, getTransformedPt(srcBox.maxVals,fwdTrans));
	}
	public void expandBoxByBox(myBBox tarBox, myBBox srcBox) {
		expandBoxPt(tarBox, srcBox.minVals);
		expandBoxPt(tarBox, srcBox.maxVals);
	}
	//expand bbox by delta in all dir
	public void expandBoxDel(myBBox tarBox, double delta) {
		myVector delVec = new myVector(tarBox.minVals);
		delVec._sub(delta, delta, delta);
		expandBoxPt(tarBox, delVec);
		delVec = new myVector(tarBox.maxVals);
		delVec._add(delta, delta, delta);
		expandBoxPt(tarBox, delVec);		
	}
	//point needs to be in box space(transformed via box's ctm)
	public boolean pointIsInBox(myBBox tarBox, myVector pt){return (((tarBox.minVals.x < pt.x) && ( pt.x < tarBox.maxVals.x)) && 
													((tarBox.minVals.y < pt.y) && ( pt.y < tarBox.maxVals.y)) && 
													((tarBox.minVals.z < pt.z) && ( pt.z < tarBox.maxVals.z)));}

	//get transformed/inverse transformed point - homogeneous coords
	public myVector getTransformedPt(myVector pt, gtMatrix trans){
		double[] newPtAra = trans.multVert(new double[]{pt.x, pt.y, pt.z, 1});	
		myVector newPt = new myVector(newPtAra[0],newPtAra[1],newPtAra[2]);
		return newPt;
	}	
	//get transformed/inverse transformed vector - homogeneous coords
	public myVector getTransformedVec(myVector vec, gtMatrix trans){
		double[] newVecAra = trans.multVert(new double[]{vec.x, vec.y, vec.z, 0});		
		myVector newVec = new myVector(newVecAra[0],newVecAra[1],newVecAra[2]);
		return newVec;
	}	
	//inv mat idx : 1; transpose mat idx : 2; adjoint mat idx : 3
	private gtMatrix[] buildMatExt(gtMatrix[] CTMara){CTMara[1] = CTMara[0].inverse();CTMara[2] = CTMara[0].transpose();CTMara[3] = CTMara[1].transpose();return CTMara;}
	//passing Mat so that can instance transformed prims like distorted spheres
	public gtMatrix[] buildCTMara(myScene scene, gtMatrix _mat){gtMatrix[] CTMara = new gtMatrix[4];CTMara[0] = _mat.multMat(scene.matrixStack.peek());return buildMatExt(CTMara);}	
	//rebuild mat ara such that passed matrix _mat1 to be fwd transformed by _baseMat.  
	public gtMatrix[] reBuildCTMara(gtMatrix _mat1, gtMatrix _prntMat){gtMatrix[] CTMara = new gtMatrix[4];CTMara[0] = _prntMat.multMat(_mat1);return buildMatExt(CTMara);}	
	public gtMatrix[] buildCTMara(myScene scene){gtMatrix[] CTMara = new gtMatrix[4];CTMara[0] = scene.matrixStack.peek(); return buildMatExt(CTMara);}	
	public gtMatrix[] buildIdentCTMara(){gtMatrix[] CTMara = new gtMatrix[4];CTMara[0] = new gtMatrix(); return buildMatExt(CTMara);	}
	
	
	//returns idx (0-2) of coord of max variance in array of arraylists for bvh
	public int getIDXofMaxBVHSpan(List<myGeomBase>[] _tmpCtrObjList){
		double maxSpan = -1, diff;
		int idxOfMaxSpan = -1, araSize = _tmpCtrObjList[0].size();
		for(int i = 0; i<3; ++i){
			diff = _tmpCtrObjList[i].get(araSize-1).trans_origin[i] - _tmpCtrObjList[i].get(0).trans_origin[i] ; 
			//System.out.println("Diff : arasize : " + araSize + " diff : " + diff);
			if(maxSpan < diff){		maxSpan = diff;	idxOfMaxSpan = i;}				
		}
		return idxOfMaxSpan;
	}//getIDXofMaxSpan

	//build a color value from a string array read in from a cli file.  stIdx is position in array where first color resides
	public myColor readColor(String[] token, int stIdx){return new myColor(Double.parseDouble(token[stIdx]),Double.parseDouble(token[stIdx+1]),Double.parseDouble(token[stIdx+2]));}
	
	//max and min of array of doubles
	public double max(double[] valAra) {double maxVal = -Double.MAX_VALUE;for (double val : valAra){	if(val > maxVal){maxVal = val;}	}return maxVal;}
	public double min(double[] valAra) {double minVal = Double.MAX_VALUE;for (double val : valAra){	if(val < minVal){minVal = val;}	}return minVal;}
	//return abs vals of vector as vector
	public myVector absVec(myVector _v){return new myVector(Math.abs(_v.x),Math.abs(_v.y),Math.abs(_v.z));}
	//interpolate 2 vectors, t==0 == a, t==1 == b
	public myVector interpVec(myVector a, double t, myVector b){
		myVector bMa = new myVector(a,b);
		return new myVector(a.x + t*bMa.x, a.y + t*bMa.y, a.z + t*bMa.z );
	}
	//vector functions (originally static, static member functions cannot be a part of an inner class, which processing treats all classes defined in pdes.	
	public myVector _mult(myVector p, double n){ myVector result = new myVector(p.x * n, p.y * n, p.z * n); return result;}                        //1 vec, 1 double
	public myVector _elemMult(myVector p, myVector q){ myVector result = new myVector(p.x *q.x, p.y * q.y, p.z * q.z); return result;}               //element-wise mult - 2 vec
	public void _mult(myVector p, myVector q, myVector r){ myVector result = new myVector(p.x *q.x, p.y * q.y, p.z * q.z); r.set(result);}       //2 vec src, 1 vec dest  	
	public myVector _add(myVector p, myVector q){ myVector result = new myVector(p.x + q.x, p.y + q.y, p.z + q.z); return result;}                //2 vec
	public void _add(myVector p, myVector q, myVector r){ myVector result = new myVector(p.x + q.x, p.y + q.y, p.z + q.z); r.set(result);}       //2 vec src, 1 vec dest  
	public myVector _sub(myVector p, myVector q){ myVector result = new myVector(p.x - q.x, p.y - q.y, p.z - q.z); return result;}                //2 vec
	public void _sub(myVector p, myVector q, myVector r){ myVector result = new myVector(p.x - q.x, p.y - q.y, p.z - q.z); r.set(result);}       //2 vec src, 1 vec dest  
	public myVector _normalize(myVector v){double magn = v._mag(); myVector newVec = (magn == 0) ? (new myVector(0,0,0)) : (new myVector( v.x /= magn, v.y /= magn, v.z /= magn)); newVec._mag(); return newVec;}

	public double _sqDist(myVector q, myVector r){  return ((r.x - q.x) *(r.x - q.x)) + ((r.y - q.y) *(r.y - q.y)) + ((r.z - q.z) *(r.z - q.z));}
	public double _dist(myVector q, myVector r){  return (double)Math.sqrt(_sqDist(q,r));}
	public double _angleBetween(myVector v1, myVector v2) {
		double 	_v1Mag = v1._mag(), 
				_v2Mag = v2._mag(), 
				dotProd = v1._dot(v2),
				cosAngle = dotProd/(_v1Mag * _v2Mag),
				angle = Math.acos(cosAngle);
		return angle;
	}//_angleBetween
	
	//get orthogonal vector to passed vector
	public myVector getOrthoVec(myVector vec){
		myVector tmpVec = new myVector(1,1,0);
		tmpVec._normalize();
		if(Math.abs(tmpVec._dot(vec) - 1) < epsVal ){		tmpVec.set(0,0,1); }	//colinear - want non-colinear vector for xprod
		myVector tmpRes = vec._cross(tmpVec);												//surface tangent vector - rotate this around normal by random amt, and extend from origin by random dist 0->radius			
		tmpRes._normalize();
		return tmpRes;
	}	
	
	//double versions of float disp functions
	public void ellipse(double a, double b, double c, double d){ellipse((float)a,(float)b,(float)c,(float)d);}
	//returns one of a set of predefined colors or a random color as an array of 0-1 doubles based on tag passed
	public myColor getClr(String colorVal){
		switch (colorVal.toLowerCase()){
			case "clr_rnd"				: { return new myColor(random(0,1),random(0,1),random(0,1));}
	    	case "clr_gray"   		    : { return new myColor(0.47,0.47,0.47);}
	    	case "clr_white"  		    : { return new myColor(1.0,1.0,1.0);}
	    	case "clr_yellow" 		    : { return new myColor(1.0,1.0,0);}
	    	case "clr_cyan"			    : { return new myColor(0,1.0,1.0);} 
	    	case "clr_magenta"		    : { return new myColor(1.0,0,1.0);}  
	    	case "clr_red"    		    : { return new myColor(1.0,0,0);} 
	    	case "clr_blue"			    : { return new myColor(0,0,1.0);}
	    	case "clr_purple"		    : { return new myColor(0.6,0.2,1.0);}
	    	case "clr_green"		    : { return new myColor(0,1.0,0);}  
	    	//lower idxs are darker
	    	case "clr_ltwood1"			: { return new myColor(0.94, 0.47, 0.12);}  
	    	case "clr_ltwood2"			: { return new myColor(0.94, 0.8, 0.4);}  
	    	
	    	case "clr_dkwood1"			: { return new myColor(0.2, 0.08, 0.08);}  
	    	case "clr_dkwood2"			: { return new myColor(0.3, 0.20, 0.16);}  
	    	
	    	case "clr_mortar1"			: { return new myColor(0.2, 0.2, 0.2);}
	    	case "clr_mortar2"			: { return new myColor(0.7, 0.7, 0.7);}

	    	case "clr_brick1_1"			: { return new myColor(0.6, 0.18, 0.22);}
	    	case "clr_brick1_2"			: { return new myColor(0.8, 0.26, 0.33);}
	    	
	    	case "clr_brick2_1"			: { return new myColor(0.6, 0.32, 0.16);}
	    	case "clr_brick2_2"			: { return new myColor(0.8, 0.45, 0.25);}
	    	
	    	case "clr_brick3_1"			: { return new myColor(0.3, 0.01, 0.07);}
	    	case "clr_brick3_2"			: { return new myColor(0.6, 0.02, 0.13);}
	    	
	    	case "clr_brick4_1"			: { return new myColor(0.4, 0.1, 0.17);}
	    	case "clr_brick4_2"			: { return new myColor(0.6, 0.3, 0.13);}
	    	
	    	case "clr_darkgray"   	    : { return new myColor(0.31,0.31,0.31);}
	    	case "clr_darkred"    	    : { return new myColor(0.47,0,0);}
	    	case "clr_darkblue"  	 	: { return new myColor(0,0,0.47);}
	    	case "clr_darkpurple"		: { return new myColor(0.4,0.2,0.6);}
	    	case "clr_darkgreen"  	    : { return new myColor(0,0.47,0);}
	    	case "clr_darkyellow" 	    : { return new myColor(0.47,0.47,0);}
	    	case "clr_darkmagenta"	    : { return new myColor(0.47,0,0.47);}
	    	case "clr_darkcyan"   	    : { return new myColor(0,0.47,0.47);}	  
	    	
	    	case "clr_lightgray"   	    : { return new myColor(0.78,0.78,0.78);}
	    	case "clr_lightred"    	    : { return new myColor(1.0,.43,.43);}
	    	case "clr_lightblue"   	    : { return new myColor(0.43,0.43,1.0);}
	    	case "clr_lightgreen"  	    : { return new myColor(0.43,1.0,0.43);}
	    	case "clr_lightyellow"	    : { return new myColor(1.0,1.0,.43);}
	    	case "clr_lightmagenta"	    : { return new myColor(1.0,.43,1.0);}
	    	case "clr_lightcyan"   	    : { return new myColor(0.43,1.0,1.0);}
	    	
	    	case "clr_black"		    : { return new myColor(0,0,0);}
	    	case "clr_nearblack"		: { return new myColor(0.05,0.05,0.05);}
	    	case "clr_faintgray" 		: { return new myColor(0.43,0.43,0.43);}
	    	case "clr_faintred" 	 	: { return new myColor(0.43,0,0);}
	    	case "clr_faintblue" 	 	: { return new myColor(0,0,0.43);}
	    	case "clr_faintgreen" 	    : { return new myColor(0,0.43,0);}
	    	case "clr_faintyellow" 	    : { return new myColor(0.43,0.43,0);}
	    	case "clr_faintcyan"  	    : { return new myColor(0,0.43,0.43);}
	    	case "clr_faintmagenta"  	: { return new myColor(0.43,0,0.43);}    	
	    	case "clr_offwhite"			: { return new myColor(0.95,0.98,0.92);}
	    	default         		    : { System.out.println("Color not found : " + colorVal + " so using white.");	return new myColor(1.0,1.0,1.0);}    
		}//switch
	}//getClr

}//rayParser class
