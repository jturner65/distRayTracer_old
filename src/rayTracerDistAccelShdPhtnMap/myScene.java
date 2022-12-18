package rayTracerDistAccelShdPhtnMap;

import java.util.*;
import java.util.Map.Entry;
import java.util.concurrent.*;

import processing.core.*;


//class to hold all objects within a desired scene
public abstract class myScene {
	public DistRayTracer p;
	
	//multi-threaded stuff
	public ExecutorService th_exec;
	
	////
	//constants and control variables for a particular scene - make changeable via .cli reader
	//maximum height/depth of matrix stack
	public int matStackMaxHeight = 20;

	//row and col values of scene - scene dimensions in pxls
	public int sceneCols = 300;
	//public static final int sceneRows = 240;
	public int sceneRows = 300;
	//number of rays in play - including initial ray (recursive depth)
	public int numRays = 8;

	//origin of eye rays
	public myVector eyeOrigin;
	
	//halfway point in rows and cols
	protected double rayYOffset, rayXOffset;
	
	//epsilon value for double calcs; 
	public int objCnt = 0;
	
	/////////////////////////////
	//refining index list
	public int[] RefineIDX;
	public int curRefineStep;
	public ArrayDeque<String> srcFileNames;
	
	public String saveName, folderName;										//name of cli file used to describe this scene, save file, name of containing folder
	
	//the current texture to be used for subsequent objects for either their "top" or "bottom" as determined by their normal, the texture for the background, result image of rendering
	public PImage currTextureTop, currTextureBottom, currBkgTexture, rndrdImg;
	
	//an array list of all the objects in the scene : objList does not include lights, objAndLightList includes all lights
	public ArrayList<myGeomBase> allObjsToFind;
	//the following to facilitate faster light lookup and instancing
	public ArrayList<myGeomBase> objList;
	public ArrayList<myGeomBase> lightList; 	
	
	//objects to be put in accel structure - temp storage as list is being built
	public ArrayList<myGeomBase> tmpObjList;

	//named objects - to be instanced later
	public TreeMap<String, myGeomBase> namedObjs;
	public TreeMap<String, Integer> numInstances;
	
	//background texture-holding sphere.
	public mySphere mySkyDome;	
	
	//current number of lights, number of objects built, used for ID field in mySceneObject constructor and addNewLights, num non-light objects
	public int numLights, objCount, numNonLights, numPxls, numNamedObjs;
	//debug - count rays refracted, reflected
	public long refrRays = 0, reflRays = 0, globRayCount = 0;
	
	public boolean[] scFlags;			//boolean flags describing characteristics of this scene
	
	public static final int	
		//scene/global level flags
		debugIDX			= 0,		//enable debug functionality
		renderedIDX 		= 1,		//this scene has been rendered since any changes were made
		saveImageIDX		= 2,		//whether or not to save an image
		saveImgInDirIDX		= 3,		//save image inside specific image directory, rather than root
		simpleRefrIDX		= 4,		//whether scene should use simplified refraction (intended to minimize refactorization requirement in mySceneObject)
		
		flipNormsIDX		= 5,		//whether or not we should flip the normal directions in this scene
		hasDpthOfFldIDX			= 6,		//using depth of field
		showObjInfoIDX 		= 7,		//print out object info after rendering image
		addToTmpListIDX		= 8,		//add object to the temp list, so that the objects will be added to some accel structure.
		timeRndrIDX			= 9,		//time the length of rendering and display the results.
		glblTxtrdBkgIDX		= 10,		//whether the background is to be textured
		glblRefineIDX 		= 11,		//whether scene should be rendered using iterative refinement technique
		useFGColorIDX		= 12,		//whether or not to use a foregroundcolor in this scene

		//currenlty-loading object level flags
		glblTxtrdTopIDX		= 13,		//whether the currently loading object should be txtred on the top
		glblTxtrdBtmIDX		= 14,		//whether the currently loading object should be txtred on the bottom
	
		usePhotonMapIDX		= 15,		//whether to use a photon map for caustics, indirect illumination
		isCausticPhtnIDX	= 16,		//whether to use caustic photons or indirect illumination photons (hacky)
		isPhtnMapRndrdIDX	= 17,		//whether or not photons have been cast yet
		doFinalGatherIDX	= 18;		//whether to do final gather
	
	public static final int numFlags = 19;
	
	////////////////
	//indirect illumination/caustics stuff
	public myKD_Tree photonTree;		//TODO make two of these, one for caustics, one for indirect illum
	public int numGatherRays;
	public int numPhotons, kNhood;
	public float ph_max_near_dist;
	//recursive depth for photons
	public int numPhotonRays = 4;	
	//to correct for light power oddness
	public double causticsLightPwrMult  = 40.0,
			diffuseLightPwrMult = 8.0;
	
	//end indirect illumination/caustics stuff
	////////
	
	//////////////////
	// proc txtrs - colors and constants 	
	public int txtrType;				//set to be 0 : none; 1 : image; 2 : noise; 3 : wood; 4 : marble; 5 : stone/brick
	public double noiseScale;			//set if using noise-based proc txtr
	
	public myColor[] noiseColors = new myColor[]{new myColor(.7,.7,.7),
												new myColor(.2,.2,.2)};
	//how much weight each color should have - TODO
	public Double[] clrWts = new Double[]{.5,.5};
	//turbulence values
	public int numOctaves = 8;
	public double //turbScale = 128.0,						// equivalent to max scale/size value of turbulence
					turbMult = 1.0;							//multiplier of turbulence effect
	public double colorScale = 10.0, colorMult = .2;
	public myVector pdMult = new myVector(10,10,10);		//vector to hold multipliers for periodic textures
	public boolean rndColors,								//randomize colors a little bit; 
					useCustClrs,							//colors specified in cli file
					useFwdTrans;							//use fwd-transformed hit location (for meshes, don't use for prims)		
	public int numOverlays = 1;								//number of parallel txtrs used as overlays (semi transparent)
	//worley noise
	public double avgNumPerCell = 1, mortarThresh = .04;						//avg # of particles per cell
	public int numPtsDist = 2;								//# of closest points we take distance of to get color/effect
	public int distFunc = 1;								//distance function used for worley noise : 0-> manhattan, 1->euclid
	public int roiFunc = 1;									//function for determining region of interest in worley noise - 0 : nearest x sum, 1 : alternating linear sum, 2: ...
	//end proc txtrs
	/////////////////////////////	
		
	public Calendar now;
	//global values set by surface parameter
	public myColor currDiffuseColor, currAmbientColor, currSpecularColor,globCurPermClr, currKReflClr, backgroundColor;//, foregroundColor;
	public double currPhongExp, currKRefl,globRfrIdx,currKTrans, currDepth, lens_radius, lens_focal_distance;   //value representing the visible depth of a colloid (i.e. subsurface experiment) 
	
	public myPlane focalPlane;							//if Depth of Field scene, this is the plane where the lens is in focus
	
	//replaced by num rays per pixel
	public int numRaysPerPixel;
	
	//constant color for mask of fisheye
	public myColor blkColor = new myColor(0,0,0);
	
	public double maxDim, yStart, xStart, fishMult;			//compensate for # rows or # cols not being max - make sure projection is centered in non-square images
	
	//project 3 timer stuff - length of time to render
	public float renderTime;
	
	///////
	//transformation stack stuff
	//
	public myMatStack matrixStack;
	//current depth in matrix stack - starts at 0;
	public int currMatrixDepthIDX;	
	
	public myScene(DistRayTracer _p) {
		p = _p;
		now = Calendar.getInstance();
		folderName = "pics." +getDateTimeString(); 
		setImageSize(p.sceneCols,p.sceneRows);		
		
		allObjsToFind = new ArrayList<myGeomBase>();
		lightList = new ArrayList<myGeomBase>();
		objList = new ArrayList<myGeomBase>();
		
		srcFileNames = new ArrayDeque<String>();
		eyeOrigin = new myVector(0,0,0);
		gtInitialize();       															 //sets up matrix stack
		initFlags();
		scFlags[saveImageIDX] = true;    												//default to saving image
		scFlags[saveImgInDirIDX] = true;    											//save to timestamped directories, to keep track of changing images
		scFlags[showObjInfoIDX] = true;    												//default to showing info
		namedObjs = new TreeMap<String,myGeomBase>();
		numInstances = new TreeMap<String,Integer>();

		tmpObjList = new ArrayList<myGeomBase>();										//list to hold objects being built to be put into acceleration structure	
		initVars();
	}
	
	public myScene(myScene _old){//copy ctor, for when scene type is being set - only use when old scene is being discarded (shallow copy)
		p = _old.p;
		now = _old.now;
		folderName = _old.folderName;
		setImageSize(_old.sceneCols,_old.sceneRows);

		allObjsToFind = _old.allObjsToFind;
		lightList = _old.lightList;
		objList = _old.objList;
		
		srcFileNames = _old.srcFileNames;
		eyeOrigin = new myVector(0,0,0);
		eyeOrigin.set(_old.eyeOrigin);
		namedObjs = _old.namedObjs;
		numInstances = _old.numInstances;

		tmpObjList = _old.tmpObjList;
		
		gtInitialize();       															 //sets up matrix stack
		scFlags=new boolean[numFlags];for(int i=0;i<numFlags;++i){scFlags[i]=_old.scFlags[i];}
		initVars(_old);	
		matrixStack = _old.matrixStack;
	}//myScene
	
	public void initFlags(){scFlags=new boolean[numFlags];for(int i=0;i<numFlags;++i){scFlags[i]=false;}}
	
	//scene-wide variables set during loading of scene info from .cli file
	public void initVars(){
		currTextureTop = null;
		currTextureBottom = null;
		currBkgTexture = null;
		photonTree = null;
		numGatherRays = 0;
		numPhotons =0;
		kNhood = 0;
		ph_max_near_dist = 0;
		
		saveName = "";
		txtrType = 0;
		noiseScale = 1;
		numNonLights = 0;
		numLights = 0;
		objCount = 0;
		numNamedObjs = 0;
		
		backgroundColor = new myColor(0,0,0);
//		foregroundColor = new myColor(1,1,1);
		currDiffuseColor = new myColor(0,0,0);
		currAmbientColor = new myColor(0,0,0);
		currSpecularColor = new myColor(0,0,0);
		currPhongExp = 0;
		currKRefl = 0;
		currKReflClr = new myColor(0,0,0);
		currDepth = 0;		 
		currKTrans = 0;
		globRfrIdx = 0.0;
		globCurPermClr = new myColor(0,0,0);

		curRefineStep = 0;
		reflRays = 0;
		refrRays = 0;
		globRayCount = 0;
		focalPlane = new myPlane(this);
	}//initVars method
	
	//scene-wide variables set during loading of scene info from .cli file
	public void initVars(myScene _old){
		currTextureTop = _old.currTextureTop;
		currTextureBottom = _old.currTextureBottom;
		currBkgTexture = _old.currBkgTexture;
		saveName = _old.saveName;
		txtrType = _old.txtrType;
		noiseScale = _old.noiseScale;
		
		numGatherRays = _old.numGatherRays;
		photonTree = _old.photonTree;	
		numPhotons = _old.numPhotons;
		kNhood = _old.kNhood;
		ph_max_near_dist = _old.ph_max_near_dist;
		
		numNonLights = _old.numNonLights;
		numLights = _old.numLights;
		objCount = _old.objCount;
		numNamedObjs = _old.numNamedObjs;
		
		backgroundColor = _old.backgroundColor;
//		foregroundColor = _old.foregroundColor;
		currDiffuseColor = _old.currDiffuseColor;
		currAmbientColor = _old.currAmbientColor;
		currSpecularColor = _old.currSpecularColor;
		currPhongExp = _old.currPhongExp;
		currKRefl = _old.currKRefl;
		currKReflClr = _old.currKReflClr;
		currDepth = _old.currDepth;		 
		currKTrans = _old.currKTrans;
		globRfrIdx = _old.globRfrIdx;
		globCurPermClr = _old.globCurPermClr;

		curRefineStep = _old.curRefineStep;
		reflRays = _old.reflRays;
		refrRays = _old.refrRays;
		globRayCount = _old.globRayCount;

		lens_radius = _old. lens_radius;
		lens_focal_distance = _old.lens_focal_distance;		
		focalPlane = _old.focalPlane;		
	}//initVars from old scene method	
	public abstract void setSceneParams(double[] args);	

	public void startTmpObjList(){
		tmpObjList = new ArrayList<myGeomBase>();
		scFlags[addToTmpListIDX] = true;
	}
	//end building the accel struct - if is list just build arraylist object, otherwise build acceleration struct
	public void endTmpObjList(int lstType){
		scFlags[addToTmpListIDX] = false;
		myAccelStruct accelObjList = null;
		if(lstType == 0){//flat geomlist
			accelObjList = new myGeomList(this);
			//int objAdded = 1;
			for(myGeomBase obj : tmpObjList){		((myGeomList)accelObjList).addObj(obj);	}
		} else if(lstType == 1){//bvh tree
			System.out.println("begin adding to BVH structure - # objs : " + tmpObjList.size());
			accelObjList = new myBVH(this);
			List<myGeomBase>[] _tmpCtr_ObjList = ((myBVH)accelObjList).buildSortedObjAras(tmpObjList,-1);
			((myBVH)accelObjList).addObjList(_tmpCtr_ObjList, 0, _tmpCtr_ObjList[0].size()-1);
			System.out.println("");
			System.out.println("Done Adding to BVH structure");
		} else if(lstType == 2){//KDtree/octree TODO
			System.out.println("begin adding to octree structure TODO");
			System.out.println("Done Adding to octree structure TODO");
		}
		addObjectToScene(accelObjList);
	}//endTmpObjList
	
	//build miniTet - call recursively to build successively smaller tets
	//set shader for object based on current level - pastel bunnies
	private void setSierpShdr(int level, int maxLevel){
		float bVal = 1.0f - Math.min(1,(1.5f*level/maxLevel)), 
				rVal = 1.0f - bVal, 
				tmp = Math.min((1.2f*(level-(maxLevel/2)))/(1.0f*maxLevel),1), 
				gVal = (tmp*tmp);
		myColor cDiff = new myColor(Math.min(1,rVal+.5f), Math.min(1,gVal+.5f), Math.min(1,bVal+.5f));
		scFlags[glblTxtrdTopIDX]  = false;
		scFlags[glblTxtrdBtmIDX] = false;
		setSurface(cDiff,new myColor(0,0,0),new myColor(0,0,0),0,0);
	}
	//recursively build sierpenski tetrahedron - dim is relative dimension, decreases at each recursive call
	private void buildSierpSubTri(float dim, float scVal, String instName, int level, int maxLevel, boolean addShader){
		if(level>=maxLevel){return;}
		float newDim = scVal*dim;
		
		gtPushMatrix();
		gtTranslate(0,.1f*dim,0);
		gtRotate(70, 0, 1, 0);	
		if(addShader){	setSierpShdr(level, maxLevel);	}
		addInstance(instName,addShader);
		gtPopMatrix();
		//up bunny
		float newTrans = p.sqrt66*dim;
		gtPushMatrix();		
		gtTranslate(0,newTrans,0);
		gtScale(scVal,scVal,scVal);
		buildSierpSubTri(newDim, scVal, instName,level+1,maxLevel,addShader);
		gtPopMatrix();
		//front bunny
		gtPushMatrix();	
		sierpShiftObj(newTrans);
		gtScale(scVal,scVal,scVal);
		buildSierpSubTri(newDim, scVal,instName,level+1,maxLevel,addShader);
		gtPopMatrix();
		//left bunny
		gtPushMatrix();		
		gtRotate(120,0,1,0);
		sierpShiftObj(newTrans);
		gtRotate(-120,0,1,0);
		gtScale(scVal,scVal,scVal);
		buildSierpSubTri(newDim, scVal,instName,level+1,maxLevel,addShader);
		gtPopMatrix();		
		//right bunny
		gtPushMatrix();		
		gtRotate(-120,0,1,0);
		sierpShiftObj(newTrans);
		gtRotate(120,0,1,0);
		gtScale(scVal,scVal,scVal);
		buildSierpSubTri(newDim,scVal,instName,level+1,maxLevel,addShader);
		gtPopMatrix();		
	}
	
	private void sierpShiftObj(float newTrans){
		gtRotate(120,1,0,0);
		gtTranslate(0,newTrans,0);
		gtRotate(-120,1,0,0);		
	}
	//build a sierpinski tet arrangement using instances of object name
	//depth is how deep to build the tetrahedron, useShdr is whether or not to use changing shader based on depth
	public void buildSierpinski(String name, float scVal, int depth, boolean useShdr){
		startTmpObjList();
		buildSierpSubTri(8,scVal, name,0,depth,useShdr);
		endTmpObjList(1);			//bvh of sierp objs
		System.out.println("total buns : "+((PApplet.pow(4, depth)-1)/3.0f));
	}	
	//remove most recent object from list of objects and instead add to instance object struct.
	public void setObjectAsNamedObject(String name){
		myGeomBase _obj = allObjsToFind.remove(allObjsToFind.size()-1);
		--objCount;
		if(_obj instanceof myLight){												lightList.remove(lightList.size()-1);	numLights--;	} 
		else { 																		objList.remove(objList.size()-1); 		numNonLights--;	}	
		namedObjs.put(name, _obj);
		numInstances.put(name, 0);			//keep count of specific instances
		++numNamedObjs;
	}//setObjectAsNamedObject
	
	public void addInstance(String name, boolean addShdr){
		myGeomBase baseObj = namedObjs.get(name);
		myInstance _inst = new myInstance(this, baseObj);
		if(addShdr){		_inst.useInstShader();	}
		addObjectToScene(_inst, baseObj);			
		numInstances.put(name, numInstances.get(name)+1);
	}//
	
	// adds a new pointlight to the array of lights :   @params rgb - color, xyz - location
	public void addMyPointLight(String[] token){
		System.out.println("Point Light : current # of lights : " + numLights);
		myPointLight tmp = new myPointLight(this, numLights, 
				Double.parseDouble(token[4]),Double.parseDouble(token[5]),Double.parseDouble(token[6]),
				Double.parseDouble(token[1]),Double.parseDouble(token[2]),Double.parseDouble(token[3]));
		addObjectToScene(tmp);
	}//addMyPointLight method
	 
		// adds a new spotlight to the array of lights :   @params rgb - color, xyz - location, dx,dy,dz direction, innerThet, outerThet - angle bounds
	public void addMySpotLight(String[] token){
		double _inThet = Double.parseDouble(token[7]);
		double _outThet = Double.parseDouble(token[8]);
		System.out.println("Spotlight : current # of lights : " + numLights + " inner angle : " + _inThet + " outer angle : " + _outThet);
		mySpotLight tmp = new mySpotLight(this, numLights,
				Double.parseDouble(token[9]),Double.parseDouble(token[10]),Double.parseDouble(token[11]),
				Double.parseDouble(token[1]),Double.parseDouble(token[2]),Double.parseDouble(token[3]),
				Double.parseDouble(token[4]),Double.parseDouble(token[5]),Double.parseDouble(token[6]),
				_inThet,_outThet);  
		addObjectToScene(tmp);
	}//addMySpotLight method
	
	// adds a new disklight to the array of lights : disk_light x y z radius dx dy dz r g b
	public void addMyDiskLight(String[] token){
		double radius = Double.parseDouble(token[4]);
		System.out.println("Disk Light : current # of lights : " + numLights + " radius : " + radius);
		myDiskLight tmp = new myDiskLight(this, numLights,
				Double.parseDouble(token[8]),Double.parseDouble(token[9]),Double.parseDouble(token[10]),
				Double.parseDouble(token[1]),Double.parseDouble(token[2]),Double.parseDouble(token[3]),
				Double.parseDouble(token[5]),Double.parseDouble(token[6]),Double.parseDouble(token[7]),
				radius);  
		addObjectToScene(tmp);
	}//addMyDiskLight method	
	
	//read in prim data from txt file and create object
	public void readPrimData(String[] token){
		mySceneObject tmp = null;
		switch(token[0]){
		    case "box" : {//box xmin ymin zmin xmax ymax zmax :
		    	double minX = p.min(new double[]{Double.parseDouble(token[1]),Double.parseDouble(token[4])}),
		    		maxX = p.max(new double[]{Double.parseDouble(token[1]),Double.parseDouble(token[4])}),
		    		ctrX = (minX + maxX)*.5,
					minY = p.min(new double[]{Double.parseDouble(token[2]),Double.parseDouble(token[5])}),
		    		maxY = p.max(new double[]{Double.parseDouble(token[2]),Double.parseDouble(token[5])}),
		    		ctrY = (minY + maxY)*.5,
					minZ = p.min(new double[]{Double.parseDouble(token[3]),Double.parseDouble(token[6])}),
		    		maxZ = p.max(new double[]{Double.parseDouble(token[3]),Double.parseDouble(token[6])}),
		    		ctrZ = (minZ + maxZ)*.5;
		    	//putting box as a rendered bbox to minimize size of pure bboxes - rendered bbox is a bbox + shdr ref + some shdr-related functions and vars.
		    	tmp = new myRndrdBox(this,ctrX, ctrY, ctrZ, new myVector(minX, minY, minZ),	new myVector(maxX, maxY, maxZ));
		    	break;}			    
		    case "plane" : {			//infinite plane shape
		    	tmp = new myPlane(this);
		    	((myPlane)tmp).setPlaneVals(Double.parseDouble(token[1]), Double.parseDouble(token[2]),Double.parseDouble(token[3]),Double.parseDouble(token[4]));
		    	break;}
		    case "cyl" : {//old cylinder code
		    	double rad = Double.parseDouble(token[1]), hght = Double.parseDouble(token[2]);
		    	double xC = Double.parseDouble(token[3]), yC = Double.parseDouble(token[4]),zC = Double.parseDouble(token[5]);
		    	double xO = 0, yO = 1, zO = 0;
		    	try {
		    		xO = Double.parseDouble(token[6]);yO = Double.parseDouble(token[7]);zO = Double.parseDouble(token[8]);
		    	} catch (Exception e){	        			    	}
		    	tmp = new myCylinder(this, rad,hght,xC,yC,zC,xO,yO,zO);
		    	break;}			   
		    case "cylinder" : { //new reqs : cylinder radius x z ymin ymax
		    	double rad = Double.parseDouble(token[1]), xC = Double.parseDouble(token[2]), zC = Double.parseDouble(token[3]);
		    	double yMin  = Double.parseDouble(token[4]), yMax = Double.parseDouble(token[5]);
		    	double hght = yMax - yMin;
		    	
		    	double xO = 0,yO = 1,zO = 0;
		    	tmp = new myCylinder(this, rad,hght,xC,yMin,zC,xO,yO,zO);
		    	break;}			    

		    case "hollow_cylinder" : {//hollow_cylinder radius x z ymin ymax
		    	double rad = Double.parseDouble(token[1]), xC = Double.parseDouble(token[2]), zC = Double.parseDouble(token[3]);
		    	double yMin  = Double.parseDouble(token[4]), yMax = Double.parseDouble(token[5]);
		    	double hght = yMax - yMin;
		    	
		    	double xO = 0,yO = 1,zO = 0;
		    	tmp = new myHollow_Cylinder(this, rad,hght,xC,yMin,zC,xO,yO,zO);		    	
		    	break;}
		    case "sphere" : {
		    	//create sphere
		    	tmp = new mySphere(this, Double.parseDouble(token[1]),Double.parseDouble(token[2]),Double.parseDouble(token[3]),Double.parseDouble(token[4]),true);
		    	break;}			    
		    case "moving_sphere": {//moving_sphere radius x1 y1 z1 x2 y2 z2
		    	tmp = new myMovingSphere(this, 
					Double.parseDouble(token[1]),Double.parseDouble(token[2]),Double.parseDouble(token[3]),Double.parseDouble(token[4]), 
					Double.parseDouble(token[5]),Double.parseDouble(token[6]),Double.parseDouble(token[7]), true);
		    	break;}			    
		    case "sphereIn" : {
		    	//create sphere with internal reflections - normals point in
		    	tmp = new mySphere(this, Double.parseDouble(token[1]),Double.parseDouble(token[2]),Double.parseDouble(token[3]),Double.parseDouble(token[4]),true);
		    	tmp.rFlags[mySceneObject.invertedIDX] = true;
		    	break;}	
		    case "ellipsoid" : {//create elliptical sphere with 3 radii elements in each of 3 card directions			    	
		    	tmp = new mySphere(this, 
		    			Double.parseDouble(token[1]),Double.parseDouble(token[2]),Double.parseDouble(token[3]),
		    			Double.parseDouble(token[4]),Double.parseDouble(token[5]),Double.parseDouble(token[6]));
		    	break;}
		    default :{
		    	System.out.println("Object type not handled : "+ token[0]);
		    	return;
		    }
		}//switch
		//set shader and texture
		tmp.shdr = getCurShader();
		//add object to scene
    	addObjectToScene(tmp);		
	}//
	
	//return a shader built with the current settings
	public myObjShader getCurShader(){
		myObjShader tmp = (scFlags[simpleRefrIDX]) ? new mySimpleReflObjShdr(this) : new myObjShader(this);
		tmp.txtr = getCurTexture(tmp);				
		return tmp;
	}//getCurShader
	
	//return appropriate texture handler
	public myTextureHandler getCurTexture(myObjShader tmp){
		switch (txtrType){
			case 0 : {	return new myNonTexture(this,tmp);}						//diffuse/shiny only
			case 1 : { 	return new myImageTexture(this,tmp);}					//has an image texture
			case 2 : {	return new myNoiseTexture(this,tmp,noiseScale);}		//perlin noise
			case 3 : {	return new myBaseWoodTexture(this,tmp,noiseScale);}
			case 4 : { 	return new myMarbleTexture(this,tmp,noiseScale);}
			case 5 : { 	return new myCellularTexture(this,tmp,noiseScale);}	
			case 6 : { 	return new myWoodTexture(this,tmp,noiseScale);}
			default : { return new myNonTexture(this,tmp);}
		}
	}//myTextureHandler
	//return appropriate texture handler
	public String getTxtrName(){
		switch (txtrType){
			case 0 : {	return "No Texture";}						//diffuse/shiny only
			case 1 : { 	return "Image Texture";}					//has an image texture
			case 2 : {	return "Pure Noise Texture";}		//perlin noise
			case 3 : {	return "Base Wood Texture";}
			case 4 : { 	return "Marble Texture";}
			case 5 : { 	return "Cellular Texture";}	
			case 6 : { 	return "Alt Wood Texture";}
			default : { return "No Texture";}
		}
	}//myTextureHandler
	
	//entry point
	public void addObjectToScene(myGeomBase _obj){addObjectToScene(_obj,_obj);}
	public void addObjectToScene(myGeomBase _obj, myGeomBase _cmpObj){
		if(scFlags[addToTmpListIDX]){tmpObjList.add(_obj); return;}
		if(_cmpObj instanceof myLight){			lightList.add(_obj);	numLights++;} 
		else {									objList.add(_obj);		numNonLights++;}
		allObjsToFind.add(_obj);
		objCount++;
	}//addObjectToScene

	/////////
	//setting values from reader
	////
	
	public void setNoise(double scale, String[] vals){
		resetDfltTxtrVals();//reset values for next call
		txtrType = 2;
		System.out.println("Setting Noise to scale : " + scale);	
		noiseScale = scale;		//only specified value is scale currently
	}//setNoise
	
	//call after a proc texture is built, to reset values to defaults
	public void resetDfltTxtrVals(){
		setProcTxtrVals(new int[]{0, 4, 1, 2, 1, 1}, 
				new boolean[]{false,false,false}, 
				new double[]{1.0, 1.0, 5.0, .1, 1.0, 0.05}, 
				new myVector[]{new myVector(1.0,1.0,1.0)}, 
				new myColor[]{ p.getClr("clr_nearblack"),p.getClr("clr_white")},
				new Double[]{.5,.5});
	}//resetDfltTxtrVals

	//int[] ints =  txtrType,  numOctaves,  numOverlays,  numPtsDist,distFunc, roiFunc
	//boolean[] bools = rndColors , useCustClrs, useFwdTrans
	//double[] dbls = noiseScale, turbMult , colorScale , colorMult, avgNumPerCell, mortarThresh
	//myVector[] vecs = pdMult
	//myColor[] clrs = noiseColors
	//Double[] wts = clrWts
	private void setProcTxtrVals(int[] ints, boolean[] bools, double[] dbls, myVector[] vecs, myColor[] clrs, Double[] wts){
		txtrType = ints[0];	numOctaves = ints[1];numOverlays = ints[2];	numPtsDist = ints[3];	distFunc = ints[4];	roiFunc = ints[5];	
		rndColors = bools[0];useCustClrs = bools[1];useFwdTrans = bools[2];		
		noiseScale = dbls[0];turbMult = dbls[1];colorScale = dbls[2];colorMult = dbls[3];avgNumPerCell = dbls[4];	mortarThresh = dbls[5];		
		pdMult = vecs[0];		
		noiseColors = clrs;
//		clrWts = wts;		
	}//setProcTxtrVals

	//set colors used by proc texture
	public void setTxtrColor(String[] clrs){
		//get current noise color array
		if(!useCustClrs){
			noiseColors = new myColor[0];
			clrWts = new Double[0];
			useCustClrs = true;
		}
		ArrayList<myColor> tmpAra = new ArrayList<myColor>(Arrays.asList(noiseColors));
		ArrayList<Double> tmpWtAra = new ArrayList<Double>(Arrays.asList(clrWts));
		//<noise color spec tag> (<'named'> <clr name>) or  (<color r g b>)  <wt> <-specify once for each color
		try{	
			myColor tmp = null;
			int wtIdx;
			//name has format "clr_<colorname>"
			if(clrs[1].equals("named")){	tmp = p.getClr(clrs[2]);	wtIdx = 3;} 
			else {							tmp = p.readColor(clrs, 1);	wtIdx = 4;}
			tmpAra.add(tmp);
			try{
				double wt = Double.parseDouble(clrs[wtIdx]);
				tmpWtAra.add(wt);	//normalize at end, when all colors have been added and txtr being built		
			}
			catch(Exception e) {//no weight specified, just add avg of existing weights
				double tmpWt = 0;
				if(tmpWtAra.size() == 0){tmpWt = 1.0;} 
				else {
					int cnt = 0;
					for(Double wt : tmpWtAra){	tmpWt += wt;cnt++;}
					tmpWt /= (1.0*cnt);
				}
				tmpWtAra.add(tmpWt);
			}//catch		
			System.out.println("Finished loading color : " + tmp + " for txtr " + getTxtrName());
		}
		catch (Exception e) {String res = "Invalid color specification : " ;	for(int i =0; i<clrs.length;++i){res+=" {"+clrs[i]+"} ";}res+=" so color not added to array";System.out.println(res);}	 		
		noiseColors = tmpAra.toArray(new myColor[0]);
	}//setTxtrColors
	
	//read in constants configured for perlin noise
	private boolean readProcTxtrPerlinVals(String[] vals){
		boolean useDefaults;
		//may just have <typ> or may have up to color scale, or may have all values - use defaults for all values not specified
		//<typ> <noise scale> <numOctaves> <turbMult> <pdMult x y z> <multByPI 1/0 1/0 1/0> <useFwdTransform 0/1> <?rndomize colors colorScale - if present then true> <color mult> <num overlays - if present, otherwise 1>
		try{
			noiseScale = Double.parseDouble(vals[1]);
			numOctaves = Integer.parseInt(vals[2]);
			turbMult = Double.parseDouble(vals[3]);
			pdMult = new myVector(Double.parseDouble(vals[4]),Double.parseDouble(vals[5]),Double.parseDouble(vals[6]));
			myVector pyMult = new myVector(Double.parseDouble(vals[7]),Double.parseDouble(vals[8]),Double.parseDouble(vals[9]));
			if(pyMult.sqMagn > 0){									//if vector values specified
				pyMult._mult(PConstants.TWO_PI-1.0);				//change 0|1 to 0|2pi-1 vector
				pyMult._add(1.0,1.0,1.0);							//change 0|2pi-1 vector to 1|2pi vector
				pdMult = p._elemMult(pdMult, pyMult);
			}
			useFwdTrans = (Double.parseDouble(vals[10]) == 1.0);		//whether or not to use fwd transform on points
			try{
				colorScale = Double.parseDouble(vals[11]);
				colorMult = Double.parseDouble(vals[12]);
				rndColors = true;
				try{	numOverlays = Integer.parseInt(vals[13]);	}	catch (Exception e) {numOverlays = 1;	}		
			}
			catch (Exception e) {	
				System.out.println("Proc Perlin-based txtr not specifying randomize colors or overlays so defaults are used for txtr type : " + getTxtrName());	
				rndColors = false;	colorScale = 25.0;colorMult = .1;numOverlays = 1;}	 
			useDefaults = false;
			System.out.println("Finished loading custom values for texture type : " + getTxtrName());
		}
		catch (Exception e) {System.out.println("No Proc Texture values specified for texture type : " + getTxtrName() + " so using defaults.");	useDefaults = true;	}	 
		return useDefaults;	
	}//readProcTxtrPerlinVals
	
	//read in constants configured for worley noise
	private boolean readProcTxtrWorleyVals(String[] vals){
		boolean useDefaults;
		//different format than perlin
		//may just have <typ> or may have up to color scale, or may have all values - use defaults for all values not specified
		//<typ> <noise scale> <distfunction 0=man/1=euc> <roiFunc 0=altLinSum/1=?><num pts for dist func - should be even> <avg # pts per cell> <mortar threshold 0.0-1.0> <useFwdTransform 0/1> <?rndomize colors colorScale - if present then true> <color mult> <num overlays - if present, otherwise 1>
		try{
			int parseIDX = 1;
			noiseScale = Double.parseDouble(vals[parseIDX++]);
			distFunc = Integer.parseInt(vals[parseIDX++]);
			roiFunc = Integer.parseInt(vals[parseIDX++]);	//function for determining region of interest in worley noise - 0 : nearest lin sum 1 : alternating linear sum, 2+: ....
			numPtsDist = Integer.parseInt(vals[parseIDX++]);
			avgNumPerCell = Double.parseDouble(vals[parseIDX++]);
			mortarThresh = Double.parseDouble(vals[parseIDX++]);
			useFwdTrans = (Double.parseDouble(vals[parseIDX++]) == 1.0);		//whether or not to use fwd transform on points
			try{
				colorScale = Double.parseDouble(vals[parseIDX++]);
				colorMult = Double.parseDouble(vals[parseIDX++]);
				rndColors = true;
				try{	numOverlays = Integer.parseInt(vals[parseIDX++]);	}	catch (Exception e) {numOverlays = 1;	}		
			}
			catch (Exception e) {	
				System.out.println("Proc txtr not specifying randomize colors or overlays so defaults are used for txtr type : " + getTxtrName());	
				rndColors = false;	colorScale = 25.0;colorMult = .1;numOverlays = 1;}	 
			useDefaults = false;
			System.out.println("Finished loading custom values for texture type : " + getTxtrName());
		}
		catch (Exception e) {System.out.println("No Proc Texture values specified for texture type : " + getTxtrName() + " so using defaults.");	useDefaults = true;	}	 
		return useDefaults;	
	}//readProcTxtrWorleyVals
	
	//read in procedural texture values for perlin noise and populate globals used to build txtr
	public boolean readProcTxtrVals(String[] vals, boolean isPerlin){
		if(isPerlin){return readProcTxtrPerlinVals(vals);}
		else { return readProcTxtrWorleyVals(vals);}
	}//readProcTxtrVals	
	
	//proc texture components
	public void setTexture(String[] vals){
		//TODO get all these values from CLI
		resetDfltTxtrVals();//reset values for next call
		String _typ = vals[0];
		switch(_typ){
			case "wood":{
				txtrType = 3;
				boolean useDefaults = readProcTxtrVals(vals, true);
				//may be overwritten by later color commands in cli
				if(!useCustClrs){noiseColors = new myColor[]{ p.getClr("clr_dkwood1"),	p.getClr("clr_ltwood1")};clrWts = new Double[] {1.0,1.0};}
				if(useDefaults){					
					setProcTxtrVals(new int[]{txtrType, 4, 1, 2, 1,1}, 				//int[] ints =  txtrType,  numOctaves,  numOverlays,  numPtsDist,distFunc (not used for perlin), roiFunc (not used for perlin)
							new boolean[]{true,useCustClrs,false}, 					//boolean[] bools = rndColors , useCustClrs, useFwdTrans -> useFwdTrans needs to be specified in command in cli
							new double[]{2.0, .4, 25.0, .2, 1.0, 0.05}, 			//double[] dbls = noiseScale, turbMult , colorScale , colorMult, avgNumPerCell,mortarThresh
							new myVector[]{new myVector(PConstants.TWO_PI* 2.7,3.6,4.3)}, 	//myVector[] vecs = pdMult
							noiseColors, 											//myColor[] clrs = noiseColors
							clrWts);												//Double[] wts = clrWts					
				} break;}		
			case "wood2"  : {//yellow-ish by default
				txtrType = 6;
				boolean useDefaults = readProcTxtrVals(vals, true);
				if(!useCustClrs){noiseColors = new myColor[]{ p.getClr("clr_dkwood2"),	p.getClr("clr_ltwood2")};}// 		clrWts = new Double[] {1.0,1.0};}
				//turbulence values
				if(useDefaults){
					setProcTxtrVals(new int[]{txtrType, 8, 1, 2, 1, 1}, 			//int[] ints =  txtrType,  numOctaves,  numOverlays,  numPtsDist,distFunc (not used), roiFunc (not used)
							new boolean[]{true,useCustClrs,false}, 					//boolean[] bools = rndColors , useCustClrs, useFwdTrans -> useFwdTrans needs to be specified in command in cli
							new double[]{1.0, .4, 25.0, .3, 1.0, 0.05}, 			//double[] dbls = noiseScale, turbMult , colorScale , colorMult, avgNumPerCell
							new myVector[]{new myVector(PConstants.TWO_PI*3.5,7.9,6.2)}, 	//myVector[] vecs = pdMult
							noiseColors, 											//myColor[] clrs = noiseColors
							clrWts);												//Double[] wts = clrWts					
				} break;}		
			case "marble":{
				txtrType = 4;
				boolean useDefaults = readProcTxtrVals(vals, true);
				if(!useCustClrs){noiseColors = new myColor[]{ p.getClr("clr_nearblack"),	p.getClr("clr_offwhite")}; clrWts = new Double[] {1.0,1.0};}
				//turbulence values
				if(useDefaults){
				setProcTxtrVals(new int[]{txtrType, 16, 1, 2, 1, 1}, 				//int[] ints =  txtrType,  numOctaves,  numOverlays,  numPtsDist,distFunc (not used), roiFunc (not used)
							new boolean[]{true,useCustClrs,false}, 					//boolean[] bools = rndColors , useCustClrs, useFwdTrans -> useFwdTrans needs to be specified in command in cli
							new double[]{1.0, 15.0, 24.0, .1, 1.0, 0.05}, 			//double[] dbls = noiseScale, turbMult , colorScale , colorMult, avgNumPerCell
							new myVector[]{new myVector(PConstants.TWO_PI * 0.1,PConstants.TWO_PI * 31.4,PConstants.TWO_PI *4.1)}, 	//myVector[] vecs = pdMult
							noiseColors, 											//myColor[] clrs = noiseColors
							clrWts);												//Double[] wts = clrWts					
				} break;}			
			//this uses 
			case "stone":{
				txtrType = 5;
				boolean useDefaults = readProcTxtrVals(vals, false);
				if(!useCustClrs){
					noiseColors = new myColor[]{p.getClr("clr_mortar1"), p.getClr("clr_mortar2"),
							p.getClr("clr_brick1_1"),p.getClr("clr_brick1_2"),p.getClr("clr_brick2_1"),p.getClr("clr_brick2_2"),				//"brick2" color 1,2
							p.getClr("clr_brick3_1"),p.getClr("clr_brick3_2"),p.getClr("clr_brick4_1"),p.getClr("clr_brick4_2")};
					clrWts = new Double[] {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};//normalized in txtr
				}
				if(useDefaults){
					setProcTxtrVals(new int[]{txtrType, 8, 1, 2, 1, 1}, 		//int[] ints =  txtrType,  numOctaves,  numOverlays,  numPtsDist, distfunc
							new boolean[]{true,useCustClrs,false}, 				//boolean[] bools = rndColors , useCustClrs, useFwdTrans -> useFwdTrans needs to be specified in command in cli
							new double[]{4.0, 1.0, 12.0, .2, 1.0, 0.05},		//double[] dbls = noiseScale, turbMult , colorScale , colorMult, avgNumPerCell, mortar thresh
							new myVector[]{new myVector(10.0,10.0,10.0)}, 		//myVector[] vecs = pdMult - not used here yet
							noiseColors, 											//myColor[] clrs = noiseColors
							clrWts);											//Double[] wts = clrWts					
				} break;}
			default : {	System.out.println("Unknown Texture type : " + _typ); txtrType = 0; return;}
		}
		System.out.println("Set Texture type : " + _typ); 
	}//setTexture
	
	
	public void setImageSize(int numCols, int numRows){//set size and all size-related variables, including image dims
		sceneCols = numCols;
		sceneRows = numRows;
		numPxls = sceneRows * sceneCols;
		rayYOffset = sceneRows/2.0;
		rayXOffset = sceneCols/2.0;
		
		maxDim = Math.max(sceneRows,sceneCols);
		yStart = ((maxDim - sceneRows)/2.0) - rayYOffset;	
		xStart = ((maxDim - sceneCols)/2.0) - rayXOffset;			//compensate for # rows or # cols not being max - make sure projection is centered in non-square images
		fishMult = 2.0/maxDim; 
			
		rndrdImg = p.createImage(sceneCols,sceneRows,PConstants.RGB);		
	}
	
	//refining
	public void setRefine(String refState){
    	curRefineStep = 0;
		scFlags[myScene.glblRefineIDX] = refState.toLowerCase().equals("on");
		//build refinement #pxls array dynamically by finding average dim of image and then math.
		int refIDX = (int)(Math.log10(.5*(this.sceneCols + this.sceneRows)/16.0)/p.log10_2);
		RefineIDX = new int[(refIDX+1)];
		for(int i =refIDX; i >=0; --i){	RefineIDX[refIDX-i]=p.pow2[i];}
	}//setRefine
		
	public void setDpthOfFld(double lRad, double lFD){//depth of field effect
		lens_radius = lRad;
		lens_focal_distance = lFD;
		scFlags[hasDpthOfFldIDX] = true;
		focalPlane.setPlaneVals(0, 0, 1, lens_focal_distance);      //needs to be modified so that d = viewZ + lens_focal_distance once FOV has been specified.
	}
	
	public void setNumRaysPerPxl(int _num){
		System.out.println("Num Rays Per Pixel : " + _num);
		this.numRaysPerPixel = _num;
	}
	
	public void setSurface(myColor Cdiff, myColor Camb, myColor Cspec, double phongExp, double newKRefl){
		txtrType = 0;			//set to be no texture
		currDiffuseColor.set(Cdiff);
		currAmbientColor.set(Camb);
		currSpecularColor.set(Cspec);
		currPhongExp = phongExp;
		setKRefl(newKRefl,newKRefl,newKRefl,newKRefl);
		setRfrIdx(0, 0, 0, 0);	
		currKTrans = 0;
		globRfrIdx = 0;
		globCurPermClr.set(0,0,0);
	}//setSurface method

	public void setSurface(myColor Cdiff, myColor Camb, myColor Cspec, double phongExp, double KRefl, double KTrans){
		setSurface(Cdiff,Camb, Cspec,phongExp, KRefl);
		currKTrans = KTrans;
		setRfrIdx(0, 0, 0, 0);	
	}//setSurface method with refractance

	public void setSurface(myColor Cdiff, myColor Camb, myColor Cspec, double phongExp, double KRefl, double KTrans, double rfrIdx){
		setSurface(Cdiff,Camb, Cspec,phongExp, KRefl, KTrans);
		//set permiability of object to light
		setRfrIdx(rfrIdx, rfrIdx, rfrIdx,rfrIdx);
	}//setSurface method with refractance
	public void setDepth(double depth){  currDepth = depth;}//subsurface depth - TODO
	
	public void setPhong(double newPhong){ currPhongExp = newPhong;}
	public void setKTrans(double newKTrans){  currKTrans = newKTrans;}
	public void setRfrIdx(double rfrIdx){setRfrIdx(rfrIdx, rfrIdx, rfrIdx, rfrIdx);}
	public void setRfrIdx(double rfrIdx, double pR, double pG, double pB){
		globRfrIdx = rfrIdx;
		//can control color "tint" of transmitted ray through object.  for now all the same
		globCurPermClr.set(pR,pG,pB);
	}
	
	public void setKRefl(double newKRefl){setKRefl(newKRefl,newKRefl,newKRefl,newKRefl);}
	public void setKRefl(double newKRefl, double kr, double kg, double kb){
	  currKRefl  = newKRefl;
	  //can control color "tint" of transmitted ray through object		
	  currKReflClr.set(kr,kg,kb);
	}
	public void setBackgroundColor(double r, double g, double b){  backgroundColor.set(r,g,b);}//initBackground 3 color method
//	public void setForegroundColor(double r, double g, double b){  foregroundColor.set(r,g,b);  scFlags[useFGColorIDX] = true;}//initBackground 3 color method
	////////
	///end setting values
	////////
	
	///////
	//RT functionality
	///////	
	//get random location within "lens" for depth of field calculation - consider loc to be center, pick random point in constant z plane within some radius of loc point
	public myVector getDpthOfFldEyeLoc(myVector loc){
		myVector tmp = p.rotVecAroundAxis(new myVector(0,1,0),new myVector(0,0,-1),ThreadLocalRandom.current().nextDouble(0,PConstants.TWO_PI));				//rotate surfTangent by random angle
		tmp._normalize();
		double mult = ThreadLocalRandom.current().nextDouble(0,lens_radius);			//find displacement radius from origin
		tmp._mult(mult);
		tmp._add(loc);																														//find displacement point on origin
		return tmp;		
	}
	//determines if a light source is blocked by another object for shadow detection
	//currently only returns 1 or 0 if light is blocked
	
	public int calcShadow(myRay _ray, double distToLight){
		//for each object in scene, check if intersecting any objects before hitting light
		for (myGeomBase obj : objList){
			if(obj.calcShadowHit(_ray, _ray.getTransformedRay(_ray, obj.CTMara[obj.invIDX]), obj.CTMara, distToLight) == 1){	return 1;}
		}//for each object in scene
		return 0;
	}//findLight method
	
	//eventually multithread/shdr per object?
	public rayHit findClosestRayHit(myRay _ray){
		//objList does not hold lights - no need to check pointlights - TODO need to check lights for non-point lights- ?	
		TreeMap<rayHit, myGeomBase>objsAtRayHits = new TreeMap<rayHit,myGeomBase>();
		objsAtRayHits.put(new rayHit(false), null);
		//myRay transRay;
		for (myGeomBase obj : objList){	
			rayHit _hit = null;
			try{
				_hit = obj.intersectCheck(_ray,_ray.getTransformedRay(_ray, obj.CTMara[obj.invIDX]),obj.CTMara);		
			} catch (Exception e){
				System.out.println("find closest ray hit exception :"+e);
			}
			if(_hit.isHit){			objsAtRayHits.put(_hit, _hit.obj);		}
		}//for obj in scenelist
		return objsAtRayHits.firstKey();
	}//findClosestRayHit
	
	
	//determine color of a reflected ray - careful with recursive depth  
	public myColor reflectRay(myRay _ray){
		rayHit hitChk = findClosestRayHit(_ray);
		//if ((hitChk.isHit)) {												return(hitChk.obj.getColorAtPos(hitChk));}//to debug BVH use this - displays colors of leaf boxes (red/blue)
		if ((hitChk.isHit)) {												return(hitChk.shdr.getColorAtPos(hitChk));}
		else if (scFlags[glblTxtrdBkgIDX]) {								return getBackgroundTextureColor(_ray);	} 	//using skydome
//		else if ((_ray.direction.z > epsVal) && (scFlags[useFGColorIDX])){	return foregroundColor;	} 					//for getting color reflected from behind viewer
		else {																return backgroundColor;	}
	}//reflectRay	
	
	
	///////////////
	//setup photon handling kdtree
	public void setPhotonHandling(String[] token){
		scFlags[usePhotonMapIDX] = true;
		scFlags[isPhtnMapRndrdIDX] = false;
		String type = token[0];
		scFlags[isCausticPhtnIDX] = type.contains("caustic");
		//type is "caustic_photons" or "diffuse_photons"
		//caustic_photons 80000000  80 0.05
		numPhotons = Integer.parseInt(token[1]);
		kNhood = Integer.parseInt(token[2]);
		ph_max_near_dist = Float.parseFloat(token[3]);
		//build photon tree
		photonTree = new myKD_Tree(this, numPhotons, kNhood, ph_max_near_dist);
	}
	//final gather procedure
	public void setFinalGather(String[] token){
		numGatherRays = Integer.parseInt(token[1]);
	}
	/////////////

	//from a hit on an object, return reflection dir
  	public myVector getReflDir(rayHit hit){
  		myVector backToEyeDir = new myVector(hit.fwdTransRayDir);  		
  		backToEyeDir._mult(-1);  		
 		double dotProd = 2 * (backToEyeDir._dot(hit.objNorm));
  		myVector tempVect = new myVector(hit.objNorm.x * dotProd,hit.objNorm.y * dotProd,hit.objNorm.z * dotProd);
  		myVector reflDir = p._sub(tempVect,backToEyeDir);
  		reflDir._normalize();
  		return reflDir;  
  	}//getReflDir
  				
  	//TODO set  up seperate maps for caustic and diffuse photons
	//send out photons for all lights if we are using photon mapping - put these photons in kd tree
	//need different photon configurations/processes for caustic photons and indirect illumination (diffuse)photons
	protected void sendCausticPhotons(){
		int numDiv = 100,lastCastCnt = 0, starCount = 0, pctCastCnt = numDiv/10;
		int numCastPerDisp = photonTree.num_Cast/numDiv;
		myLight tmpLight; 
		double pwrMult = causticsLightPwrMult/photonTree.num_Cast;
		myRay reflRefrRay;
		double[] tmpPwr,photon_pwr, d;
		rayHit hitChk;
		myPhoton phn;
		//TODO scale # of photons sent into scene by light intensity
		for(myGeomBase light : lightList){//either a light or an instance of a light
			tmpLight = (light instanceof myLight) ? tmpLight = (myLight)light : ((myLight)((myInstance)light).obj);
			System.out.print("Casting " + photonTree.num_Cast + " Caustic photons for light ID " + tmpLight.ID + ": Progress:");			
			for(int i =0; i<photonTree.num_Cast; ++i){
				photon_pwr = new double[]{tmpLight.lightColor.RGB.x * pwrMult,tmpLight.lightColor.RGB.y * pwrMult,tmpLight.lightColor.RGB.z * pwrMult };
				hitChk = findClosestRayHit(tmpLight.genRndPhtnRay());//first hit
				if((!hitChk.isHit) || (!hitChk.shdr.shdrFlags[hitChk.shdr.hasCaustic])){continue;}			//either hit background or 1st hit is diffuse object - caustic has to hit spec first
				//System.out.println("hit obj : " + hitChk.obj.ID);
				//we have first hit here at caustic-generating surface.  need to propagate through to first diffuse surface
				hitChk.phtnPwr = photon_pwr;
				do{
		 			reflRefrRay = hitChk.shdr.findCausticRayHit(hitChk,hitChk.phtnPwr);
		 			if(reflRefrRay != null){	
		 				tmpPwr = new double[]{hitChk.phtnPwr[0],hitChk.phtnPwr[1],hitChk.phtnPwr[2]};
		 				hitChk = findClosestRayHit(reflRefrRay);
		 				hitChk.phtnPwr = tmpPwr;
		 			} else {					hitChk.isHit = false;}
				} while((hitChk.isHit) && (hitChk.shdr.shdrFlags[hitChk.shdr.hasCaustic]) && (reflRefrRay.gen <= numPhotonRays));				//keep going while we have a hit and we are hitting a caustic
				if((!hitChk.isHit) || (reflRefrRay.gen > numPhotonRays)){continue;}																//bounced off into space
				
				//d = hitChk.fwdTransRayDir._normalized().getAsAra();				
				//phn = new myPhoton(photonTree, hitChk.phtnPwr, hitChk.fwdTransHitLoc, Math.acos(d[2]), PConstants.PI + Math.atan2(d[1], d[0])); 	
				//phn = new myPhoton(photonTree, photon_pwr, hitChk.fwdTransHitLoc.x,  hitChk.fwdTransHitLoc.y,  hitChk.fwdTransHitLoc.z); 	
				phn = new myPhoton(photonTree, hitChk.phtnPwr, hitChk.fwdTransHitLoc.x,  hitChk.fwdTransHitLoc.y,  hitChk.fwdTransHitLoc.z); 	
				photonTree.add_photon(phn);
				//this just calcs when to display progress bar, can be deleted
				if(i > lastCastCnt){
					lastCastCnt += numCastPerDisp;					
					System.out.print((starCount % pctCastCnt == 0) ? (10.0*starCount/pctCastCnt) + "%" : "*");
					starCount++;
				}
			}//for each photon of light
			System.out.println("100.0%");
			starCount = 0;lastCastCnt=0;
		}//for each light
		photonTree.build_tree();
	}//sendCausticPhotons
	
	protected void sendDiffusePhotons(){
		//for every light
		//for each photon of n, 
		//for every object
		//check if hit, save where lands
		int numDiv = 100,lastCastCnt = 0, starCount = 0, pctCastCnt = numDiv/10;
		int numCastPerDisp = photonTree.num_Cast/numDiv;
		myLight tmpLight; 
		double pwrMult = diffuseLightPwrMult/photonTree.num_Cast;
		myRay reflRefrRay;
		double[] tmpPwr,photon_pwr, d;
		rayHit hitChk;
		myPhoton phn;
		//TODO scale # of photons sent into scene by light intensity
		for(myGeomBase light : lightList){//either a light or an instance of a light
			tmpLight = (light instanceof myLight) ? tmpLight = (myLight)light : ((myLight)((myInstance)light).obj);
			System.out.print("Casting " + photonTree.num_Cast + " Diffuse (indirect) photons for light ID " + tmpLight.ID + ": Progress:");			
			for(int i =0; i<photonTree.num_Cast; ++i){
				photon_pwr = new double[]{tmpLight.lightColor.RGB.x * pwrMult,tmpLight.lightColor.RGB.y * pwrMult,tmpLight.lightColor.RGB.z * pwrMult };
				hitChk = findClosestRayHit(tmpLight.genRndPhtnRay());//first hit
				if(!hitChk.isHit){continue;}							//hit background - ignore
				//now we hit an object, spec or diffuse - if specular, bounce without storing, if diffuse store and bounce with prob based on avg color				
				//System.out.println("hit obj : " + hitChk.obj.ID);
				//we have first hit here at caustic-generating surface.  need to propagate through to first diffuse surface
				hitChk.phtnPwr = photon_pwr;
				boolean done = false, firstDiff = true;
				do{
					if(hitChk.shdr.KRefl == 0){//diffuse, store and maybe bounce
						double prob = 0;
						if(!firstDiff){//don't store first
							phn = new myPhoton(photonTree, hitChk.phtnPwr, hitChk.fwdTransHitLoc.x,  hitChk.fwdTransHitLoc.y,  hitChk.fwdTransHitLoc.z); 	
							photonTree.add_photon(phn);
							prob = ThreadLocalRandom.current().nextDouble(0,1.0);//russian roulette to see if casting
						}
						firstDiff = false;
						if(prob < hitChk.shdr.avgDiffClr){	//reflect in new random dir, scale phtn power by diffClr/avgClr
							//get new bounce dir
							myVector hitLoc = hitChk.fwdTransHitLoc;
					  		//first calc random x,y,z
					  		double x=0,y=0,z=0, sqmag;
							do{
								x = ThreadLocalRandom.current().nextDouble(-1.0,1.0);
								y = ThreadLocalRandom.current().nextDouble(-1.0,1.0);			
								sqmag = (x*x) + (y*y);
							}
							while ((sqmag >= 1.0) || (sqmag < p.epsVal));
							z = Math.sqrt(1 - (sqmag));							//cosine weighting preserved by projecting up to sphere
					  		//then build ortho basis from normal - n' , q' , r' 
					  		myVector n = new myVector(hitChk.objNorm),_p = new myVector(),_q = new myVector(); 
					  				//tmpV = (((n.x > n.y) && (n.x > n.z)) || ((-n.x > -n.y) && (-n.x > -n.z))  ? new myVector(0,0,1)  : new myVector(1,0,0));//find vector not close to n or -n to use to find tangent
							double nxSq = n.x * n.x, nySq = n.y * n.y, nzSq = n.z * n.z;
							myVector tmpV = (((nxSq > nySq) && (nxSq > nzSq))  ? new myVector(0,0,1)  : new myVector(1,0,0));//find vector not close to n or -n to use to find tangent
					  		//set _p to be tangent, _q to be binorm
					  		_p = n._cross(tmpV);	_q = _p._cross(n);
					  		//if(_p.sqMagn < p.epsVal){System.out.println("bad _p : " + _p + " | n : " + n + " | tmpV : " + tmpV);}
					  		//lastly multiply ortho basis vectors by x,y,z : x * p, y * q', z*n', and then sum these products - z is projection/hemisphere dir, so should coincide with normal
					  		n._mult(z);	_p._mult(x);_q._mult(y);
					  		myVector bounceDir = new myVector(n.x + _p.x + _q.x,n.y + _p.y + _q.y,n.z + _p.z + _q.z);
					  		bounceDir._normalize();
					 		//save power before finding ray hit, to reset it after ray hit
					  		tmpPwr = new double[]{hitChk.phtnPwr[0]*hitChk.shdr.phtnDiffScl.x,hitChk.phtnPwr[1]*hitChk.shdr.phtnDiffScl.y,hitChk.phtnPwr[2]*hitChk.shdr.phtnDiffScl.z};			//hitChk changes below, we want to propagate tmpPwr
					 		//new photon ray - photon power : 
			 				reflRefrRay = new myRay(this, hitLoc, bounceDir, hitChk.transRay.gen+1);
			 				hitChk = findClosestRayHit(reflRefrRay);
			 				hitChk.phtnPwr = tmpPwr;
						} else {	done = true;}						
					} else {					//specular, just bounce, scale power by avg spec power (?)
			 			reflRefrRay = hitChk.shdr.findCausticRayHit(hitChk,hitChk.phtnPwr);
			 			if(reflRefrRay != null){	
			 				tmpPwr = new double[]{hitChk.phtnPwr[0],hitChk.phtnPwr[1],hitChk.phtnPwr[2]};//save updated photon power
			 				hitChk = findClosestRayHit(reflRefrRay);
			 				hitChk.phtnPwr = tmpPwr;
			 			} else {					hitChk.isHit = false;}
					}
				} while((hitChk.isHit) && (!done) && (hitChk.transRay.gen <= numPhotonRays));				//keep going while we have a hit and we are less than ray recursion depth
				if((!hitChk.isHit) || (hitChk.transRay.gen > numPhotonRays)){continue;}		//bounced off into space
				
				//d = hitChk.fwdTransRayDir._normalized().getAsAra();			//direction of ray hit, for anisotropic materials (TODO)			
				//phn = new myPhoton(photonTree, hitChk.phtnPwr, hitChk.fwdTransHitLoc, Math.acos(d[2]), PConstants.PI + Math.atan2(d[1], d[0])); 	
				//phn = new myPhoton(photonTree, photon_pwr, hitChk.fwdTransHitLoc.x,  hitChk.fwdTransHitLoc.y,  hitChk.fwdTransHitLoc.z); 	
				//this just calcs when to display progress bar, can be deleted
				if(i > lastCastCnt){
					lastCastCnt += numCastPerDisp;					
					System.out.print((starCount % pctCastCnt == 0) ? (10.0*starCount/pctCastCnt) + "%" : "*");
					starCount++;
				}
			}
			System.out.println("100.0%");
			starCount = 0;lastCastCnt=0;
		}
		photonTree.build_tree();
	}//sendDiffusePhotons

	
	
	//initialize drawing routine - build photon map if it exists
	protected void initRender(){
		rndrdImg.loadPixels();
		if ((scFlags[usePhotonMapIDX]) && (!scFlags[isPhtnMapRndrdIDX])){	if (scFlags[isCausticPhtnIDX]) {sendCausticPhotons(); } else {sendDiffusePhotons();scFlags[isPhtnMapRndrdIDX] = true;}}		
	}//initRender	
	
	/////////////
	////skydome stuff - move to sphere code, set flag for internal normals TODO
	// find corresponding u and v values for background texture	
	public double findSkyDomeT(myRay transRay){
		//similar to intersection of known direction vectors to lights
		//this code finds t of where passed ray hits mySkyDome edge
		double t = -Double.MAX_VALUE;  //this t is the value of the ray equation where it hits the dome - init to bogus value	  
		//find t for intersection : 
		double a = mySkyDome.getAVal(transRay), b = mySkyDome.getBVal(transRay), c = mySkyDome.getCVal(transRay);	    
		double discr = ((b * b) - (4 * a * c));
		//quadratic - check first if imaginary - if so then no intersection
		if (discr > 0){     
			double discr1 = Math.pow(discr,.5), t1 = (-1*b + discr1)/ (2*a), t2 = (-1*b - discr1)/ (2*a), tVal = Math.min(t1,t2);
			if (tVal < p.epsVal){tVal = Math.max(t1,t2);}//if the min t val is less than 0
			t = tVal;
		}//if positive t
		else {System.out.println ("error - non-colliding ray doesn't hit sky dome.  b^2 - 4ac , eval : "+ b + "^2 - 4 " + a + "*" + c + " : " + discr);}//should never get here - any ray that is sent to this function hasn't intersected with any other
		return t;
	}//findSkjDomeT func
	
	// find corresponding u and v values for background texture
	public double findBkgTextureV(myVector isctPt, double t){
		double v = 0.0;
		double a0 = isctPt.y - mySkyDome.origin.y;
		double a1 = a0 /(mySkyDome.radY);
		a1 = (a1 > 1)? 1 : (a1 < -1) ? -1 : a1;   
		v = (((myImageTexture)mySkyDome.shdr.txtr).myTextureBottom.height-1) * Math.acos(a1)/ Math.PI;
	  return v;
	}

	public double findBkgTextureU(myVector isctPt, double v, double t){
		double u = 0.0, q,a0, a1, a2, shWm1 = ((myImageTexture)mySkyDome.shdr.txtr).myTextureBottom.width-1, z1 = (isctPt.z - mySkyDome.origin.z);	  
		q = v/(((myImageTexture)mySkyDome.shdr.txtr).myTextureBottom.height-1);//normalize v to be 0-1
		a0 = (isctPt.x - mySkyDome.origin.x)/ (mySkyDome.radX);
		a0 = (a0 > 1) ? 1 : (a0 < -1) ? -1 : a0;
		a1 = ( Math.sin(q* Math.PI));
		a2 = ( Math.abs(a1) < p.epsVal) ? 1 : a0/a1;
		u = (z1 <= p.epsVal) ? ((shWm1 * ( Math.acos(a2))/ (PConstants.TWO_PI)) + shWm1/2.0f) : 
					shWm1 - ((shWm1 * ( Math.acos(a2))/ (PConstants.TWO_PI)) + shWm1/2.0f);
		u = (u < 0) ? 0 : (u > shWm1) ? shWm1 : u;
		return u;
	} 

	public myColor getBackgroundTextureColor(myRay ray){
		double t = findSkyDomeT(ray);
		myVector isctPt = ray.pointOnRay(t);
		double v = findBkgTextureV(isctPt, t), u = findBkgTextureU(isctPt, v, t);
		return new myColor(((myImageTexture)mySkyDome.shdr.txtr).myTextureBottom.pixels[(int)v * ((myImageTexture)mySkyDome.shdr.txtr).myTextureBottom.width + (int)u]);
	}//getBackgroundTexturecolor
	
	//////////////////
	//flip the normal directions for this scene
	public void flipNormal(){
		scFlags[flipNormsIDX] = !scFlags[flipNormsIDX];
		scFlags[renderedIDX] = false;				//so will be re-rendered
		scFlags[saveImageIDX] =  true;				//save image with flipped norm
		curRefineStep = 0;
		reflRays = 0;
		refrRays = 0;
		globRayCount = 0;
		for (myGeomBase obj : objList){//set for all scene objects or instances of sceneobjects
			if(obj instanceof mySceneObject){((mySceneObject)obj).setFlags(mySceneObject.invertedIDX, scFlags[flipNormsIDX]);}//either a scene object or an instance of a scene object
			else {if(obj instanceof myInstance && ((myInstance)obj).obj instanceof mySceneObject){((mySceneObject)((myInstance)obj).obj).setFlags(mySceneObject.invertedIDX, scFlags[flipNormsIDX]);}}
		}
	}//flipNormal	
	
	//////////////////////
	//draw utility functions
	//////////////////////
	//write span of pixels with same value, for iterative refinement
	public int writePxlSpan(int clrInt, int row, int col, int _numPxls, int[] pxls){
		int pxlIDX = (row*sceneCols) + col;								//idx of pxl in pxls ara
		int rowPxlCnt, rowStart = row, rowEnd = Math.min(rowStart + _numPxls, sceneRows ),//dont try to go beyond pxl array dims
			colStart = col, colEnd = Math.min(colStart + _numPxls, sceneCols);
		for(int pxlR = rowStart; pxlR < rowEnd; ++pxlR){rowPxlCnt = (pxlR*sceneCols);	for(int pxlC = colStart; pxlC < colEnd; ++pxlC){pxls[(rowPxlCnt + pxlC)] = clrInt;}}		
		return pxlIDX;
	}//writePxlSpan
	
	//instance scene-specific 
	//public abstract myColor calcAAColor(double pRayX, double pRayY, double xIncr, double yIncr);
	public abstract myColor shootMultiRays(double pRayX, double pRayY);
	public abstract void draw(); 

	//file
	private void saveFile(){
		now = Calendar.getInstance();
		String tmpSaveName;
		if(saveName.equals("")){saveName=p.gCurrentFile;}				//use cli file name as default save name
		String[] tmp = saveName.split("\\.(?=[^\\.]+$)");				//remove extension from given savename
		//if (scFlags[saveImgInDirIDX]){	tmpSaveName = folderName.toString() + "\\"  + tmp[0]+(scFlags[myScene.flipNormsIDX] ? "_normFlipped" : "")+"_"+getDateTimeString(false,true,"-") + ".png";} //rebuild name to include directory and image name including render time
		if (scFlags[saveImgInDirIDX]){	tmpSaveName = folderName.toString() + "\\"  + tmp[0]+(scFlags[myScene.flipNormsIDX] ? "_normFlipped" : "")+ ".png";} //rebuild name to include directory and image name including render time
		else {							tmpSaveName = tmp[0]+(scFlags[myScene.flipNormsIDX] ? "_normFlipped" : "")+".png";		}
		System.out.println("File saved as  : "+ saveName);
		rndrdImg.save(tmpSaveName);
		scFlags[saveImageIDX] =  false;//don't keep saving every frame
	}//save image
	  
	//common finalizing for all rendering methods
	protected void finishImage(){
		if (scFlags[saveImageIDX]){		saveFile();	}//if savefile is true, save the file
		if (scFlags[showObjInfoIDX]){
			for (myGeomBase obj : allObjsToFind){
	     		System.out.println(obj.toString());
	     		System.out.println();
	     		if(obj instanceof mySceneObject){
		     		if (((mySceneObject)obj).shdr.txtr.txtFlags[myTextureHandler.txtrdTopIDX]){
		     			System.out.println("" + ((mySceneObject)obj).showUV());
		     		}//if textured
	     		}
		     	else if(obj instanceof myInstance){
		     		myInstance inst = (myInstance)obj;
		     		if ((inst.obj instanceof mySceneObject) && (((mySceneObject)inst.obj).shdr.txtr.txtFlags[myTextureHandler.txtrdTopIDX])){			//TODO need to modify this when using instanced polys having textures - each instance will need a notion of where it should sample from
		     			System.out.println("" + ((mySceneObject)inst.obj).showUV());
		     		}	     	
		     	}
	     	}
		}//for objects and instances, to print out info
		if (scFlags[glblTxtrdBkgIDX]){
			System.out.println("\nBackground : \n");
			System.out.println("" + mySkyDome.showUV());  
		}
		System.out.println("total # of rays : " + globRayCount + " | refl/refr rays " + reflRays +"/" + refrRays);
		p.DispEnd();
	}
	
	/**
	*  build translate, scale and rotation matricies to use for ray tracer
	*  need to implement inversion for each matrix - will apply inverses of these matricies to generated ray so that object is rendered in appropriate manner :
	*
	*  so if object A is subjected to a translate/rotate/scale sequence to render A' then to implement this we need to 
	*  subject the rays attempting to intersect with it by the inverse of these operations to find which rays will actually intersect with it.
	*/
	 public void gtDebugStack(String caller){ System.out.println("Caller : "+caller + "\nCurrent stack status : \n"+matrixStack.toString()); }//gtdebugStack method

	 public void gtInitialize() {
		 currMatrixDepthIDX = 0;
		 matrixStack = new myMatStack(this.matStackMaxHeight);
		 matrixStack.initStackLocation(0);
	 }//gtInitialize method

	public void gtPushMatrix() {
		if (currMatrixDepthIDX < matStackMaxHeight){
	    	matrixStack.push();
	    	currMatrixDepthIDX++;
		} else {	System.out.println("Error, matrix depth maximum " + matStackMaxHeight + " exceeded");	}	  
	}//gtPushMatrix method

	public void gtPopMatrix() { 
		if (matrixStack.top == 0){System.out.println("Error : Cannot pop the last matrix in the matrix stack");} 
		else {		//temp was last matrix at top of stack - referencing only for debugging purposes
			myMatrix temp = matrixStack.pop();
			currMatrixDepthIDX--;
		}
	}//gtPopMatrix method

	public void gtTranslate(double tx, double ty, double tz) { 
		//build and push onto stack the translation matrix
		myMatrix TransMat = new myMatrix();
		//set the 4th column vals to be the translation coordinates
		TransMat.setValByIdx(0,3,tx);
		TransMat.setValByIdx(1,3,ty);
		TransMat.setValByIdx(2,3,tz);
		updateCTM(TransMat);
	}//gtTranslate method

	public void gtScale(double sx, double sy, double sz) {
		//build and push onto stack the scale matrix
		myMatrix ScaleMat = new myMatrix();
		//set the diagonal vals to be the scale coordinates
		ScaleMat.setValByIdx(0,0,sx);
		ScaleMat.setValByIdx(1,1,sy);
		ScaleMat.setValByIdx(2,2,sz);
		updateCTM(ScaleMat);
	}//gtScale method

	/**
	*  sets a rotation matrix to be in "angle" degrees CCW around the axis given by ax,ay,az
	*  and multiples this matrix against the CTM
	*/
	public void gtRotate(double angle, double ax, double ay, double az) { 
		// build and add to top of stack the rotation matrix
		double angleRad = (double)(angle * Math.PI)/180.0;
		myMatrix RotMat = new myMatrix();
		myMatrix RotMatrix1 = new myMatrix();      //translates given axis to x axis
		myMatrix RotMatrix2 = new myMatrix();      //rotation around x axis by given angle
		myMatrix RotMatrix1Trans = new myMatrix();
	  
		myVector axisVect, axisVectNorm, bVect, bVectNorm, cVect, cVectNorm, normVect;
		//first build rotation matrix to rotate ax,ay,az to lie in line with x axis		
		axisVect = new myVector(ax,ay,az);
		axisVectNorm = axisVect._normalized();
	  
		if (ax == 0) { 	normVect = new myVector(1,0,0);} 
		else {			normVect = new myVector(0,1,0);}
		bVect = axisVectNorm._cross(normVect);
		bVectNorm = bVect._normalized();
	  
		cVect = axisVectNorm._cross(bVectNorm);
		cVectNorm = cVect._normalized();
	  
		RotMatrix1.setValByRow(0,axisVectNorm);
		RotMatrix1.setValByRow(1,bVectNorm);
		RotMatrix1.setValByRow(2,cVectNorm);
		
		RotMatrix1Trans = RotMatrix1.transpose();
		//second build rotation matrix to rotate around x axis by angle
		//need to set 1,1 ; 1,2 ; 2,1 ; and 2,2 to cos thet, neg sine thet, sine thet, cos thet, respectively
	 
		RotMatrix2.setValByIdx(1,1,(Math.cos(angleRad)));
		RotMatrix2.setValByIdx(1,2,(-Math.sin(angleRad)));
		RotMatrix2.setValByIdx(2,1,(Math.sin(angleRad)));
		RotMatrix2.setValByIdx(2,2,(Math.cos(angleRad)));
		//lastly, calculate full rotation matrix

		myMatrix tmp = RotMatrix2.multMat(RotMatrix1);
		RotMat = RotMatrix1Trans.multMat(tmp);
		updateCTM(RotMat);
	}//gtrotate
	
	public void updateCTM(myMatrix _mat){		
		myMatrix CTM = matrixStack.peek();
		matrixStack.replaceTop(CTM.multMat(_mat));
	}
	/////
	//end matrix stuff
	/////	

	//build a date with each component separated by token
	public String getDateTimeString(){return getDateTimeString(true, false,".");}
	public String getDateTimeString(boolean useYear, boolean toSecond, String token){
		String result = "";
		int val;
		if(useYear){val = now.get(Calendar.YEAR);		result += ""+val+token;}
		val = now.get(Calendar.MONTH)+1;				result += (val < 10 ? "0"+val : ""+val)+ token;
		val = now.get(Calendar.DAY_OF_MONTH);			result += (val < 10 ? "0"+val : ""+val)+ token;
		val = now.get(Calendar.HOUR);					result += (val < 10 ? "0"+val : ""+val)+ token;
		val = now.get(Calendar.MINUTE);					result += (val < 10 ? "0"+val : ""+val);
		if(toSecond){val = now.get(Calendar.SECOND);	result += token + (val < 10 ? "0"+val : ""+val);}
		return result;
	}
	//describe scene
	public String toString(){
		String res = "";
		return res;
	}
}//myScene

class myFOVScene extends myScene {
	//current field of view
	public double fov, fovRad, viewZ;			//degrees, radians,distance from eye the current view plane exists at, in -z direction
	//public List<Future<Boolean>> callFOVFutures;
	//public List<myFOVCall> callFOVCalcs;

	public myFOVScene(DistRayTracer _p) {
		super(_p);
		viewZ = -1;
		//callFOVCalcs= new ArrayList<myFOVCall>();
		//callFOVFutures = new ArrayList<Future<Boolean>>(); 
	}
	public myFOVScene(myScene _scene) {
		super( _scene);
		//callFOVCalcs= new ArrayList<myFOVCall>();
		//callFOVFutures = new ArrayList<Future<Boolean>>(); 
	}

	@Override
	public void setSceneParams(double[] args){
		fov = args[0];
		fovRad = Math.PI*fov/180.0;
		if (Math.abs(fov - 180) < .001){//if illegal fov value, modify to prevent divide by 0
			fov -= .001;
			fovRad -= .0001;
		}
		//virtual view plane exists at some -z so that fov gives sceneCols x sceneRows x and y;-1 for negative z, makes coord system right handed
		//if depth of field scene, then use lens focal distance as z depth (??)
		viewZ = -1 *  (Math.max(sceneRows,sceneCols)/2.0)/Math.tan(fovRad/2);
		if(scFlags[hasDpthOfFldIDX]){//depth of field variables already set from reader - build focal plane
			focalPlane.setPlaneVals(0, 0, 1,  lens_focal_distance);  		
			System.out.println("View z : " + viewZ  + "\nfocal plane : "+ focalPlane);
		}
	}//setSceneParams

	//getDpthOfFldEyeLoc()
	//no anti aliasing for depth of field, instead, first find intersection of ray with focal plane, then 
	//then find start location of ray via getDpthOfFldEyeLoc(), and build multiple rays 
	protected myColor shootMultiDpthOfFldRays(double pRayX, double pRayY) {
		myColor result,aaResultColor;
		double redVal = 0, greenVal = 0, blueVal = 0;//, rayYOffset = sceneRows/2.0, rayXOffset = sceneCols/2.0;
		myVector lensCtrPoint = new myVector(pRayX,pRayY,viewZ);
		lensCtrPoint._normalize();
		myRay ray = new myRay(this, eyeOrigin, lensCtrPoint, 0);					//initial ray - find intersection with focal plane
		//find intersection point with focal plane, use this point to build lens rays
		rayHit hit = focalPlane.intersectCheck( ray, ray.getTransformedRay(ray, focalPlane.CTMara[focalPlane.invIDX]),focalPlane.CTMara);						//should always hit
		myVector rayOrigin,														//
			focalPt = hit.hitLoc;
		for(int rayNum = 0; rayNum < numRaysPerPixel; ++rayNum){
			rayOrigin = this.getDpthOfFldEyeLoc(lensCtrPoint);										//get some random pt within the lens to use as the ray's origin
			ray = new myRay(this, rayOrigin, new myVector(rayOrigin, focalPt),0);
			aaResultColor = reflectRay(ray);
			redVal += aaResultColor.RGB.x; //(aaResultColor >> 16 & 0xFF)/256.0;//gets red value
			greenVal += aaResultColor.RGB.y; // (aaResultColor >> 8 & 0xFF)/256.0;//gets green value
			blueVal += aaResultColor.RGB.z;//(aaResultColor & 0xFF)/256.0;//gets blue value
		}//rayNum
		result = new myColor ( redVal/numRaysPerPixel, greenVal/numRaysPerPixel, blueVal/numRaysPerPixel); 
		return result;	  
	}//shootMultiRays	
	
	protected void drawDpthOfFld(){
		if (!scFlags[renderedIDX]){	
			initRender();
			//index of currently written pixel
			int pixIDX = 0;
			int progressCount = 0;
			double rayY, rayX;
			myColor showColor;
			//myRay ray;
			boolean skipPxl = false;
			int stepIter = 1;
			if(scFlags[glblRefineIDX]){
				stepIter = RefineIDX[curRefineStep++];
				skipPxl = curRefineStep != 1;			//skip 0,0 pxl on all sub-images except the first pass
			} 
			if(stepIter == 1){scFlags[renderedIDX] = true;			}
			for (int row = 0; row < sceneRows; row+=stepIter){
				rayY = (-1 * (row - rayYOffset));        
				for (int col = 0; col < sceneCols; col+=stepIter){
					if(skipPxl){skipPxl = false;continue;}			//skip only 0,0 pxl		
					rayX = col - rayXOffset;      
					showColor = shootMultiDpthOfFldRays(rayX,rayY);
					pixIDX = writePxlSpan(showColor.getInt(),row,col,stepIter,rndrdImg.pixels);
					if ((1.0 * pixIDX)/(numPxls) > (progressCount * .02)){System.out.print("-|");progressCount++;}//progressbar  
				}//for col
			}//for row  
			System.out.println("-");
			//update the display based on the pixels array
			rndrdImg.updatePixels();
			if(scFlags[renderedIDX]){//only do this stuff when finished				
				finishImage();
			}	
		}
		p.imageMode(PConstants.CORNER);
		p.image(rndrdImg,0,0);			
	}//drawDpthOfFld
	
	
	@Override//calculates color based on multiple rays shot into scene
	public myColor shootMultiRays(double xBseVal, double yBseVal) {
		myColor result,aaResultColor;
		double redVal = 0, greenVal = 0, blueVal = 0, rayY, rayX;//, rayYOffset = sceneRows/2.0, rayXOffset = sceneCols/2.0;
		myRay ray;		
		for(int rayNum = 0; rayNum < numRaysPerPixel; ++rayNum){//vary by +/- .5
			rayY = yBseVal + ThreadLocalRandom.current().nextDouble(-.5,.5);
			rayX = xBseVal + ThreadLocalRandom.current().nextDouble(-.5,.5);
			ray = new myRay(this, this.eyeOrigin, new myVector(rayX,rayY,viewZ),0);
			aaResultColor = reflectRay(ray);
			redVal += aaResultColor.RGB.x; //(aaResultColor >> 16 & 0xFF)/256.0;//gets red value
			greenVal += aaResultColor.RGB.y; // (aaResultColor >> 8 & 0xFF)/256.0;//gets green value
			blueVal += aaResultColor.RGB.z;//(aaResultColor & 0xFF)/256.0;//gets blue value
		}//rayNum
		result = new myColor ( redVal/numRaysPerPixel, greenVal/numRaysPerPixel, blueVal/numRaysPerPixel); 
		return result;	  
	}//shootMultiRays	
//		
//	public void drawTmp (){
//		callInitBoidCalcs.clear();
//		myBoid[] tmpList;
//		for(int c = 0; c < boidFlock.length; c+=mtFrameSize){
//			int finalLen = (c+mtFrameSize < boidFlock.length ? mtFrameSize : boidFlock.length - c);
//			tmpList = new myBoid[finalLen];
//			System.arraycopy(boidFlock, c, tmpList, 0, finalLen);			
//			callInitBoidCalcs.add(new myInitPredPreyMaps(p, this, preyFlock, predFlock, fv, tmpList));
//		}							//find next turn's motion for every creature by finding total force to act on creature
//		try {callInitFutures = th_exec.invokeAll(callInitBoidCalcs);for(Future<Boolean> f: callInitFutures) { f.get(); }} catch (Exception e) { e.printStackTrace(); }			
//
//		
//	}
	
	
	@Override
	//distribution draw
	public void draw(){
		if(scFlags[hasDpthOfFldIDX]){drawDpthOfFld(); return;}
		if (!scFlags[renderedIDX]){	
			initRender();
			//index of currently written pixel
			int pixIDX = 0;
			int progressCount = 0;
			double rayY, rayX;
			myColor showColor;
			//myRay ray;
			boolean skipPxl = false;
			int stepIter = 1;
			if(scFlags[glblRefineIDX]){
				stepIter = RefineIDX[curRefineStep++];
				skipPxl = curRefineStep != 1;			//skip 0,0 pxl on all sub-images except the first pass
			} 
			if(stepIter == 1){scFlags[renderedIDX] = true;			}
			if (numRaysPerPixel == 1){//only single ray shot into scene for each pixel
				for (int row = 0; row < sceneRows; row+=stepIter){
					rayY = (-1 * (row - rayYOffset));         
					for (int col = 0; col < sceneCols; col+=stepIter){
						if(skipPxl){skipPxl = false;continue;}			//skip only 0,0 pxl					
						rayX = col - rayXOffset;
						showColor = reflectRay(new myRay(this,this.eyeOrigin, new myVector(rayX,rayY,viewZ),0)); 
						pixIDX = writePxlSpan(showColor.getInt(),row,col,stepIter,rndrdImg.pixels);
						if ((1.0 * pixIDX)/(numPxls) > (progressCount * .02)){System.out.print("-|");progressCount++;}//progressbar         
					}//for col
				}//for row	     
			} else{    //multiple rays shot into scene per pxl
				for (int row = 0; row < sceneRows; row+=stepIter){
					rayY = (-1 * (row - rayYOffset));        
					for (int col = 0; col < sceneCols; col+=stepIter){
						if(skipPxl){skipPxl = false;continue;}			//skip only 0,0 pxl		
						rayX = col - rayXOffset;      
						showColor = shootMultiRays(rayX,rayY);
						pixIDX = writePxlSpan(showColor.getInt(),row,col,stepIter,rndrdImg.pixels);
						if ((1.0 * pixIDX)/(numPxls) > (progressCount * .02)){System.out.print("-|");progressCount++;}//progressbar  
					}//for col
				}//for row  
			}//if antialiasing
			
			System.out.println("-");
			//update the display based on the pixels array
			rndrdImg.updatePixels();
			if(scFlags[renderedIDX]){//only do this stuff when finished				
				finishImage();
			}	
		}
		p.imageMode(PConstants.CORNER);
		p.image(rndrdImg,0,0);		
	}
	
}//myFOVScene

class myFishEyeScene extends myScene{
	//current field of view
	public double fishEye, fishEyeRad;			//fisheye in degrees, radians	
	private double aperatureHlf;
	//public List<Future<Boolean>> callFishFutures;
	//public List<myFishCall> callFishCalcs;
	
	public myFishEyeScene(DistRayTracer _p) {
		super(_p);		
		//callFishCalcs= new ArrayList<myFishCall>();
		//callFishFutures = new ArrayList<Future<Boolean>>(); 

	}
	
	public myFishEyeScene(myScene _s){
		super(_s);
		//callFishCalcs= new ArrayList<myFishCall>();
		//callFishFutures = new ArrayList<Future<Boolean>>(); 
	}
	
	@Override
	public void setSceneParams(double[] args) {
		fishEye = args[0];
		fishEyeRad = Math.PI*fishEye/180.0;
		aperatureHlf = fishEyeRad/2.0;
	}//setSceneParams

	@Override
	public myColor shootMultiRays(double xBseVal, double yBseVal) {
		myColor result,aaResultColor;
		double xVal, yVal, r, rSqTmp, theta, phi, redVal = 0, greenVal = 0, blueVal = 0;		
		myRay ray;
		for(int rayNum = 0; rayNum < numRaysPerPixel; ++rayNum){//vary by +/- .5
			yVal = (yBseVal + ThreadLocalRandom.current().nextDouble(-.5,.5)) *fishMult;
			xVal = (xBseVal + ThreadLocalRandom.current().nextDouble(-.5,.5)) *fishMult; 
			rSqTmp = yVal * yVal + xVal*xVal;
			if(rSqTmp <= 1){
				r = Math.sqrt(rSqTmp);
				theta = r * aperatureHlf;
				phi = Math.atan2(-yVal,xVal); 					
				double sTh = Math.sin(theta);
				ray = new myRay(this,this.eyeOrigin, new myVector(sTh * Math.cos(phi),sTh * Math.sin(phi),-Math.cos(theta)),0);
				aaResultColor = reflectRay(ray);
				redVal += aaResultColor.RGB.x; //(aaResultColor >> 16 & 0xFF)/256.0;//gets red value
				greenVal += aaResultColor.RGB.y; // (aaResultColor >> 8 & 0xFF)/256.0;//gets green value
				blueVal += aaResultColor.RGB.z;//(aaResultColor & 0xFF)/256.0;//gets blue value
			}			
		}//rayNum
		result = new myColor ( redVal/numRaysPerPixel, greenVal/numRaysPerPixel, blueVal/numRaysPerPixel); 
		return result;	  
	}//shootMultiRays
	
	@Override
	//distribution draw
	public void draw(){
		if (!scFlags[renderedIDX]){
			initRender();
			//index of currently written pixel
			int pixIDX = 0;
			int progressCount = 0;
			double r, phi, theta,yVal, xVal, ySq, sTh, rTmp;
			//fishRad2 = .5*fishEyeRad;
			myColor showColor;
			boolean skipPxl = false;
			int stepIter = 1;
			if(scFlags[glblRefineIDX]){
				stepIter = RefineIDX[curRefineStep++];
				skipPxl = curRefineStep != 1;			//skip 0,0 pxl on all sub-images except the first pass
			} 
			if(stepIter == 1){scFlags[renderedIDX] = true;			}
			//fisheye assumes plane is 1 away from eye
			if (numRaysPerPixel == 1){											//single ray into scene per pixel
				for (int row = 0; row < sceneRows; row+=stepIter){
					yVal = (row + yStart) * fishMult;
					ySq = yVal * yVal;
					for (int col = 0; col < sceneCols; col+=stepIter){
						if(skipPxl){skipPxl = false;continue;}			//skip only 0,0 pxl	
						xVal = (col + xStart)* fishMult;
						rTmp = xVal*xVal+ySq;
						if(rTmp > 1){	pixIDX = writePxlSpan(blkColor.getInt(),row,col,stepIter,rndrdImg.pixels);	} 
						else {
							r = Math.sqrt(rTmp);
							theta = r * aperatureHlf;
							phi = Math.atan2(-yVal,xVal); 					
							sTh = Math.sin(theta);
							showColor = reflectRay(new myRay(this,this.eyeOrigin, new myVector(sTh * Math.cos(phi),sTh * Math.sin(phi),-Math.cos(theta)),0)); 
							pixIDX = writePxlSpan(showColor.getInt(),row,col,stepIter,rndrdImg.pixels);
						}
						if ((1.0 * pixIDX)/(numPxls) > (progressCount * .02)){System.out.print("-|");progressCount++;}//progressbar  
					}//for col
				}//for row	     
			} else{    
				for (int row = 0; row < sceneRows; row+=stepIter){
					yVal = (row + yStart);
					for (int col = 0; col < sceneCols; col+=stepIter){
						if(skipPxl){skipPxl = false;continue;}			//skip only 0,0 pxl		
						xVal = (col + xStart);
						showColor = shootMultiRays(xVal, yVal); 			//replace by base radian amt of max(x,y) 
						pixIDX = writePxlSpan(showColor.getInt(),row,col,stepIter,rndrdImg.pixels);			
						if ((1.0 * pixIDX)/(numPxls) > (progressCount * .02)){System.out.print("-|");progressCount++;}//progressbar  
					}//for col
				}//for row  
			}//if antialiasing
			System.out.println("-");
			if(scFlags[renderedIDX]){			finishImage();	}	
		}		
		p.imageMode(PConstants.CORNER);
		p.image(rndrdImg,0,0);			
	}//draw

}//myFishEyeScene

class myOrthoScene extends myScene{
	//width and height of view - for ortho projection. for perspective will be screen width and height
	public double orthoWidth, orthoHeight;	
	private double orthPerRow, orthPerCol;			//normalizers for ortho projection
	//public List<Future<Boolean>> callOrthoFutures;
	//public List<myOrthoCall> callOrthoCalcs;

	public myOrthoScene(DistRayTracer _p) {
		super(_p);
		orthoWidth = p.sceneCols;
		orthoHeight = p.sceneRows;
		double div = Math.min(p.sceneCols, p.sceneRows);
		orthPerRow = orthoHeight/div;
		orthPerCol = orthoWidth/div;
		//callOrthoCalcs= new ArrayList<myOrthoCall>();
		//callOrthoFutures = new ArrayList<Future<Boolean>>(); 
		
		
	}
	public myOrthoScene(myScene _s){
		super(_s);
		orthoWidth = p.sceneCols;
		orthoHeight = p.sceneRows;
		double div = Math.min(p.sceneCols, p.sceneRows);
		orthPerRow = orthoHeight/div;
		orthPerCol = orthoWidth/div;		
		//callOrthoCalcs= new ArrayList<myOrthoCall>();
		//callOrthoFutures = new ArrayList<Future<Boolean>>(); 
	}

	@Override
	public void setSceneParams(double[] args) {
		orthoWidth = args[0];
		orthoHeight = args[1];	
		double div = Math.min(p.sceneCols, p.sceneRows);
		orthPerRow = orthoHeight/div;
		orthPerCol = orthoWidth/div;
	}//setSceneParams
	
	@Override
	public myColor shootMultiRays(double xBseVal, double yBseVal) {
		myColor result,aaResultColor;
		double redVal = 0, greenVal = 0, blueVal = 0, rayY, rayX;//,rayYOffset = 1.0/sceneRows, rayXOffset = 1.0/sceneCols;
		myRay ray;
		for(int rayNum = 0; rayNum < numRaysPerPixel; ++rayNum){//vary by +/- .5
			rayY = yBseVal + (orthPerRow*ThreadLocalRandom.current().nextDouble(-.5,.5));
			rayX = xBseVal + (orthPerCol*ThreadLocalRandom.current().nextDouble(-.5,.5));				
			aaResultColor = reflectRay(new myRay(this, new myVector(rayX,rayY,0), new myVector(0,0,-1),0));
			redVal += aaResultColor.RGB.x; //(aaResultColor >> 16 & 0xFF)/256.0;//gets red value
			greenVal += aaResultColor.RGB.y; // (aaResultColor >> 8 & 0xFF)/256.0;//gets green value
			blueVal += aaResultColor.RGB.z;//(aaResultColor & 0xFF)/256.0;//gets blue value	      
		}//aaliasR
		result = new myColor ( redVal/numRaysPerPixel, greenVal/numRaysPerPixel, blueVal/numRaysPerPixel); 
		return result;
	}//shootMultiRays	
	
	@Override
	public void draw(){
		if (!scFlags[renderedIDX]){	
			initRender();
			//we must shoot out rays and determine what is being hit by them
			//get pixels array, to be modified by ray-shooting
			//index of currently written pixel
			int pixIDX = 0;
			int progressCount = 0;
			//  double redVal, greenVal, blueVal, divVal;
			double rayY, rayX;
			double rayYOffset = sceneRows/2.0, rayXOffset = sceneCols/2.0;
			myColor showColor;
			boolean skipPxl = false;
			int stepIter = 1;
			if(scFlags[glblRefineIDX]){
				stepIter = RefineIDX[curRefineStep++];
				skipPxl = curRefineStep != 1;			//skip 0,0 pxl on all sub-images except the first pass
			} 
			if(stepIter == 1){scFlags[renderedIDX] = true;			}
			if (numRaysPerPixel == 1){											//single ray into scene per pixel
				for (int row = 0; row < sceneRows; row+=stepIter){
					rayY = orthPerRow * (-1 * (row - rayYOffset));         
					for (int col = 0; col < sceneCols; col+=stepIter){
						if(skipPxl){skipPxl = false;continue;}			//skip only 0,0 pxl					
						rayX = orthPerCol * (col - rayXOffset);
						showColor = reflectRay(new myRay(this,new myVector(rayX,rayY,0), new myVector(0,0,-1),0)); 
						pixIDX = writePxlSpan(showColor.getInt(),row,col,stepIter,rndrdImg.pixels);
						if ((1.0 * pixIDX)/(numPxls) > (progressCount * .02)){System.out.print("-|");progressCount++;}//progressbar         
					}//for col
				}//for row	     
			} else{    //anti aliasing
				for (int row = 0; row < sceneRows; row+=stepIter){
					rayY = orthPerRow * ((-1 * (row - rayYOffset)) - .5);         
					for (int col = 0; col < sceneCols; col+=stepIter){
						if(skipPxl){skipPxl = false;continue;}			//skip only 0,0 pxl		
						rayX = orthPerCol * (col - rayXOffset - .5);      
						showColor = shootMultiRays(rayX,rayY); 
						pixIDX = writePxlSpan(showColor.getInt(),row,col,stepIter,rndrdImg.pixels);
						if ((1.0 * pixIDX)/(numPxls) > (progressCount * .02)){System.out.print("-|");progressCount++;}//progressbar  
					}//for col
				}//for row  
			}//if antialiasing			
			System.out.println("-");
			//update the display based on the pixels array
			rndrdImg.updatePixels();
			if(scFlags[renderedIDX]){	finishImage();	}	
		}
		p.imageMode(PConstants.CORNER);
		p.image(rndrdImg,0,0);			
	}//draw	

}//myOrthoScene
