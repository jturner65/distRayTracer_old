package rayTracerDistAccelShdPhtnMap;
import processing.core.PImage;


public abstract class mySceneObject extends myGeomBase{	
	
	public boolean[] rFlags;					//various state-related flags for this object
	public static final int 
			invertedIDX			= 0,				//normals point in or out
			isLightIDX			= 1//,
			//isTemplateObjIDX	= 2
			;				//this object is used to instance multiple objects (doesn't exist in scene itself)
	public static final int numFlags = 2;	

	public mySceneObject(myScene _scn, double _x, double _y, double _z){
		super(_scn, _x, _y, _z);	    
	    initFlags();
	    shdr = new myObjShader(scene);//sets currently defined colors for this object also    
	}//constructor 6 var
		
	public void initFlags(){rFlags = new boolean[numFlags];for(int i=0; i<numFlags;++i){rFlags[i]=false;}}
	
	public void setFlags(int idx, boolean val){
		rFlags[idx]=val;
		switch(idx){
		case invertedIDX		:{break;}		
		case isLightIDX			:{break;}
	//	case isTemplateObjIDX	:{break;}			//this object is used to instance multiple objects (doesn't exist in scene itself)
		}
	}
	//determine if shadow ray is a hit or not
	@Override
	public int calcShadowHit(myRay _ray, myRay transRay, gtMatrix[] _ctAra, double distToLight){		
		rayHit hitChk = intersectCheck(_ray,transRay,_ctAra);			
		if (hitChk.isHit && (distToLight - hitChk.t) > scene.p.epsVal){	return 1;}  
		//if (distToLight - hitChk.t > scene.epsVal){	return 1;}  
		return 0;
	}
	
	public String showUV(){
		String result = "ObjectType : " + type + " ID : " + ID ;//+ "\nMIN u&v :" + minU + " | " + minV + " | MAX u&v " + maxU + " | " + maxV;
		result += "\ncenter : " + origin ;
		result += shdr.showUV();
		return result;
	}//showUV
	@Override
	public myColor getColorAtPos(rayHit hit) {	return shdr.getColorAtPos(hit);}

	public String toString(){
		String result = super.toString() +"\n"+shdr.toString();
//		if(rFlags[isTemplateObjIDX]){
//			result += "Template object not instanced\n";
//		}
		return result;
	}  
	
}//class mySceneObject
//a scene object representing a box - basically a bounding box + a shader
class myRndrdBox extends mySceneObject{

	public myRndrdBox(myScene _scn, double _x, double _y, double _z, myVector _minVals, myVector _maxVals) {
		super(_scn, _x, _y, _z);
	    minVals.set(_minVals);
	    maxVals.set(_maxVals);
	    type = objType.RenderedBBox;
		postProcBBox();				//cnstrct and define bbox
	}
	public double[] findTxtrCoords(myVector isctPt, PImage myTexture, double time){
		double v = findTextureV(isctPt,myTexture,time);	
		return new double[]{findTextureU(isctPt,v,myTexture,time),v};
	}

	protected double findTextureU(myVector isctPt, double v, PImage myTexture, double time) {		//TODO
		// TODO Auto-generated method stub
		return 0;
	}

	protected double findTextureV(myVector isctPt, PImage myTexture, double time) {					//TODO
		// TODO Auto-generated method stub
		return 0;
	}
	@Override
	public myVector getOrigin(double t) {return origin;}
	@Override
	public myVector getMaxVec() {		return _bbox.getMaxVec();}
	@Override
	public myVector getMinVec() {		return _bbox.getMinVec();}
	@Override
	public rayHit intersectCheck(myRay _ray,myRay transRay, gtMatrix[] _ctAra) {return _bbox.intersectCheck(_ray, transRay, _ctAra);	}
	@Override
	public myVector getNormalAtPoint(myVector point, int[] args) {return _bbox.getNormalAtPoint(point, args);}
}

//an instance of an instanced object - used to minimize memory footprint - only have a reference of the object, and then the relevant transformation matrices
class myInstance extends myGeomBase{
	public myGeomBase obj;						//object this is instance of
	public boolean useShader, isAccel;					//whether to use instance shader or base object shader
	public myInstance(myScene scene, myGeomBase _obj){
		super(scene, 0,0,0);
		obj = _obj;										//owning object
		isAccel = (obj instanceof myAccelStruct);
		CTMara = scene.p.buildCTMara(scene, obj.CTMara[glblIDX]);//build this object's transformation matrix - since this is instancing the owning object, pass the owning object's matrix
	    //CTMara = scene.p.buildCTMara(scene);//build this object's transformation matrix		    
	    type = objType.Instance;//"Instance of "+obj.objType;
	    this.minVals = getMinVec();
	    this.maxVals = getMaxVec();
		postProcBBox();				//cnstrct and define bbox
		shdr = null;
		useShader = false;
	}
	
	//return vector with maximum x/y/z coords of this object
	public myVector getMaxVec(){return scene.p.getTransformedPt(obj.getMaxVec(), obj.CTMara[glblIDX]);}
	//return vector with minimum x/y/z coords of this object
	public myVector getMinVec(){return scene.p.getTransformedPt(obj.getMinVec(), obj.CTMara[glblIDX]);}
	@Override
	public int calcShadowHit(myRay _ray,myRay _trans, gtMatrix[] _ctAra, double distToLight) {
		return obj.calcShadowHit(_trans, _trans, _ctAra, distToLight);
	}
	@Override
	public myVector getOrigin(double t) {	return scene.p.getTransformedPt(obj.getOrigin(t), obj.CTMara[glblIDX]);}
	@Override
	public rayHit intersectCheck(myRay _ray,myRay transRay, gtMatrix[] _ctAra) {
		rayHit _hit = obj.intersectCheck(transRay, transRay, _ctAra);		//copy trans ray over so that ctm-transformed accel structs will still register hits appropriately TODO make this better
		if(useShader){_hit.shdr = shdr;}
		return _hit;	
	}
	@Override
	public myVector getNormalAtPoint(myVector point, int[] args) {		return obj.getNormalAtPoint(point, args);		}	
	//set to override base object's shader
	public void useInstShader(){	useShader = true; shdr = scene.getCurShader();}//new myObjShader(scene);	}	
	//this is probably not the best way to do this - each instance needs to have its own UV coords.  TODOscene.p.getTransformedPt(isctPt, CTMara[invIDX]);
	@Override
	//public double[] findTxtrCoords(myVector isctPt, PImage myTexture, double time){ return obj.findTxtrCoords(isctPt, myTexture, time);}
	public double[] findTxtrCoords(myVector isctPt, PImage myTexture, double time){ 
		return obj.findTxtrCoords(scene.p.getTransformedVec(isctPt, CTMara[invIDX]), myTexture, time);}

	@Override
	public myColor getColorAtPos(rayHit hit) {		return (useShader) ? shdr.getColorAtPos(hit) : hit.obj.shdr.getColorAtPos(hit);}//return (hit.obj instanceof mySceneObject) ? ((mySceneObject)hit.obj).shdr.getColorAtPos(hit) : new myColor(-1,-1,-1,1);	}
	public String toString(){
		String result = super.toString() + " which is an instance of : "+obj.type;
	    return result;
	}

}//myInstance

