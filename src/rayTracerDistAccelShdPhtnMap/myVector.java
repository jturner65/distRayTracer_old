package rayTracerDistAccelShdPhtnMap;

/**
*  replacement class for pvector 
*  keeps sqmagnitude, not magnitude - need to calculate magnitude
*/
public class myVector{
	public double x,y,z;
	public double sqMagn;
		//vector constants available to all consumers of myVector
	myVector(double x, double y, double z){this.x = x; this.y = y; this.z = z;  this._mag();}         //constructor 3 args  
	myVector(double [] v){this.x = v[0]; this.y = v[1]; this.z = v[2];  this._mag();}         //constructor ara arg
	myVector(myVector p){ this.x = p.x; this.y = p.y; this.z = p.z;  this._mag(); }                                                                                                           //constructor 1 arg  
	myVector(myVector p, myVector q){ this(q.x-p.x, q.y-p.y, q.z-p.z); }    	//constructor 2 arg  vector from p to q
	myVector(){ this(0,0,0);}                                                                                                                               //constructor 0 args
	
	public void set(double x, double y, double z){ this.x = x;  this.y = y;  this.z = z; this._mag(); }                                                                     //set 3 args 
	public void set(myVector p){ this.x = p.x; this.y = p.y; this.z = p.z;  this._mag();}                                                                                //set 1 args
	
	public myVector _mult(double n){ this.x *= n; this.y *= n; this.z *= n; this._mag(); return this; }                                                                                    //_mult 3 args  
	public void _add(double x, double y, double z){ this.x += x; this.y += y; this.z += z;  this._mag(); }                                                                   //_add 3 args
	public void _add(myVector v){ this.x += v.x; this.y += v.y; this.z += v.z;  this._mag(); this._mag(); }                                                                           //_add 1 arg  

	public void _sub(double x, double y, double z){ this.x -= x; this.y -= y; this.z -= z;  this._mag(); }                                                                   //_sub 3 args
	public void _sub(myVector v){ this.x -= v.x; this.y -= v.y; this.z -= v.z;  this._mag(); }                                                                           //_sub 1 arg 
	
	public double _mag(){ return Math.sqrt(this._SqMag());}  
	public double _SqMag(){ this.sqMagn = ((this.x*this.x) + (this.y*this.y) + (this.z*this.z)); return this.sqMagn; }  							//squared magnitude
	
	public void _normalize(){double magn = _mag();if(magn == 0){return;} _div(magn);}
	
	public myVector _normalized(){double magn = this._mag(); myVector newVec = (magn == 0) ? (new myVector(0,0,0)) : (new myVector( this.x /= magn, this.y /= magn, this.z /= magn)); newVec._mag(); return newVec;}

	public myVector cloneMe(){myVector retVal = new myVector(this.x, this.y, this.z); retVal._mag(); return retVal;}  
	
	public double _L1Dist(myVector q){return Math.abs(this.x - q.x) + Math.abs(this.y - q.y) + Math.abs(this.z - q.z); }
	public double _dist(myVector q){ return (double)Math.sqrt( ((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z)) ); }
	
	public void _div(double q){this.x /= q; this.y /= q; this.z /= q; this._mag();}  
	
	public myVector _cross(myVector b){ return new myVector((this.y * b.z) - (this.z*b.y), (this.z * b.x) - (this.x*b.z), (this.x * b.y) - (this.y*b.x));}		//cross product 
	
	public double _dot(myVector b){return ((this.x * b.x) + (this.y * b.y) + (this.z * b.z));}																	//dot product
	
	public double[] getAsAra(){return new double[]{this.x, this.y, this.z};}
	public double[] getAsHAraPt(){return new double[]{this.x, this.y, this.z,1};}
	public double[] getAsHAraVec(){return new double[]{this.x, this.y, this.z,0};}
	/**
	 * returns if this vector is equal to passed vector
	 * @param b vector to check
	 * @return whether they are equal
	 */
	public boolean equals(Object b){
		if (this == b) return true;
		if (!(b instanceof myVector)) return false;
		myVector v = (myVector)b;
		return ((this.x == v.x) && (this.y == v.y) && (this.z == v.z));		
	}	
	
	public String toStrBrf(){return "|(" + String.format("%.4f", x) + "," + String.format("%.4f", y) + "," + String.format("%.4f", z) + ")";}	
	public String toString(){return "|(" + this.x + ", " + this.y + ", " + this.z+")| sqMag:" + this.sqMagn;}

}//myVector

class myMatrix {	  
	public double[][] m;
	
	public myMatrix(){  m = new double[4][4]; initMat();}	
	public void initMat(){  this.initMat(false);}	
	//initialize this matrix to be identity matrix
	//clear= true gives all 0's | clear= false gives identity
	void initMat(boolean clear){for (int row = 0; row < 4; ++row){  for (int col = 0; col < 4; ++col){m[row][col] =  ((row == col) && !clear)  ? 1 : 0 ;}}} 
	
	//multiplies this matrix by b, in order: [this] x [b]
	//returns result in result
	public myMatrix multMat(myMatrix b){
		double resultVal = 0;
		myMatrix result = new myMatrix();
		for (int row = 0; row < 4; ++row){for (int col = 0; col < 4; ++col){for (int k = 0; k < 4; k++){resultVal += this.m[row][k] * b.getValByIdx(k,col);}result.setValByIdx(row,col,resultVal); resultVal = 0;}}
		return result;  
	}//mult method
		
	//multiplies this matrix by vertex b, in order: [this] x [b]
	//returns result vertex in result
	public double[] multVert(double[] b){
		double resultVal;
		double[] result = new double[]{0,0,0,0};
		for (int row = 0; row < 4; ++row){resultVal = 0;for (int col = 0; col < 4; ++col){resultVal += this.m[row][col] * b[col];}	result[row] = resultVal;}//for row
		return result;  
	}//mult method
	
	//returns the transpose of this matrix - also inverse if rotation matrix
	public myMatrix transpose(){ 
		myMatrix result = new myMatrix();
		for (int row = 0; row < 4; ++row){for (int col = 0; col < 4; ++col){result.m[col][row] = this.m[row][col];}}
		return result;   
	}//transpose method

	public myMatrix inverse(){
		myMatrix result = this.InvertMe();
	  return result;    
	}//method invert
	public myMatrix adjoint(){
		myMatrix result = this.InvertMe();
		
		return result.transpose();
	}
	
	//------------- inversion code
	
	private myMatrix InvertMe(){
		double[] tmp = new double[12];   // temp array for pairs
		double[] src = new double[16];   // array of transpose source matrix
		double[] dst = new double[16];   //destination matrix, in array form
		double det = 0;       // determinant 
		myMatrix dstMat = new myMatrix();
		//convert this matrix to array form and set up source vector
		for(int row = 0; row < 4; ++row){ for(int col = 0; col < 4; ++col){ src[(4*col) + row] = this.m[row][col];	}}
		
		// calculate pairs for first 8 elements (cofactors)
		tmp[0] = src[10] * src[15];
		tmp[1] = src[11] * src[14];
		tmp[2] = src[9] * src[15];
		tmp[3] = src[11] * src[13];
		tmp[4] = src[9] * src[14];
		tmp[5] = src[10] * src[13];
		tmp[6] = src[8] * src[15];
		tmp[7] = src[11] * src[12];
		tmp[8] = src[8] * src[14];
		tmp[9] = src[10] * src[12];
		tmp[10] = src[8] * src[13];
		tmp[11] = src[9] * src[12];
		
		// calculate first 8 elements (cofactors) 	  
		dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
		dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
		dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
		dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
		dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
		dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
		dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
		dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
		dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
		dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
		dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
		dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
		dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
		dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
		dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
		dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
		
		// calculate pairs for second 8 elements (cofactors) 
		
		tmp[0] = src[2]*src[7];
		tmp[1] = src[3]*src[6];
		tmp[2] = src[1]*src[7];
		tmp[3] = src[3]*src[5];
		tmp[4] = src[1]*src[6];
		tmp[5] = src[2]*src[5];
		tmp[6] = src[0]*src[7];
		tmp[7] = src[3]*src[4];
		tmp[8] = src[0]*src[6];
		tmp[9] = src[2]*src[4];
		tmp[10] = src[0]*src[5];
		tmp[11] = src[1]*src[4];
		
		// calculate second 8 elements (cofactors)
		
		dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
		dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
		dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
		dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
		dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
		dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
		dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
		dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
		dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
		dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];
		dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];
		dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
		dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];
		dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];
		dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];
		dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
		
		// calculate determinant 		
		det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];
		if(Math.abs(det) > .0000001){		
			for (int j = 0; j < 16; j++){ dst[j] /= det; }	    
			 //convert dst array to matrix
			for(int row = 0; row < 4; ++row){ for(int col = 0; col < 4; ++col){  dstMat.m[row][col] = dst[(4*row) + col];}}
		} else {
			System.out.println("uninvertible matrix -> det == 0");
		}
		return dstMat;
	}//invertme code	
	//end------------inversion code
	
	public double getValByIdx(int row, int col){   return m[row][col]; }  
	public void setValByIdx(int row, int col, double val){	   m[row][col] = val; }
	
		//writes over the first 3 cols of row given by row var with vector vals  
	public void setValByRow(int row, myVector vect){ this.m[row][0] = vect.x; this.m[row][1] = vect.y;      this.m[row][2] = vect.z; }//setValByRow
	 
	//makes a deep copy of this matrix, which it returns
	public myMatrix clone(){
		myMatrix newMat = new myMatrix();
		for (int row = 0; row < 4; ++row){ for (int col = 0; col < 4; ++col){newMat.m[row][col] = this.m[row][col];}}
		return newMat;
	}
	
	public String toString(){
	    String result = "", tmp2str = "",tmpString;
	    for (int row = 0; row < 4; ++row){
	    	result += "[";
	    	for (int col = 0; col < 4; ++col){   tmp2str = "" + m[row][col]; if (col != 3) {tmp2str += ", ";} result += tmp2str;}
	    	tmpString = "]";  if (row != 3) { tmpString += "\n"; }
	    	result += tmpString;
    	}//for row	   
	    return result;
	}//toString method
  
}//class matrix

class myMatStack {
	myMatrix[] s;
	int top;
	 
	public myMatStack(int matStackMaxHeight){
		s = new myMatrix[matStackMaxHeight];
		for (int row = 0; row < 10; ++row){s[row] = new myMatrix();}  	
		top = 0;        //point top of stack at index of base matrix
	}//stack constructor	
	public void initStackLocation(int idx){  this.initStackLocation(idx, false);}	
	public void initStackLocation(int idx, boolean clear){  s[idx].initMat(clear);}	
	//add the current top of the matrix stack to the matrix stack in a higher position
	public void push(){ ++top; initStackLocation(top);  for (int row = 0; row < 4; ++row){ for (int col = 0; col < 4; ++col){  s[top].m[row][col] = s[top - 1].m[row][col]; }}}//push     	
	//replace the current top of the matrix stack with a new matrix
	public void replaceTop(myMatrix newTopMatrix){for (int row = 0; row < 4; ++row){ for (int col = 0; col < 4; ++col){ s[top].m[row][col] = newTopMatrix.m[row][col]; }}}//replaceTop	
	//return the top of the matrix stack without popping
	public myMatrix peek(){ return s[top].clone();	}//peek	
	//remove and return top matrix on stack
	public myMatrix pop(){
		myMatrix oldTop = new myMatrix();
		oldTop = s[top].clone();
		if (top > 0) {	top--;		initStackLocation(top+1,true);    }  //reinitialize stack			
		else {		System.out.println("stack pop error");}
		return oldTop;
	}	
	//returns a string representation of this stack
	public String toString(){
		String result = "";
		for (int si = 0; si < top+1; si++){	result += "Stack[" + si + "] =\n[" + s[si].toString() + "]\n";}//for si
		return result;
	}//to String method
}//gtStack