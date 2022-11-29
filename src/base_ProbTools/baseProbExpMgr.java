package base_ProbTools;

import java.time.Instant;
import java.util.*;

import base_ProbTools.randGenFunc.transform.*;
import base_StatsTools.summary.myProbSummary_Dbls;
import base_ProbTools.quadrature.myGaussLegenQuad;
import base_ProbTools.quadrature.base.baseQuadrature;
import base_ProbTools.randGenFunc.funcs.*;
import base_ProbTools.randGenFunc.funcs.base.baseRandVarFunc;
import base_ProbTools.randGenFunc.gens.myBoundedRandGen;
import base_ProbTools.randGenFunc.gens.myFleishUniVarRandGen;
import base_ProbTools.randGenFunc.gens.myZigRandGen;
import base_ProbTools.randGenFunc.gens.base.myRandGen;
import base_Utils_Objects.dataAdapter.Base_UIDataUpdater;
import base_Utils_Objects.io.messaging.MessageObject;

public abstract class baseProbExpMgr {
	/**
	 * handle communicating via console or logs
	 */
	public MessageObject msgObj;
	
	/**
	 * Data adapter that receives and sends UI-sourced data with the window that owns this experiment
	 */
	protected Base_UIDataUpdater uiDataValues;
	
	////////////////////////////////////////
	// gauss quadrature solver structures	
	/**
	 * integral solvers for gaussian quadrature method used by this experiment
	 */
	protected baseQuadrature[] quadSlvrs;	
	/**
	 * types of solvers for comparison
	 */
	public static final int
		GL_QuadSlvrIDX 		= 0;			//gaussian legendre solver
	public static final int numGaussSlvrs = 1;
	public static final String[] quadSlvrNames = new String[] {"Gauss-Legendre"};
	
	//index of current solver
	protected int curQuadSolverIDX = 0;	
	//# of points/wts used for gaussian quad
	public int numQuadPoints = 50;
	//convergence tolerance to use for iterative derivation methods
	public double quadConvTol = 1.0e-16;	
	//quadrature calculation big decimal scale to use
    public int quadBDScale = 18;
    
	////////////////////////////////////////
    // random var generators - take uniform input and return desired distribution
	
	/**
	 * types of random number generators implemented/supported so far
	 */
	public static final int 
		boundedRandGen			= 0,	//uses bounded pdf function (i.e. cosine)
		ziggRandGen 			= 1,	//ziggurat method to generate distribution
		fleishRandGen_UniVar	= 2,	//uses fleishman algorithm for univariate - needs first 4 moments
		linearTransformMap		= 3,	//just performs linear transformation mapping - does not actually represent a random function
		uniformTransformMap		= 4;	//just performs an order-based transformation from original grade to #/n where # is order in rank (from lowest to highest) and n is total # of grades
    public static final String[] randGenAlgNames = new String[] {"Bounded PDF Algorithm", "Ziggurat Algorithm", "Fleishman Univariate Polynomial Algorithm", "Linear Transformation Mapping", "Uniform Transformation Mapping"};
	//per algorithm type settings
    protected int[][] randAlgOptSettings;
    
    /**
     * type of random variable function to use with zig or bounded rand gen algorithm
     */
    public static final int
    	normRandVarIDX			= 0,
    	gaussRandVarIDX			= 1,
    	raisedCosRandVarIDX		= 2,
		cosCDFRandVarIDX 		= 3,
		fleishRandVarIDX		= 4;
    //public static final String[] zigRandVarFuncNames = new String[] {"Normal Distribution", "Gaussian Distribution", "Cosine-PDF Distribution", "Cosine-PDF via sample-derived CDF"};
    	
	/**
	 * type of experiment to conduct
	 */
	public static final String[] expType = new String[] {"Gaussian","Linear", "Uniform Spaced","Fleishman Poly","Raised Cosine PDF","Cosine CDF derived"};
	/**
	 * type of plots to show
	 */
	public static final String[] plotType = new String[] {"PDF", "Histogram","CDF (integral)","Inverse CDF"};
	
  
	////////////////////////////////////////
	// internal functionality
	
	//time of current process start, from initial construction of mapmgr - TODO use this to monitor specific process time elapsed.  set to 0 at beginning of a particular process, then measure time elapsed in process
	protected long curProcStartTime;
	//time BaseProbExpMgr built, in millis - used as offset for instant to provide smaller values for timestamp
	protected final long expMgrBuiltTime;	
	//state flags used by instancing experiments
	protected int[] stFlags;	
	
	//internal to base class state flags - bits in array holding relevant process info restricted to base class
	private int[] _baseStFlags;						
	private static final int
			_BaseDebugIDX 				= 0;		
	private static final int _BaseNumStFlags = 1;	
	
	////////////////////////////////////////
	// visualization functionality
	
	//width of display space - used to shrink things if necessary when right sidebar menu is added
	protected float visScreenWidth;
	
	//current visible screen width and height, from owning window
	protected float[] curVisScrDims;	
	
	public baseProbExpMgr(MessageObject _msgObj, float[] _curVisScrDims) {
		msgObj = _msgObj;
		setVisibleScreenDims(_curVisScrDims);
		//base class-specific flags, isolated to within this code only
		initBaseFlags();
		//init experiment-specific flags
		initFlags();
		//initialize structure to hold settings and options for random variable functions
		randAlgOptSettings = new int[randGenAlgNames.length][];
		for(int i=0;i<randAlgOptSettings.length;++i) {			randAlgOptSettings[i]=new int[] {0};		}//will change settings based on individual random variable functions and UI input
	
		//for display of time since experiment was built 
		Instant now = Instant.now();
		expMgrBuiltTime = now.toEpochMilli();//milliseconds since 1/1/1970 when this exec was built.
		//run once - all solvers for this experiment should be built
		buildSolvers();
		//initialize exp - in instance class
		initExp();
	}//ctor
	/**
	 * Sets the initial uiDataValues UIDataUpdater from owner of experiment. This is used to pass values between experiment and UI. Builds a copy of the passed updater, appropriately cast
	 */
	protected abstract Base_UIDataUpdater initUIDataUpdater(Base_UIDataUpdater dataUpdate);
	
	/**
	 * Sets new values for the UIDataUpdater. Should only be called if values have changed.
	 * @param dataUpdate New data from UI or other data source
	 */
	public final void updateUIDataValues(Base_UIDataUpdater dataUpdate) {
		if(uiDataValues == null) {
			//set the Base_UIDataUpdater - in instance class. Should be a copy of owning window's UIDataUpdater, with all dta from dataUpdate
			uiDataValues = initUIDataUpdater(dataUpdate);
		} else {
			uiDataValues.setAllVals(dataUpdate);
		}
		//instance-class specific handling of new data values
		updateUIDataValues_Priv();
	}
	
	/**
	 * Instance-class specific handling of new UI data values, using this classe's uiDataValues variable
	 */
	protected abstract void updateUIDataValues_Priv();

	/**
	 * (re)init this experiment - specific functionality for each instance class
	 */
	public abstract void initExp();
	/**
	 * update randGen variables with appropriate values, including RVF Opts
	 */
	protected abstract void updateAllRandGens_Priv();
	protected void updateAllRandGens() {		
		updateAllRandGens_Priv();		
	}//
	
	
	/**
	 * build various quadrature solvers to experiment with different formulations
	 */
	private final void buildSolvers() {
		//build quadrature solvers if not built already
		quadSlvrs = new baseQuadrature[numGaussSlvrs];
		quadSlvrs[GL_QuadSlvrIDX] = new myGaussLegenQuad(quadSlvrNames[GL_QuadSlvrIDX],numQuadPoints, quadConvTol, quadBDScale);
		curQuadSolverIDX = 0;
		//instance specific funcitionality
		buildSolvers_Indiv();
	}//buildSolvers
	
	/**
	 * set individual random variable function options from UI - resize array if necessary
	 * _randVarFuncIDX corresponds to static int idxs normRandVarIDX, gaussRandVarIDX, etc specified above
	 * @param _randVarFuncIDX
	 * @param _settingIDX
	 * @param _settingVal
	 */
	private void _setRandVarFuncOpts(int _randVarFuncIDX, int _settingIDX, int _settingVal) {
		if(_settingIDX >= randAlgOptSettings[_randVarFuncIDX].length) {
			int[] newArray =  new int[_settingIDX+1];//large enough to hold new idx being set
			for(int i=0;i<randAlgOptSettings[_randVarFuncIDX].length;++i) {		newArray[i]=randAlgOptSettings[_randVarFuncIDX][i];		}
			randAlgOptSettings[_randVarFuncIDX] = newArray;			
		}
		if(randAlgOptSettings[_randVarFuncIDX][_settingIDX] == _settingVal) {return;}//don't update if same value
		randAlgOptSettings[_randVarFuncIDX][_settingIDX] = _settingVal;	
		//instance-class specific control to update all rvfs
		updateAllRandGens();
	}//setRandVarFuncOpts
	
	//public facing, should have one of these for each random variable function type supported
	public void setNorm_RVFOpts(int _settingIDX, int _settingVal) {_setRandVarFuncOpts(normRandVarIDX, _settingIDX, _settingVal);}
	public void setGauss_RVFOpts(int _settingIDX, int _settingVal) {_setRandVarFuncOpts(gaussRandVarIDX, _settingIDX, _settingVal);}
	public void setRaisedCos_RVFOpts(int _settingIDX, int _settingVal) {_setRandVarFuncOpts(raisedCosRandVarIDX, _settingIDX, _settingVal);}
	public void setCosCDF_RVFOpts(int _settingIDX, int _settingVal) {_setRandVarFuncOpts(cosCDFRandVarIDX, _settingIDX, _settingVal);}
	public void setFleish_RVFOpts(int _settingIDX, int _settingVal) {_setRandVarFuncOpts(fleishRandVarIDX, _settingIDX, _settingVal);}

	
	//this will build the appropriate rand variable function option flags to be used when constructing or modifying the random variables based on UI input
	//different random variable functions will have different options (if any)
	private int[] getRandVarFuncOpts(int _randVarFuncType) {
		if(_randVarFuncType == -1) {return new int[] {0};}
		int[] res = new int[randAlgOptSettings[_randVarFuncType].length];
		System.arraycopy(randAlgOptSettings[_randVarFuncType], 0, res, 0, randAlgOptSettings[_randVarFuncType].length);		
		return res;
	}//buildRandVarFuncOpts
	
	public baseRandVarFunc buildRandVarType (int _pdfType,  int _quadSlvrIdx, myProbSummary_Dbls _summaryObj) {
		//System.out.println("buildRandVarType : " + _pdfType);
		baseRandVarFunc rvf;
		switch (_pdfType) {		
		   	case normRandVarIDX			: { rvf = new myNormalFunc(quadSlvrs[_quadSlvrIdx]); 		   									break;}
	    	case gaussRandVarIDX		: { rvf = new myGaussianFunc(quadSlvrs[_quadSlvrIdx], _summaryObj);	    					break;}
	    	case raisedCosRandVarIDX	: { rvf = new myCosFunc(quadSlvrs[_quadSlvrIdx], _summaryObj);	    						break;}			
	    	case cosCDFRandVarIDX		: { rvf = new myCosFuncFromCDF(quadSlvrs[_quadSlvrIdx], _summaryObj);	    					break;}
	    	case fleishRandVarIDX		: { rvf = new myFleishFunc_Uni(quadSlvrs[_quadSlvrIdx], _summaryObj, randGenAlgNames[fleishRandGen_UniVar]);	    		break;}
			default : {return null;}		
		}
		rvf.setOptionFlags(getRandVarFuncOpts(_pdfType));
		return rvf;
	}//myRandVarFunc
	
	//using default GaussLengendre and 256 zig's for ziggurat alg
	public myRandGen buildAndInitRandGen(int _type, int _pdfType, myProbSummary_Dbls _summaryObj) {return buildAndInitRandGen(_type, _pdfType, GL_QuadSlvrIDX, 256, _summaryObj);}
	//must build rand gen through this method
	public myRandGen buildAndInitRandGen(int _randGenType, int _pdfType, int _quadSlvrIdx, int _numZigRects, myProbSummary_Dbls _summaryObj) {
		
		switch (_randGenType) {
			case boundedRandGen : {
				//need to build a random variable generator function
				baseRandVarFunc func = buildRandVarType(_pdfType, _quadSlvrIdx, _summaryObj);				
				return new myBoundedRandGen(func, randGenAlgNames[_randGenType]);}
			case ziggRandGen : {//ziggurat alg solver - will use zigg algorithm to generate a gaussian of passed momments using a uniform source of RVs
				//need to build a random variable generator function
				baseRandVarFunc func = buildRandVarType(_pdfType, _quadSlvrIdx, _summaryObj);				
				//_numZigRects must be pwr of 2 - is forced to be if is not.  Should be 256
				return new myZigRandGen(func, _numZigRects, randGenAlgNames[_randGenType]);}
			
			case fleishRandGen_UniVar : {
				//specify fleishman rand function with either moments or data - if only moments given, then need to provide hull as well
				baseRandVarFunc func = buildRandVarType(_pdfType, _quadSlvrIdx, _summaryObj);		
				return new myFleishUniVarRandGen(func,  randGenAlgNames[_randGenType]);	}

			//these are just transformations and are not described by an underlying pdf
			case linearTransformMap : {	
				//does not actually represent a random function generator, just represents a transformation, to be used by class rosters in ClassGradeExperiment				
				return new linearTransform(_summaryObj);}
			
			case uniformTransformMap : {
				//does not represent random func generator/model, just represents a mapping from grade to order-based 1/n..n/n equal partition of grade space
				return new uniformCountTransform(_summaryObj);}
			
			default	:	{		
				msgObj.dispWarningMessage("BaseProbExpMgr","buildAndInitRandGen","Unknown random generator type : " + _randGenType + ".  Aborting.");
				return null;
			}
		}//switch
		
	}//buildAndInitRandGen
	
	//set values specific to solver
	public void setSolverVals(int _numPoints, double _tol, int _BDScale) {
		numQuadPoints = _numPoints;
		quadConvTol = _tol;	
		quadBDScale = _BDScale;
		quadSlvrs[curQuadSolverIDX].setSolverVals(numQuadPoints, quadConvTol, quadBDScale);
		//set this solver for all random generators
		
	}//setSolverVals
	
	//called whenever screen width changes due to showing/hiding the right side menu
	public void setVisibleScreenDims(float[] _curVisScrDims) {
		curVisScrDims = _curVisScrDims;
		visScreenWidth = curVisScrDims[0];
		setVisWidth_Priv();
	}//setVisibleScreenWidth
	public float getVisibleSreenWidth() {return visScreenWidth;}

	//set the experiment-specific quantities dependent on visible width of the display area - use this to shrink visualizations if necessary
	protected abstract void setVisWidth_Priv();
	//check mouse over/click in 2d experiment; if btn == -1 then mouse over
	public abstract boolean checkMouseClickInExp2D(int msx, int msy, int btn);
	//check mouse over/click in 2d experiment; if btn == -1 then mouse over
	public abstract boolean checkMouseDragMoveInExp2D(int msx, int msy, int btn);
	//notify all exps that mouse has been released
	public abstract void setMouseReleaseInExp2D();
	//experiment instance-specific solver building funcitionality
	protected abstract void buildSolvers_Indiv();


	//sets which solver will be used for integration for this experiment
	public void setCurSolver(int idx) {
		if(curQuadSolverIDX == idx) {return;}
		curQuadSolverIDX = idx;
		quadSlvrs[curQuadSolverIDX].setSolverVals(numQuadPoints, quadConvTol, quadBDScale);
		//set this solver for all random generators
	}
	
	//conduct a simple test on the passed random number generator - it will get a sample based on the pdf function given to generator
	//uses 64 bit random uniform val
	protected void smplTestRandNumGen(myRandGen gen, int numVals) {
		msgObj.dispInfoMessage("BaseProbExpMgr","testRandGen","Start synthesizing " + numVals+ " values using Gen : \n\t" + gen.getFuncDataStr());
		double[] genVals = new double[numVals];
		for(int i=0;i<genVals.length;++i) {	
			//msgObj.dispMessage("BaseProbExpMgr","testRandGen","Generating val : " + i);
			genVals[i] = gen.getSample();	
		}
		//now calculate mean value and
		msgObj.dispInfoMessage("BaseProbExpMgr","testRandGen","Finished synthesizing " + numVals+ " values using Gen : " + gen.name + " | Begin analysis of values.");
		myProbSummary_Dbls analysis = new myProbSummary_Dbls(genVals);
		msgObj.dispInfoMessage("BaseProbExpMgr","testRandGen","Analysis res of " + gen.name + " : " + analysis.getMomentsVals());
	}//testGen
	
	//conduct a simple test on the passed random number generator - it will get a sample based on the pdf function given to generator
	//uses the "fast" implementation - 32 bit random uniform value
	protected void smplTestFastRandNumGen(myRandGen gen, int numVals) {
		msgObj.dispInfoMessage("BaseProbExpMgr","testRandGen","Start synthesizing " + numVals+ " values using Gen : \n\t" + gen.getFuncDataStr());
		double[] genVals = new double[numVals];
		for(int i=0;i<genVals.length;++i) {	
			//msgObj.dispMessage("BaseProbExpMgr","testRandGen","Generating val : " + i);
			genVals[i] = gen.getSampleFast();	
		}
		//now calculate mean value and
		msgObj.dispInfoMessage("BaseProbExpMgr","testRandGen","Finished synthesizing " + numVals+ " values using Gen : " + gen.name + " | Begin analysis of values.");
		myProbSummary_Dbls analysis = new myProbSummary_Dbls(genVals);
		msgObj.dispInfoMessage("BaseProbExpMgr","testRandGen","Analysis res of " + gen.name + " : " + analysis.getMomentsVals());
	}//testGen
	
	
	//test method for calculating r for rand var function (used in ziggurat algorithm)
	@SuppressWarnings("unused")
	public void testRCalc() {
		msgObj.dispInfoMessage("BaseProbExpMgr","testRCalc","Start test of r var calc");
		baseRandVarFunc randVar = new myNormalFunc(quadSlvrs[GL_QuadSlvrIDX]);
		myProbSummary_Dbls analysis = new myProbSummary_Dbls( new double[] {2.0, 3.0}, 2);
		//myRandVarFunc randGaussVar = new myGaussianFunc(this, quadSlvrs[GL_QuadSlvrIDX], analysis);
		randVar.dbgTestCalcRVal(256);		
		msgObj.dispInfoMessage("BaseProbExpMgr","testRCalc","End test of r var calc");
		
	}//testRCalc
	
	//this will test the zig function sweeping from 0->1 prob (treating as invCDF)
	//we want to make sure that feeding a uniform value in order will generate increasing random variables
	public void testSeqZigGen() {
		double numVals = 200.0;
		msgObj.dispInfoMessage("BaseProbExpMgr","testSeqZigGen","Start synthesizing " + numVals+ " sequential values to test zig gen.");
		myZigRandGen ziggen = new myZigRandGen(new myNormalFunc(quadSlvrs[GL_QuadSlvrIDX]), "Ziggurat Algorithm");
		
		boolean done = false;
		TreeMap<Long, Double> resMap = new TreeMap<Long, Double>();
		for(int i=0;i<numVals;++i) {
			done = false;
			long val = 0L;
			double prob = 0.0, res = 0.0;
			while (!done) {
				prob = (i/numVals) + 1e-10;
				val = (long) (prob * Long.MAX_VALUE);
				res = ziggen.getFuncValFromLong(val);
				done = (res==res);
				if(!done) {
					msgObj.dispInfoMessage("BaseProbExpMgr","testSeqZigGen","i: " + i+ " possible error with : " + prob + " | in : " + val + " | rnd out : " + res +" returns not done." );
					System.exit(0);
				}
			}
			msgObj.dispInfoMessage("BaseProbExpMgr","testSeqZigGen","i: " + i+ " ratio/prob : " + String.format("%.8f",prob) + " | in : " + val + " | rnd out : " + res );
			resMap.put(val,  res);
		}
		msgObj.dispInfoMessage("BaseProbExpMgr","testSeqZigGen","Finished synthesizing " + numVals+ " values using Gen : " + ziggen.name);
	}//testZigGen
	
	//////////////////////////////
	// drawing
	
	public abstract void drawExp();	
	
	
	///////////////////////////
	// end message display functionality

	/////////////////////////////
	// state flags specific to each experiment
	protected abstract void initFlags();
	public abstract void setFlag(int idx, boolean val);
	public final void setAllFlags(int[] idxs, boolean val) {for (int idx : idxs) {setFlag(idx, val);}}
	public final boolean getFlag(int idx){int bitLoc = 1<<(idx%32);return (stFlags[idx/32] & bitLoc) == bitLoc;}
	
	/////////////////////////////
	//init and manage state flags internal to base class computation
	private void initBaseFlags(){_baseStFlags = new int[1 + _BaseNumStFlags/32]; for(int i = 0; i<_BaseNumStFlags; ++i){setBaseFlag(i,false);}}
	@SuppressWarnings("unused")
	private void setAllBaseFlags(int[] idxs, boolean val) {for (int idx : idxs) {setBaseFlag(idx, val);}}
	private void setBaseFlag(int idx, boolean val){
		int flIDX = idx/32, mask = 1<<(idx%32);
		_baseStFlags[flIDX] = (val ?  _baseStFlags[flIDX] | mask : _baseStFlags[flIDX] & ~mask);
		switch (idx) {//special actions for each flag
			case _BaseDebugIDX : {break;}	
			
		}
	}//setFlag		
	@SuppressWarnings("unused")
	private boolean getBaseFlag(int idx){int bitLoc = 1<<(idx%32);return (_baseStFlags[idx/32] & bitLoc) == bitLoc;}		


}//base class for probability experiments

