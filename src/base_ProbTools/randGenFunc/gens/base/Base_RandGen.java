package base_ProbTools.randGenFunc.gens.base;

import java.util.concurrent.ThreadLocalRandom;

import base_ProbTools.randGenFunc.randGenDesc;
import base_ProbTools.randGenFunc.funcs.base.Base_RandVarFunc;
import base_ProbTools.randGenFunc.gens.myBoundedRandGen;
import base_StatsTools.summary.myProbSummary_Dbls;
import base_StatsTools.visualization.myDistFuncHistVisMgr;

/**
 * Provides generation of random variables from prob distributions given a uniform distribution
 */
public abstract class Base_RandGen implements Comparable<Base_RandGen> {
	public final int ObjID;
	private static int IDCnt = 0;
	//original data and analysis of it - fl polynomial needs to be built from a sample distribution or from a set of moments
	protected myProbSummary_Dbls summary;
    //random generator to use to generate uniform data - threadsafe
	public final String name;	
	//function name for this randgen
	protected String funcName;
	//descriptor of this random generator
	public randGenDesc desc;		
	//function this rand gen uses
	protected final Base_RandVarFunc func;	
	//visualization tool for this random generator
	//protected myDistFuncHistVisMgr distVisObj; 
	
	//state flags - bits in array holding relevant info about this random variable function
	protected RandGenStateFlags stFlags;
   
	public Base_RandGen(Base_RandVarFunc _func, String _name) {
		ObjID = IDCnt++;  
		name=_name;
		stFlags = new RandGenStateFlags(this);
		func = _func;
		initRandGen();
    }//ctor
	//overriden by transforms
	protected void initRandGen() {
		stFlags.setUseInSolver(true);
		funcName= func.name; 
		desc = new randGenDesc(func.getQuadSolverName(), funcName, this);
		//func built with summary data - allow for quick access
		summary = func.getSummary();	
		 _setFuncSummaryIndiv();
	}
	
	//set summary for this object and for function
	public void setFuncSummary(myProbSummary_Dbls _summary) {
		summary = _summary;	
		func.rebuildFuncs(_summary);
		 _setFuncSummaryIndiv();
	}//setFuncSummary
	
	//when new options are specified, rebuild functions as if new summary was specified
	public void setOptionFlags(int[][] _opts) {
		func.setOptionFlags( _opts[func.getRVFType()]);
		func.rebuildFuncs(summary);
		 _setFuncSummaryIndiv();
	}
	
	//called whenever summary object is set/reset
	public abstract void _setFuncSummaryIndiv();	

    //thread-safe queries for uniform values
    protected long getNextLong() {return ThreadLocalRandom.current().nextLong();}  
    protected int getNextInt() {return ThreadLocalRandom.current().nextInt();  }
    protected double getNextDouble() {return ThreadLocalRandom.current().nextDouble();}
    protected double getUniform01() {
    	long val = ThreadLocalRandom.current().nextLong();
    	return .5+ .5 * val/Long.MAX_VALUE;
    }
    /**
     * uniformly between [min,max)
     * @param min
     * @param max
     * @return
     */
    protected int getUniInt(int min, int max) {    	
    	return ThreadLocalRandom.current().nextInt(min,max);
    }
	
    public myProbSummary_Dbls getSummary() {return summary;}
    public Base_RandVarFunc getFunc() {return func;}
	public double getMean() {return summary.mean();}
	public double getStd() {return summary.std();}
	public double getVar() {return summary.var();}
	public double getSkew() {return summary.skew();}
	public double getKurt() {return summary.kurt();}    
    
    public abstract double[] getMultiSamples(int num);
    public abstract double[] getMultiFastSamples(int num);
    /**
     * return a sample based on func  - momments defined by myRandVarFunc
     * @return
     */
	public abstract double getSample();
	public abstract double getSampleFast();
	/**
	 * mapping to go from distribution to uniform 0->1 (from val -> prob p(X<= val))
	 * @param _val
	 * @return
	 */
	public abstract double CDF(double _val);
	/**
	 * mapping to go from uniform 0->1 to distribution (from p(X<= val) -> val) 
	 * @param _val
	 * @return
	 */
	public abstract double inverseCDF(double _val);
	/**
	 * alias for CDF
	 * @param _val
	 * @return
	 */
	public double distToUniform(double _val) {return CDF(_val);}
	/**
	 * alias for inverseCDF
	 * @param _val
	 * @return
	 */
	public double uniformToDist(double _val) {return inverseCDF(_val);}
	/**
	 * test integral evaluation
	 * @param min
	 * @param max
	 * @return
	 */
	public double testInteg(double min, double max) {		return func.integral_f(min, max);	}
	
	@Override
	public int compareTo(Base_RandGen othr) {return desc.compareTo(othr.desc);}
	
	
	//private final String[] dispMultiStrsConst = new String[] {"PDF hist",}; 
	/**
	 * build dist and hist and also take passed cosine randgen and superimpose the values for its pdf
	 * @param numVals
	 * @param numBuckets
	 * @param low
	 * @param high
	 * @param cosGen
	 */
	public void buildFuncHistCosPlot(myDistFuncHistVisMgr distMdlViz,int numVals, int numBuckets, double low, double high, myBoundedRandGen cosGen) {
		//first build histogram
		calcHistValsForDisp(distMdlViz, numVals, numBuckets);
		//build and set pdf function values
		func.buildFuncPlotVals(distMdlViz, numVals, low, high, Base_RandVarFunc.queryPDFIDX);
		//now use passed cosGen but populate it into this object's distVisObj
		cosGen.func.buildFuncPlotVals(distMdlViz, numVals, low, high, Base_RandVarFunc.queryPDFIDX);
		
		String histName = funcName+" PDF hist",
				gaussName = func.getDispFuncName(Base_RandVarFunc.queryPDFIDX),
				cosName = cosGen.func.getDispFuncName(Base_RandVarFunc.queryPDFIDX);
		
		//get min and max histogram values and get min/max/diff y values for larger of two dists, either cosine or gauss
		double[][] minMaxDiffHist = distMdlViz.getSpecificMinMaxDiff(histName),//use this for x values		
				minMaxDiffCos = distMdlViz.getSpecificMinMaxDiff(cosName),
				minMaxDiffGauss = distMdlViz.getSpecificMinMaxDiff(gaussName);
		double[][] minMaxDiff = new double[2][];
		minMaxDiff[0] = minMaxDiffHist[0];
		minMaxDiff[1] = new double[3];
		minMaxDiff[1][0] = (minMaxDiffCos[1][0] < minMaxDiffGauss[1][0]) ? minMaxDiffCos[1][0] : minMaxDiffGauss[1][0];
		minMaxDiff[1][1] = (minMaxDiffCos[1][1] > minMaxDiffGauss[1][1]) ? minMaxDiffCos[1][1] : minMaxDiffGauss[1][1];
		minMaxDiff[1][2] = minMaxDiff[1][1] - minMaxDiff[1][0];
		String[] dispMultiStrs = new String[] {histName, gaussName, cosName};
		distMdlViz.setCurMultiDispVis(dispMultiStrs,minMaxDiff);
		distMdlViz.setColorVals(cosName,"stroke", new int[] {255,255,0,255});
	}//buildFuncHistCosPlot
	
	
	/**
	 * synthesize numVals values from low to high to display 
	 * @param numVals
	 * @param low
	 * @param high
	 * @param funcType
	 */
	public void calcFuncValsForDisp(myDistFuncHistVisMgr distMdlViz, int numVals,double low, double high,  int funcType ) {
		func.buildFuncPlotVals(distMdlViz, numVals, low, high, funcType);
	}//calcFValsForDisp
		
	/**
	 * build display function for histogram Num buckets should be << numVals
	 * @param numVals
	 * @param numBuckets
	 */
	public void calcHistValsForDisp(myDistFuncHistVisMgr distMdlViz, int numVals, int numBuckets) {
		//x val is distribution/max bucket value, y val is count
		double [] histVals = getMultiSamples(numVals);		
		myProbSummary_Dbls summary = new myProbSummary_Dbls(histVals);
		//build buckets : numBuckets+1 x 2 array; 2nd idxs : idx 0 is lower x value of bucket, y value is count; last entry should always have 0 count
		double[][] distBuckets = summary.calcBucketVals(numBuckets);
		System.out.println("calcHistValsForDisp : numVals : " + numVals + " | size of histVals :"+histVals.length + " | size of distBuckets : "+ distBuckets.length + " | distVisObj is null :  "+ (null==distMdlViz));
		
		//min, max and diff values for x axis (rand val) and y axis (counts)
		double[][] minMaxDiffXVals = new double[2][3];
		minMaxDiffXVals[0][0] = summary.getMin();
		minMaxDiffXVals[0][1] = summary.getMax();
		minMaxDiffXVals[0][2] = minMaxDiffXVals[0][1] - minMaxDiffXVals[0][0];
		//bucket count min max diff - y axis
		minMaxDiffXVals[1][0] = 100000;
		minMaxDiffXVals[1][1] = -100000;
		for(int i=0;i<distBuckets.length;++i) {
			minMaxDiffXVals[1][0] = (minMaxDiffXVals[1][0] > distBuckets[i][1] ? distBuckets[i][1] : minMaxDiffXVals[1][0]);
			minMaxDiffXVals[1][1] = (minMaxDiffXVals[1][1] < distBuckets[i][1] ? distBuckets[i][1] : minMaxDiffXVals[1][1]);
		}		
		minMaxDiffXVals[1][2] = minMaxDiffXVals[1][1] - minMaxDiffXVals[1][0];
		
		distMdlViz.setValuesHist(funcName+" PDF hist", new int[][] {new int[] {255,0,0,255}, new int[] {255,255,255,255}}, distBuckets, minMaxDiffXVals);
	}//calcDistValsForDisp
		
	/**
	 * Debug mode functionality. Called from flags structure
	 * @param val
	 */
	public final void handleRandGenStateFlagsDebugMode(boolean val) {
		//TODO	
	}
	
	/**
	 * return string description of rand function
	 * @return
	 */
	public String getFuncDataStr() {return func.getFuncDataStr();}
	
	/**
	 * get short string suitable for key for map
	 * @return
	 */
	public String getTransformName() {
		String res = name+"_"+ desc.quadName+"_" + func.getMinDescString();
		return res;
	}
	
	/**
	 * get short display string
	 * @return
	 */
	public String getDispTransName() {
		return name +" "+summary.getMinNumMmntsDesc();
	}

}//class myRandGen

////////////////////////////////////////////////////////////////////////////////////////
// child classes 




