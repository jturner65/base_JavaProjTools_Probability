package base_ProbTools.samples;

import java.util.*;

import base_ProbTools.randGenFunc.gens.myBoundedRandGen;
import base_ProbTools.randGenFunc.gens.base.myRandGen;
import base_StatsTools.visualization.myDistFuncHistVisMgr;
import base_Utils_Objects.io.messaging.MessageObject;

/**
 * a sample of multiple observations from a distribution
 * @author john
 *
 */
public abstract class mySampleSet implements Comparable<mySampleSet> {
	//experiment owning/using this sample set
	public static MessageObject msgObj;
	public final int ObjID;
	private static int IDCnt = 0;
	//sample set name
	public final String name;
	/**
	 * rand gen used to model underlying grade distribution for this class
	 */
	protected HashMap<String, myRandGen> baseDistModels;
	/**
	 * map keyed by same value holding visualization tools
	 */
	protected HashMap<String, myDistFuncHistVisMgr> distModelVis;
	//currently used dist model
	protected String curDistModel;

	public mySampleSet( String _name) {
		msgObj = MessageObject.buildMe();
		ObjID = IDCnt++;  name=_name;	
		curDistModel = "";
		baseDistModels = new HashMap<String, myRandGen>();
		distModelVis = new HashMap<String, myDistFuncHistVisMgr>();
	}//ctor
	
	//when new transform added, need to clear out existing transformed grades
	public void setBaseDistModel(myRandGen _randGen) {
		curDistModel = _randGen.name;
		baseDistModels.put(curDistModel, _randGen);
		distModelVis.put(curDistModel, buildVisMgr(_randGen.name));
		setBaseDistModel_Indiv();
	}//setRandGenAndType
	
	/**
	 * Build visualization mgr corresponding to current randGen model
	 * @return
	 */
	protected abstract myDistFuncHistVisMgr buildVisMgr(String _name);
	
	
	public void setCurDistModel(String desMdlName) {
		if(null==baseDistModels.get(desMdlName)) {
			msgObj.dispWarningMessage("myClassRoster", "setCurDistModel", "Desired base dist model : " + desMdlName+" has not been set/is null.  Aborting");	return;
		}
		curDistModel = desMdlName;
	}
	public String getCurDistModel() {return curDistModel;}
	//instance class specific functionality for setting base distribution model
	protected abstract void setBaseDistModel_Indiv();
	
	////////////////////////////////////////
	// underlying distribution config, evaluation and plotting functions
	public void setRVFOptionFlags(int[][] _opts) {
		myRandGen baseDistModel = baseDistModels.get(curDistModel);
		if(baseDistModel != null) {baseDistModel.setOptionFlags(_opts);}}
	
	/**
	 * This will evaluate the cosine and the gaussian functions against a histogram, showing the performance of these functions when built from a histogram data
	 * @param numVals
	 * @param numBuckets
	 * @param low
	 * @param high
	 */
	public void evalCosAndNormWithHist(int numVals, int numBuckets, double low, double high) {
		//we wish to build a histogram of current gaussian distribution, then we wish to superimpose the gaussian pdf curve over the histogram, and then superimpose the cosine pdf curve
		myRandGen baseDistModel = baseDistModels.get(curDistModel);
		myDistFuncHistVisMgr distMdlViz = distModelVis.get(curDistModel);
		myBoundedRandGen cosGen = (myBoundedRandGen) baseDistModels.get("Bounded PDF Algorithm");
		baseDistModel.buildFuncHistCosPlot(distMdlViz, numVals, numBuckets, low, high, cosGen);
		
	}//evalCosAndNormWithHist
	
	
	public void evalAndPlotFuncRes(int numVals, double low, double high, int funcType) {
		myRandGen baseDistModel = baseDistModels.get(curDistModel);
		myDistFuncHistVisMgr distMdlViz = distModelVis.get(curDistModel);
		if(baseDistModel == null) {			
			msgObj.dispWarningMessage("myClassRoster", "evalAndPlotFuncRes", "curDistModel has not been set/is null.  Aborting");	
			return;	}
		baseDistModel.calcFuncValsForDisp(distMdlViz, numVals, low, high, funcType);		
	}
	public void evalAndPlotHistRes(int numVals, int numBuckets) {
		myRandGen baseDistModel = baseDistModels.get(curDistModel);
		myDistFuncHistVisMgr distMdlViz = distModelVis.get(curDistModel);
		if(baseDistModel == null) {			
			msgObj.dispWarningMessage("myClassRoster", "evalAndPlotHistRes", "curDistModel has not been set/is null.  Aborting");	
			return;	}		
		baseDistModel.calcHistValsForDisp(distMdlViz, numVals, numBuckets);
	}
	
	public void clearPlotEval() {	distModelVis.get(curDistModel).clearEvalVals();	}	
	//draw plot results from functional histogram/evaluation of baseDistModel
	public void drawPlotRes() {		distModelVis.get(curDistModel).drawVis();	}	
	
	//incase we wish to store class sample sets in sorted mechanism
	@Override
	public int compareTo(mySampleSet othr) {
		int res = this.name.toLowerCase().compareTo(othr.name.toLowerCase());
		return (res == 0 ? Integer.compare(ObjID, othr.ObjID) : res);
	}//compareTo
	
	
	
}//class mySampleSet


