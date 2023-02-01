package base_ProbTools.randGenFunc.transform.base;

import base_ProbTools.randGenFunc.randGenDesc;
import base_ProbTools.randGenFunc.gens.base.Base_RandGen;
import base_StatsTools.summary.myProbSummary_Dbls;

//////////////////////////////////
//linear and uniform transformation classes
//	these are just mappers and will not be used to synthesize random values

public abstract class baseTransform extends Base_RandGen{
	//func will be null for these, so all functionality that is dependent on func variable needs to be overridden
	public baseTransform(String _name, myProbSummary_Dbls _summary) {
		super(null, _name);
		setFuncSummary(_summary);
	}
	//overrding base class verison to remove refs to func
	@Override
	protected final void initRandGen() {
		stFlags.setUseInSolver(func!=null);
		funcName=  "No Rand Func for Transform " + name; 
		desc = new randGenDesc("No Quad Solver", "No Rand Func", this);
		//distVisObj = null;
	}//initRandGen
	
	//override base class version to remove ref to func, which will be null
	@Override
	public void setFuncSummary(myProbSummary_Dbls _summary) {
		summary = _summary;	
		 _setFuncSummaryIndiv();
	}//setFuncSummary
	
	//for a transform this is meaningless - transforms just remap given data to affine transformations, they don't model them
	public void calcDistValsForDisp(int numVals, int numBuckets) {}
	public void calcFuncValsForDisp(int numVals, double low, double high, int funcType ) {}
	
	//return string description of rand function
	@Override
	public String getFuncDataStr() {return "No Function for Transform RandGen - only has mapping";}

	@Override
	public String getTransformName() {		return name+"_"+ _getTransformNameIndiv();	}
	
	public abstract String _getTransformNameIndiv();


}//class transform



