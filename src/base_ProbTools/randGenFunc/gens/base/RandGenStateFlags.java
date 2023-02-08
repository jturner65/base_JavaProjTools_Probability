/**
 * 
 */
package base_ProbTools.randGenFunc.gens.base;

import base_Utils_Objects.tools.flags.Base_BoolFlags;

/**
 * @author John Turner
 *
 */
public class RandGenStateFlags extends Base_BoolFlags {
	/**
	 * owning myRandGen 
	 */
	private final Base_RandGen owner;
	
	public static final int		
		funcSetIDX						= _numBaseFlags;		//whether or not this random variable will be used in a solver		
	private static final int _numStateFlags =_numBaseFlags +1;	

	/**
	 * @param _numFlags
	 */
	public RandGenStateFlags(Base_RandGen _owner) {
		super(_numStateFlags);
		owner = _owner;
	}
	/**
	 * Get whether or not this random variable will be used in a solver		
	 * @return
	 */
	public final boolean getUseInSolver() {return getFlag(funcSetIDX);}
	
	/**
	 * Set whether or not this random variable will be used in a solver		
	 * @param val
	 */
	public void setUseInSolver(boolean val) {
		setFlag(funcSetIDX, val);
	}
	
	@Override
	protected void handleFlagSet_Indiv(int idx, boolean val, boolean oldVal) {		
		switch (idx) {//special actions for each flag
			case debugIDX : 		{break;}	
			case funcSetIDX : 		{break;}	
		}
	}

	@Override
	protected void handleSettingDebug(boolean val) {
		owner.handleRandGenStateFlagsDebugMode(val);		
	}

}//class RandGenStateFlags
