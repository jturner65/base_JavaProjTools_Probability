package base_ProbTools.randGenFunc.funcs.base;

import base_Utils_Objects.tools.flags.Base_BoolFlags;

/**
 * @author John Turner
 *
 */
public class RandVarFuncStateFlags extends Base_BoolFlags {
	/**
	 * Owning random variable function
	 */
	private final Base_RandVarFunc owner;
	public static final int
		//whether to use zig alg for solving	
		useZigAlgIDX				= _numBaseFlags,		//whether or not we're using ziggurat solver		
		//quad solver set
		quadSlvrSetIDX				= _numBaseFlags+1;		//whether quadrature solver has been set or not
	private static final int _numStateFlags =_numBaseFlags +2;	
	
	/**
	 * @param _numFlags
	 */
	public RandVarFuncStateFlags(Base_RandVarFunc _owner) {
		super(_numStateFlags);
		owner = _owner;
	}
	
	/**
	 * Get whether or not this random variable will be used in a ziggurat solver		
	 * @return
	 */
	public final boolean getUseZigSolver() {return getFlag(useZigAlgIDX);}
	
	/**
	 * Set whether or not this random variable will be used in a ziggurat solver		
	 * @param val
	 */
	public void setUseZigSolver(boolean val) {
		setFlag(useZigAlgIDX, val);
	}
	
	/**
	 * Get whether quadrature solver has been set	
	 * @return
	 */
	public final boolean getQuadSolverSet() {return getFlag(quadSlvrSetIDX);}
	
	/**
	 * Set whether quadrature solver has been set		
	 * @param val
	 */
	public void setQuadSolverSet(boolean val) {
		setFlag(quadSlvrSetIDX, val);
	}
	
		
	/**
	 * Set or clear debug functionality for flag owner
	 */
	@Override
	protected void handleSettingDebug(boolean val) {owner.handleDebugMode(val);		}

	@Override
	protected void handleFlagSet_Indiv(int idx, boolean val) {
		switch (idx) {//special actions for each flag		
			case useZigAlgIDX 	: {break;}
			case quadSlvrSetIDX	: {break;}
		}	
	}//handleFlagSet_Indiv

}//class RandVarFuncStateFlags
