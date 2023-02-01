package base_ProbTools.randGenFunc;

import base_ProbTools.randGenFunc.gens.base.Base_RandGen;

/**
 * this class holds a description of a random number generator, including the momments and algorithms it uses
 * @author john
 */
public class randGenDesc implements Comparable<randGenDesc>{
	//owning rand gen
	public final Base_RandGen randGen;
	//names of quadrature algorithm and random number distribution name and algorithm name 
	public final String quadName, distName, algName;
	
	public randGenDesc(String _quadName, String _distName, Base_RandGen _randGen) {
		quadName = _quadName; 
		randGen = _randGen;		
		distName = _distName;
		algName = randGen.name;
	}//ctor

	@Override
	public int compareTo(randGenDesc othr) {
		int res = distName.toLowerCase().compareTo(othr.distName.toLowerCase());
		res = (res == 0 ? quadName.toLowerCase().compareTo(othr.quadName.toLowerCase()) : res);
		res = (res == 0 ? algName.toLowerCase().compareTo(othr.algName.toLowerCase()) : res);
		return (res == 0 ? Integer.compare(randGen.ObjID, othr.randGen.ObjID) : res);
	}//compareTo	
	
	@Override
	public String toString() {
		String res = "Alg Name : " + algName + " | Dist : " + distName + " | Quad Alg : " + quadName;
		return res;
	}
}//class RandGenDesc
