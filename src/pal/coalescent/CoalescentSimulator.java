// CoalescentSimulator.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package pal.coalescent;

import pal.tree.*;

/**
 * Simulates a set of coalescent intervals given a demographic model.
 *
 * @version $Id: CoalescentSimulator.java,v 1.5 2001/07/12 12:17:43 korbinian Exp $
 *
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 */
public class CoalescentSimulator {

	/**
	 * Simulates a set of CoalescentIntervals from a genealogy assuming
	 * contemporaneous tips.
	 * @param numLines the number of tips in the sample genealogy
	 * @param model the demographic model to use
	 */
	public CoalescentIntervals simulateIntervals(int numLines, DemographicModel model) {
	
		CoalescentIntervals ci = new CoalescentIntervals(numLines-1);
	
		double currentTime = 0.0;
		for (int i = 0; i < (numLines - 1); i++) {
			//try {
				ci.setInterval(i, 
					model.getSimulatedInterval(numLines,
					currentTime));
			//} catch (CoalescentException ce) {
			//	ce.printStackTrace();
			//}
			ci.setNumLineages(i, numLines);
			currentTime += ci.getInterval(i);
			numLines -= 1;
		}
		return ci;
	}
}

