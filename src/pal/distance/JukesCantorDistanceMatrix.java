// JukesCantorDistanceMatrix.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

/**
 * Comment by Michael Woodhams:
 * I found a **HUGELY** significant bug in here.
 * In computeDistances, const1 and const2 were set from *THE NUMBER OF TAXA*
 * when it should have been the number of *STATES*.
 * I've fixed this.
 * However, in jccorrection, there is an upper limit of  BranchLimits.MAXARC
 * on the distance, which defaults to 1.0, which is way too low.
 * As a result, I've lost confidence in this code. Testing must have
 * been non-existent.
 */
package pal.distance;

import pal.alignment.*;
import pal.misc.*;


/**
 * compute jukes-cantor corrected distance matrix
 *
 * @version $Id: JukesCantorDistanceMatrix.java,v 1.5 2002/12/05 04:27:28 matt Exp $
 *
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 */
public class JukesCantorDistanceMatrix extends DistanceMatrix
{
	//
	// Public stuff
	//

	/**
	 * compute jukes-cantor corrected distances
	 * (assumes nucleotides as underlying data)
	 *
	 * @param dist distance matrix
	 */
	public JukesCantorDistanceMatrix(DistanceMatrix dist)
	{
		this(dist, 4);
	}


	/**
	 * compute jukes-cantor corrected distances
	 *
	 * @param dist distance matrix
	 * @param numStates number of states of underlying data
	 */
	public JukesCantorDistanceMatrix(DistanceMatrix dist, int numStates)
	{
		super(computeDistances(dist, numStates), dist);

	}


	/**
	 * compute jukes-cantor corrected distances
	 *
	 * @param alignment Alignment
	 */
	public JukesCantorDistanceMatrix(Alignment alignment)
	{
		this(new SitePattern(alignment));
	}

	/**
	 * compute jukes-cantor corrected distances
	 *
	 * @param sitePattern SitePattern
	 */
	public JukesCantorDistanceMatrix(SitePattern sitePattern)
	{
		this(	new AlignmentDistanceMatrix(sitePattern),
			sitePattern.getDataType().getNumStates());
	}

	private static final double[][] computeDistances(final DistanceMatrix dist, final int numberOfStates)
	{
		final int numSeqs = dist.getSize();
		final double[][] distance = new double[numSeqs][numSeqs];
		final double[][] obsDistance = dist.getDistances();
		// MDW: Huge bug fix: used to have numSeqs instead of numberOfStates for these two statements
		final double const1 = (numberOfStates-1)/(double)numberOfStates;
		final double const2 = numberOfStates/(double)(numberOfStates-1);
		for (int i = 0; i < numSeqs-1; i++)
		{
			distance[i][i] = 0.0;
			for (int j = i+1; j < numSeqs; j++)
			{
				distance[i][j] = distance[j][i] = jccorrection(const1, const2, obsDistance[i][j]);
			}
		}
		return distance;
	}


	private static final double jccorrection(final double const1,  final double const2, double obsdist)
	{
		if (obsdist == 0.0) return 0.0;

		if (obsdist >= const1)
		{
			return BranchLimits.MAXARC;
		}

		double expDist = -const1 * Math.log(1.0 - (const2 * obsdist));

		if (expDist < BranchLimits.MAXARC)
		{
			return expDist;
		}
		else
		{
			return BranchLimits.MAXARC;
		}
	}
}
