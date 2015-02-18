// BootstrappedAlignment.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package pal.alignment;

import pal.math.*;


/**
 * generates bootstrapped alignments from a raw alignment
 *
 * @version $Id: BootstrappedAlignment.java,v 1.6 2003/03/23 00:12:57 matt Exp $
 *
 * @author Korbinian Strimmer
 */
public class BootstrappedAlignment extends AbstractAlignment
{
	//
	// Public stuff
	//
	/**
	 * Constructor
	 *
	 * @param raw original alignment
	 */	
	public BootstrappedAlignment(Alignment raw) {
		this(raw,0);
	}
	/**
	 * Constructor
	 *
	 * @param raw original alignment
	 * @param seed RNG seed
	 */
	public BootstrappedAlignment(Alignment raw, long seed)
	{
		this(raw,seed,0);
	}

	/**
	 * Constructor
	 *
	 * @param raw original alignment
	 * @param seed RNG seed
	 * @param nSites Number of sites to include in bootstrapped alignment. 
	 *        (0=same length as original alignment.)
	 */
	public BootstrappedAlignment(Alignment raw, long seed, int nSites)
	{
		rawAlignment = raw;

		numSeqs = raw.getSequenceCount();
		idGroup = raw;
		numSites = (nSites==0) ? raw.getSiteCount() : nSites;
		setDataType(raw.getDataType());

		alias = new int[numSites];
		urn = new UrnModel(numSites,seed);

		bootstrap();
	}

	// Implementation of abstract Alignment method

	/** sequence alignment at (sequence, site) */
	public char getData(int seq, int site)
	{
		return rawAlignment.getData(seq, alias[site]);
	}


	/** bootstrap alignment */
	public void bootstrap()
	{
		for (int i = 0; i < numSites; i++)
		{
			alias[i] = urn.drawPutBack();
		}
	}


	//
	// Private stuff
	//

	private UrnModel urn;
	private Alignment rawAlignment;
	private int[] alias;


}
