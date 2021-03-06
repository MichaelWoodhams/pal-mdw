// ConcatenatedAlignment.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package pal.alignment;

/**
 * concatenates a list of alignments to one single alignment,
 * increasing the number of sites
 *
 * @version $Id: ConcatenatedAlignment.java,v 1.3 2001/07/13 14:39:12 korbinian Exp $
 *
 * @author Korbinian Strimmer
 */
public class ConcatenatedAlignment extends AbstractAlignment
{
	//
	// Public stuff
	//
	
	/**
	 * concatenate alignments
	 *
	 * @param Alignment array with alignment to concatenate
	 */
	public ConcatenatedAlignment(Alignment[] list)
		throws IllegalArgumentException
	{
		alignmentList = list;
		
		numAlignments = alignmentList.length;
		if (numAlignments == 0)
		{
			throw new IllegalArgumentException("NO ALIGNMENT");
		} 
		
		numSeqs = alignmentList[0].getSequenceCount();
		idGroup = alignmentList[0];
		dataType = alignmentList[0].getDataType(); 
		
		numSites = 0;
		for (int i = 0; i < numAlignments; i++)
		{
			numSites += alignmentList[i].getSiteCount();
			
			if (alignmentList[i].getSequenceCount() != numSeqs || alignmentList[i].getDataType() != dataType)
			{
				throw new IllegalArgumentException("INCOMPATIBLE ALIGNMENTS");
			}
		}
		
		// Create indices
		alignmentIndex = new int[numSites];
		siteIndex = new int[numSites];
		
		int s = 0;
		for (int i = 0; i < numAlignments; i++)
		{
			for (int j = 0; j < alignmentList[i].getSiteCount(); j++)
			{
				alignmentIndex[s+j] = i;
				siteIndex[s+j] = j;
			}
			s += alignmentList[i].getSiteCount();
		}		
	}

	// Implementation of abstract Alignment method

	/** sequence alignment at (sequence, site) */
	public char getData(int seq, int site)
	{
		return alignmentList[alignmentIndex[site]].getData(seq, siteIndex[site]);
	}
	

	//
	// Private stuff
	//
	
	private Alignment[] alignmentList;
	private int numAlignments;
	private int[] alignmentIndex;
	private int[] siteIndex;
}
