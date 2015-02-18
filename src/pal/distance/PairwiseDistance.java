// PairwiseDistance.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package pal.distance;

import pal.alignment.*;
import pal.datatype.*;
import pal.math.*;
import pal.misc.*;
import pal.substmodel.*;


/**
 * Given an alignment (as a SitePattern), calculates any the pairwise distance.
 * Optionally the distance is based on a SubstitutionModel (corrected for multiple mutations etc.)
 * or you can take the 'observed' distance, which is the proportion of sites which do not match.
 * Distances are not cached - they are recalculated every time.
 * 
 *
 * @version $Id: PairwiseDistance.java,v 1.10 2002/09/08 03:47:01 matt Exp $
 * 2009/03/19, v1.20 by Michael Woodhams - added getLogDetDistance methods.
 * Std Error of LogDet distance taken from Massingham and Goldman, MBE 24 2277 (2007)
 * Also added divergence matrix as public data member - this is the divergence matrix
 * for the last estimated distance, like distanceSE is.
 * 
 * The definition of LogDet distance used here includes a 1/N factor (where N is the number of states).
 * This factor is included in Lockhart et al. MBE 11 605 (1994) and Massingham and Goldman but not in
 * Lake, PNAS 91 1455 (1994)
 *
 * @author Korbinian Strimmer
 */
public class PairwiseDistance implements java.io.Serializable
{
	//
	// Public stuff
	//

	/** last estimated distance */
	public double distance;

	/** last estimate standard error of a distance */
	public double distanceSE;
	
	//
	// Private stuff
	//

	private int numSites;
	private int numPatterns;
	private int numStates;
	private int[] weight;
	private double jcratio;
	private boolean modelBased;
	private SitePattern sitePattern;
	private SubstitutionModel model;
	private UnivariateMinimum um;
	private SequencePairLikelihood of;
	private DataType patternDataType_;
	/** Divergence matrix for last estimated distance */
	private IntMatrix divergence;


	//
	// Serialization code
	//
	private static final long serialVersionUID=4721003948793106145L;

	//serialver -classpath ./classes pal.distance.PairwiseDistance
	private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
		out.writeByte(1); //Version number
		out.writeDouble(distance);
		out.writeDouble(distanceSE);
		out.writeObject(divergence);
		out.writeInt(numSites);
		out.writeInt(numPatterns);
		out.writeInt(numStates);
		out.writeObject(weight);
		out.writeDouble(jcratio);
		out.writeBoolean(modelBased);
		out.writeObject(model);
		out.writeObject(sitePattern);
	}

	private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException{
		byte version = in.readByte();
		switch(version) {
			default : {
				distance  = in.readDouble();
				distanceSE = in.readDouble();
				divergence = (IntMatrix)in.readObject();
				numSites = in.readInt();
				numPatterns = in.readInt();
				numStates  =in.readInt();
				weight = (int[])in.readObject();
				jcratio =  in.readDouble();
				//is modelBased really required?
				modelBased = in.readBoolean();
				model = (SubstitutionModel)in.readObject();
				sitePattern = (SitePattern)in.readObject();

				um = new UnivariateMinimum();
				if(model!=null) {
					of = new SequencePairLikelihood(sitePattern, model);
				}
				patternDataType_ = sitePattern.getDataType();
				break;
			}
		}
	}


	/**
	 * Constructor 1 (estimate observed distances only)
	 *
	 * @param sp site pattern
	 */
	public PairwiseDistance(SitePattern sp)
	{
		updateSitePattern(sp);
		divergence = new IntMatrix(numStates);
	}

	/**
	 * Constructor 2 (uses evolutionary model)
	 *
	 * @param sp site pattern
	 * @param m evolutionary model
	 */
	public PairwiseDistance(SitePattern sp, SubstitutionModel m)
	{
		this(sp);

		modelBased = true;
		of = new SequencePairLikelihood(sp, m);
		um = new UnivariateMinimum();
	}

	/**
	 * update model of substitution
	 *
	 * @param model of substitution
	 */
	public void updateModel(SubstitutionModel m)
	{
		if(of==null) {
			modelBased = true;
			of = new SequencePairLikelihood(sitePattern, m);
			um = new UnivariateMinimum();
			this.model = m;
		} else {
			of.updateModel(m);
		}
	}

	/**
	 * update site pattern
	 *
	 * @param site pattern
	 */
	public void updateSitePattern(SitePattern sp)
	{
		sitePattern = sp;
		numSites = sp.getSiteCount();
		numPatterns = sp.numPatterns;
		numStates = sp.getDataType().getNumStates();
		weight = sp.weight;

		jcratio = ((double) numStates-1.0)/(double) numStates;

		if (modelBased) of.updateSitePattern(sp);
		patternDataType_ = sp.getDataType();

	}

	/**
	 * compute distance between two sequences in the given alignment
	 *
	 * @param s1 number of first sequence
	 * @param s2 number of second sequence
	 *
	 * @return estimated distance (observed or ML, depending on constructor used)
	 */
	public double getDistance(int s1, int s2)
	{
		return getDistance(sitePattern.pattern[s1], sitePattern.pattern[s2]);
	}

	/**
	 * compute distance between two sequences (not necessarily
	 * in the given alignment but with the same weights in the site pattern)
	 *
	 * @param s1 site pattern of first sequence
	 * @param s2 site pattern of second sequence
	 *
	 * @return estimated distance (observed or ML, depending on constructor used)
	 */
	public double getDistance(byte[] s1, byte[] s2)
	{
		double dist = getObservedDistance(s1, s2);

		if (modelBased && dist != 0.0)
		{
			// Apply generalized JC correction if possible
			double start = 1.0 - dist/jcratio;
			if (start > 0.0)
			{
				start = -jcratio*Math.log(start);
			}
			else
			{
				start = dist;
			}

			// Determine ML distance
			of.setSequences(s1, s2);
			if (start > BranchLimits.MAXARC || start < BranchLimits.MINARC)
			{
				// Don't use start value
				dist = um.findMinimum(of, BranchLimits.FRACDIGITS);
			}
			else
			{
				// Use start value
				dist = um.findMinimum(start, of, BranchLimits.FRACDIGITS);
			}
		}

		if (modelBased)
		{
			double f2x = um.f2minx;

			if (1.0/(BranchLimits.MAXARC*BranchLimits.MAXARC) < f2x)
			{
				distanceSE = Math.sqrt(1.0/f2x);
			}
			else
			{
				distanceSE = BranchLimits.MAXARC;
			}
		}
		else
		{
			distanceSE = 0.0;
		}

		distance = dist;

		return dist;
	}


	private boolean isDifferent(int s1, int s2)
	{
		// ? is considered identical to anything
		if (patternDataType_.isUnknownState(s1) || patternDataType_.isUnknownState(s2))
		{
			return false;
		}

		// Check for identity
		return (s1 != s2);
	}

	private double getObservedDistance(byte[] seqPat1, byte[] seqPat2)
	{
		int diff = 0;
		for (int i = 0; i < numPatterns; i++)
		{
			if (isDifferent(seqPat1[i], seqPat2[i]))
			{
				diff += weight[i];
			}
		}

		return (double) diff/(double) numSites;
	}
	
	/**
	 * compute LogDet distance between two sequences in the given alignment
	 *
	 * @param s1 number of first sequence
	 * @param s2 number of second sequence
	 *
	 * @return estimated distance 
	 * 
	 * Also calculates distance standard error, but this has not been tested - use at your own risk.
	 */

	public double getLogDetDistance(int s1, int s2)
	{
		return getLogDetDistance(sitePattern.pattern[s1], sitePattern.pattern[s2]);
	}
		
	/**
	 * compute LogDet distance between two sequences (not necessarily
	 * in the given alignment but with the same weights in the site pattern)
	 * Sites where one or other sequence has a gap are simply ignored.
	 *
	 * @param s1 site pattern of first sequence
	 * @param s2 site pattern of second sequence
	 *
	 * @return estimated distance 
	 * 
	 * Sets 'divergence' data member, which can then be accessed by 'getDivergence()'. 
	 * 
	 * Also calculates distance standard error, but this has not been tested - use at your own risk.
	 */
	// WARNING! I haven't tested the distanceSE calculation. Convince yourself it is correct before using.
	public double getLogDetDistance(byte[] seqPat1, byte[] seqPat2) {
		
		// divergence matrix is called "D" in Massingham and Goldman (M&G) (citation above)
		for (int i=0; i<numStates; i++) {
			for (int j=0; j<numStates; j++) {
				divergence.m[i][j] = 0;
			}
		}
		for (int site = 0; site < seqPat1.length; site++) {
			byte col = seqPat1[site];
			byte row = seqPat2[site];
			if (col >=0 && col < numStates && row >=0 && row < numStates) {
				// site is 'real' on both sequences (no gaps or unknowns)
				divergence.m[row][col] += weight[site];
			}
		}
		double[] distWithError = LogDetDistance.distanceWithError(divergence);
		distance = distWithError[0];
		distanceSE = distWithError[1];
		return distance;
	}		

	
	/**
	 * Return divergence matrix of last pairwise distance calculated
	 * WARNING!!! Implementation incomplete!!!
	 * Only gets recalculated for LogDet distance calculations!
	 * @return Divergence matrix
	 */
	public int[][] getDivergence() {
		return divergence.getMatrixAsShortInts();
	}
	public IntMatrix getDivergenceAsIntMatrix() {
		return divergence;
	}
}
