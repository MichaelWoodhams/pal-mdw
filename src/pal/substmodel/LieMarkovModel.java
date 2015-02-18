package pal.substmodel;

import java.util.Arrays;

import pal.datatype.Nucleotides;

/**
 * Provides static methods to produce GeneralLinearRateMatrix's 
 * for various Lie Markov models. (See papers by Jeremy Sumner et al.)
 * 
 * This is a static class.
 * 
 * @author Michael
 *
 */

/*
 * TODO:
 * Reconsider this being a static class. It is subclassed by 
 * EfficientLieMarkovModel, which is not static.
 * Possibly this should be a subclass of GeneralLinearRateMatrix.
 */

public class LieMarkovModel {

	private static final int N = Nucleotides.DEFAULT_INSTANCE.getNumStates(); // 4
	
	private static final int A = Nucleotides.A_STATE;
	private static final int C = Nucleotides.C_STATE;
	private static final int G = Nucleotides.G_STATE;
	private static final int T = Nucleotides.UT_STATE;
	private static final int[] INDICES = new int[]{A,G,C,T}; // When printing matrices, will use this order for indices.
	
	// permute YR-conserving models to be Weak-Strong (Watson-Crick) pair conserving models
	private static final int[] PERMUTE_WS = new int[4]; 
	// permute YR-conserving models to have AC, GT be the conserved pairs
	private static final int[] PERMUTE_MK = new int[4];
	static {
		// Using a static block to maintain independence from which indices A, G, C, T correspond to
		PERMUTE_WS[A]=A;
		PERMUTE_WS[G]=T;
		PERMUTE_WS[C]=C;
		PERMUTE_WS[T]=G;
		PERMUTE_MK[A]=A;
		PERMUTE_MK[G]=C;
		PERMUTE_MK[C]=G;
		PERMUTE_MK[T]=T;
	}
	
	// 'AG' is the rate from A to G. 
	private static final int[] AG = new int[]{A, G};
	private static final int[] AC = new int[]{A, C};
	private static final int[] AT = new int[]{A, T};
	private static final int[] GA = new int[]{G, A};
	private static final int[] GC = new int[]{G, C};
	private static final int[] GT = new int[]{G, T};
	private static final int[] CA = new int[]{C, A};
	private static final int[] CG = new int[]{C, G};
	private static final int[] CT = new int[]{C, T};
	private static final int[] TA = new int[]{T, A};
	private static final int[] TG = new int[]{T, G};
	private static final int[] TC = new int[]{T, C};
	
	protected static final double[][][] orbit101  = permute(identity(), new int[][]{AC,AT,GC,GT,CA,CG,TA,TG});
	protected static final double[][][] orbit110  = permute(identity(), new int[][]{AG,GA,CT,TC});
	protected static final double[][][] orbit201a = permute(AG(),    new int[][]{AC,CA,GT,TG});
	protected static final double[][][] orbit201b = permute(AG(),    new int[][]{AT,TG,GC,CA});
	protected static final double[][][] orbit201c = permute(YR(),    new int[][]{AC,AT,GC,GT});
	protected static final double[][][] orbit210  = permute(YR(),    new int[][]{AG,GA});
	protected static final double[][][] orbit212  = permute(YR(),    new int[][]{AC,AG,AT,GA,GC,GT});
	protected static final double[][][] orbit401a = permute(even(),  new int[][]{AC,AT,CG,TG});
	protected static final double[][][] orbit401b = permute(AG_CT(), new int[][]{AC,AT,TA,TG});
	protected static final double[][][] orbit401c = permute(even(),  new int[][]{AC,AT});
	protected static final double[][][] orbit401d = permute(AG_CT(), new int[][]{AT,CG});
	protected static final double[][][] orbit401e = permute(AG_YR(), new int[][]{AT,GC});
	protected static final double[][][] orbit401f = permute(AG_CT(), new int[][]{AT,TA});
	protected static final double[][][] orbit410a = permute(even(),  new int[][]{AG});
	protected static final double[][][] orbit410b = permute(AG_CT(), new int[][]{AG,CT});
	protected static final double[][][] orbit411a = permute(AG_CT(), new int[][]{AG,AC,CA,CT});
	protected static final double[][][] orbit411b = permute(AG_CT(), new int[][]{AG,GC,CT,TA});
	protected static final double[][][] orbit412a = permute(even(),  new int[][]{AG,AC,AT});
	protected static final double[][][] orbit412b = permute(even(),  new int[][]{AG,CA,TA});
	protected static final double[][][] orbit412c = permute(even(),  new int[][]{AG,GC,GT});
	protected static final double[][][] orbit412d = permute(AG_CT(), new int[][]{AC,AG,AT,TA,TG,TC});
	protected static final double[][][] orbit412e = permute(AG_CT(), new int[][]{AG,GC,GT,TC,CA,CG});
	protected static final double[][][] orbit412f = permute(AG_CT(), new int[][]{AG,CA,TA,AT,GT,TC});
	protected static final double[][][] orbit414a = permute(even(),  new int[][]{AG,AG,AC,AC,AT,AT,CA,CG,TA,TG});
	protected static final double[][][] orbit414b = permute(even(),  new int[][]{AG,AG,GC,GC,GT,GT,CA,CG,TA,TG});
	protected static final double[][][] orbit414c = permute(even(),  new int[][]{AG,AG,CA,CA,TA,TA,AC,AT,GC,GT});
	protected static final double[][][] orbit416  = permute(even(),  new int[][]{GA,CG,TG,AC,AT,GC,GT});
	protected static final double[][][] orbit432  = permute(even(),  new int[][]{AG,CA,CT,TA,TC});
	protected static final double[][][] orbit801a = permute(S2wrS2(),new int[][]{AT}); // just used for inefficient implementation of GMM
	protected static final double[][][] orbit801b = permute(S2wrS2(),new int[][]{AC,AT,CG,TA});
	protected static final double[][][] orbit811  = permute(S2wrS2(),new int[][]{AG,TA,TA,TC});
	protected static final double[][][] orbit812a = permute(S2wrS2(),new int[][]{AT,TG,TC});
	protected static final double[][][] orbit812b = permute(S2wrS2(),new int[][]{AG,GC,CA});
	
	// WS6.6 is the strand symmetric model
	public static GeneralLinearRateMatrix getStrandSymmetricModel() {
		return getModel_ws6_6();
	}
	
	// model 2_2a not implemented as it forbids transversions entirely.
	public static GeneralLinearRateMatrix getModel2_2b() {
		double[][][] basis = concatAll(orbit110,orbit101);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel3_3a() {
		double[][][] basis = concatAll(orbit110,orbit201a);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel3_3b() {
		double[][][] basis = concatAll(orbit110,orbit201b);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel3_3c() {
		double[][][] basis = concatAll(orbit101,orbit210);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel3_4() {
		double[][][] basis = concatAll(orbit110, orbit101, orbit212);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel4_4a() {
		double[][][] basis = concatAll(orbit412a);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel4_4b() {
		double[][][] basis = concatAll(orbit210,orbit201c);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel4_5a() {
		double[][][] basis = concatAll(orbit110,orbit201a,orbit212);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel4_5b() {
		double[][][] basis = concatAll(orbit110,orbit201b,orbit212);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel5_6a() {
		double[][][] basis = concatAll(orbit201a,orbit201b,orbit210);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel5_6b() {
		double[][][] basis = concatAll(orbit110,orbit101,orbit412a);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel5_7a() {
		double[][][] basis = concatAll(orbit110,orbit201a,orbit412d);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel5_7b() {
		double[][][] basis = concatAll(orbit110,orbit201a,orbit412e);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel5_7c() {
		double[][][] basis = concatAll(orbit110,orbit201a,orbit412f);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel5_11a() {
		double[][][] basis = concatAll(orbit101,orbit210,orbit412d,orbit414a);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel5_11b() {
		double[][][] basis = concatAll(orbit101,orbit210,orbit412e,orbit414b);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel5_11c() {
		double[][][] basis = concatAll(orbit101,orbit210,orbit412f,orbit414c);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel5_16() {
		double[][][] basis = concatAll(orbit110,orbit101,orbit412f,orbit416,orbit212);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel6_6() {
		double[][][] basis = concatAll(orbit210,orbit401e);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel6_7a() {
		double[][][] basis = concatAll(orbit110,orbit201a,orbit412a);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel6_7b() {
		double[][][] basis = concatAll(orbit110,orbit201b,orbit412a);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel6_8a() {
		double[][][] basis = concatAll(orbit210,orbit201c,orbit412a);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel6_8b() {
		double[][][] basis = concatAll(orbit210,orbit201c,orbit412b);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel6_17a() {
		double[][][] basis = concatAll(orbit110,orbit201a,orbit412f,orbit416,orbit212,orbit432);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel6_17b() {
		double[][][] basis = concatAll(orbit110,orbit201b,orbit412f,orbit416,orbit212,orbit432);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel8_8() {
		double[][][] basis = concatAll(orbit410a,orbit401c);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel8_10a() {
		double[][][] basis = concatAll(orbit210,orbit401e,orbit412a);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel8_10b() {
		double[][][] basis = concatAll(orbit210,orbit401e,orbit412b);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel8_16() {
		double[][][] basis = concatAll(orbit210,orbit201c,orbit412a,orbit412b,orbit401a);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel8_17() {
		double[][][] basis = concatAll(orbit110,orbit412a,orbit432,orbit411a,orbit401d);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel8_18() {
		double[][][] basis = concatAll(orbit201a,orbit412a,orbit401b,orbit410b,orbit412c);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel9_20a() {
		double[][][] basis = concatAll(orbit201a,orbit201b,orbit410a,orbit801b,orbit401b);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel9_20b() {
		double[][][] basis = concatAll(orbit201b,orbit210,orbit812b,orbit411b,orbit401f);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel10_12() {
		double[][][] basis = concatAll(orbit410a,orbit401c,orbit401e);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	public static GeneralLinearRateMatrix getModel10_34() {
		double[][][] basis = concatAll(orbit210,orbit412a,orbit412b,orbit401d,orbit401e,orbit812a,orbit811);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	// 12.12 is the GMM (General Markov Model) implemented in an inefficient way. 
	public static GeneralLinearRateMatrix getModel12_12() {
		double[][][] basis = concatAll(orbit410a,orbit801a);
		return new GeneralLinearRateMatrix(4, null, basis, true); 
	}
	// Due to symmetry, 3.3a = WS3.3a = MK3.3a; 4.4a = WS4.4a = MK4.4a; 6.7a = WS6.7a = MK6.7a
	public static GeneralLinearRateMatrix getModel_ws2_2b()  {return new GeneralLinearRateMatrix(getModel2_2b(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws3_3a()  {return new GeneralLinearRateMatrix(getModel3_3a(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws3_3b()  {return new GeneralLinearRateMatrix(getModel3_3b(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws3_3c()  {return new GeneralLinearRateMatrix(getModel3_3c(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws3_4()   {return new GeneralLinearRateMatrix(getModel3_4(),  PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws4_4a()  {return new GeneralLinearRateMatrix(getModel4_4a(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws4_4b()  {return new GeneralLinearRateMatrix(getModel4_4b(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws4_5a()  {return new GeneralLinearRateMatrix(getModel4_5a(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws4_5b()  {return new GeneralLinearRateMatrix(getModel4_5b(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws5_6a()  {return new GeneralLinearRateMatrix(getModel5_6a(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws5_6b()  {return new GeneralLinearRateMatrix(getModel5_6b(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws5_7a()  {return new GeneralLinearRateMatrix(getModel5_7a(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws5_7b()  {return new GeneralLinearRateMatrix(getModel5_7b(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws5_7c()  {return new GeneralLinearRateMatrix(getModel5_7c(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws5_11a() {return new GeneralLinearRateMatrix(getModel5_11a(),PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws5_11b() {return new GeneralLinearRateMatrix(getModel5_11b(),PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws5_11c() {return new GeneralLinearRateMatrix(getModel5_11c(),PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws5_16()  {return new GeneralLinearRateMatrix(getModel5_16(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws6_6()   {return new GeneralLinearRateMatrix(getModel6_6(),  PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws6_7a()  {return new GeneralLinearRateMatrix(getModel6_7a(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws6_7b()  {return new GeneralLinearRateMatrix(getModel6_7b(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws6_8a()  {return new GeneralLinearRateMatrix(getModel6_8a(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws6_8b()  {return new GeneralLinearRateMatrix(getModel6_8b(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws6_17a() {return new GeneralLinearRateMatrix(getModel6_17a(),PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws6_17b() {return new GeneralLinearRateMatrix(getModel6_17b(),PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws8_8()   {return new GeneralLinearRateMatrix(getModel8_8(),  PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws8_10a() {return new GeneralLinearRateMatrix(getModel8_10a(),PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws8_10b() {return new GeneralLinearRateMatrix(getModel8_10b(),PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws8_16()  {return new GeneralLinearRateMatrix(getModel8_16(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws8_17()  {return new GeneralLinearRateMatrix(getModel8_17(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws8_18()  {return new GeneralLinearRateMatrix(getModel8_18(), PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws9_20a() {return new GeneralLinearRateMatrix(getModel9_20a(),PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws9_20b() {return new GeneralLinearRateMatrix(getModel9_20b(),PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws10_12() {return new GeneralLinearRateMatrix(getModel10_12(),PERMUTE_WS);}
	public static GeneralLinearRateMatrix getModel_ws10_34() {return new GeneralLinearRateMatrix(getModel10_34(),PERMUTE_WS);}

	public static GeneralLinearRateMatrix getModel_mk2_2b()  {return new GeneralLinearRateMatrix(getModel2_2b(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk3_3a()  {return new GeneralLinearRateMatrix(getModel3_3a(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk3_3b()  {return new GeneralLinearRateMatrix(getModel3_3b(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk3_3c()  {return new GeneralLinearRateMatrix(getModel3_3c(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk3_4()   {return new GeneralLinearRateMatrix(getModel3_4(),  PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk4_4a()  {return new GeneralLinearRateMatrix(getModel4_4a(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk4_4b()  {return new GeneralLinearRateMatrix(getModel4_4b(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk4_5a()  {return new GeneralLinearRateMatrix(getModel4_5a(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk4_5b()  {return new GeneralLinearRateMatrix(getModel4_5b(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk5_6a()  {return new GeneralLinearRateMatrix(getModel5_6a(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk5_6b()  {return new GeneralLinearRateMatrix(getModel5_6b(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk5_7a()  {return new GeneralLinearRateMatrix(getModel5_7a(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk5_7b()  {return new GeneralLinearRateMatrix(getModel5_7b(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk5_7c()  {return new GeneralLinearRateMatrix(getModel5_7c(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk5_11a() {return new GeneralLinearRateMatrix(getModel5_11a(),PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk5_11b() {return new GeneralLinearRateMatrix(getModel5_11b(),PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk5_11c() {return new GeneralLinearRateMatrix(getModel5_11c(),PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk5_16()  {return new GeneralLinearRateMatrix(getModel5_16(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk6_6()   {return new GeneralLinearRateMatrix(getModel6_6(),  PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk6_7a()  {return new GeneralLinearRateMatrix(getModel6_7a(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk6_7b()  {return new GeneralLinearRateMatrix(getModel6_7b(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk6_8a()  {return new GeneralLinearRateMatrix(getModel6_8a(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk6_8b()  {return new GeneralLinearRateMatrix(getModel6_8b(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk6_17a() {return new GeneralLinearRateMatrix(getModel6_17a(),PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk6_17b() {return new GeneralLinearRateMatrix(getModel6_17b(),PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk8_8()   {return new GeneralLinearRateMatrix(getModel8_8(),  PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk8_10a() {return new GeneralLinearRateMatrix(getModel8_10a(),PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk8_10b() {return new GeneralLinearRateMatrix(getModel8_10b(),PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk8_16()  {return new GeneralLinearRateMatrix(getModel8_16(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk8_17()  {return new GeneralLinearRateMatrix(getModel8_17(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk8_18()  {return new GeneralLinearRateMatrix(getModel8_18(), PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk9_20a() {return new GeneralLinearRateMatrix(getModel9_20a(),PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk9_20b() {return new GeneralLinearRateMatrix(getModel9_20b(),PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk10_12() {return new GeneralLinearRateMatrix(getModel10_12(),PERMUTE_MK);}
	public static GeneralLinearRateMatrix getModel_mk10_34() {return new GeneralLinearRateMatrix(getModel10_34(),PERMUTE_MK);}




	// The permutation groups:
	// The method calculates and caches the permutation first time it is used,
	// and just returns in on subsequent calls.
	// S2wrS2 is the full symmetry group of permutations (AG), (CT), (AC)(GT).
	private static int[][] S2wrS2perm = null;
	private static int[][] S2wrS2() {
		if (S2wrS2perm == null) {
			S2wrS2perm = new int[8][4];
			// Identity permutation:
			S2wrS2perm[0][A] = A;
			S2wrS2perm[0][G] = G;
			S2wrS2perm[0][C] = C;
			S2wrS2perm[0][T] = T;
			// (AG) permutation:
			S2wrS2perm[1][A] = G;
			S2wrS2perm[1][G] = A;
			S2wrS2perm[1][C] = C;
			S2wrS2perm[1][T] = T;
			// (CT) permutation:
			S2wrS2perm[2][A] = A;
			S2wrS2perm[2][G] = G;
			S2wrS2perm[2][C] = T;
			S2wrS2perm[2][T] = C;
			// (AG)(CT) 
			S2wrS2perm[3][A] = G;
			S2wrS2perm[3][G] = A;
			S2wrS2perm[3][C] = T;
			S2wrS2perm[3][T] = C;
			// (AC)(GT)
			S2wrS2perm[4][A] = C;
			S2wrS2perm[4][G] = T;
			S2wrS2perm[4][C] = A;
			S2wrS2perm[4][T] = G;
			// (AC)(GT) . (AG):
			S2wrS2perm[5][A] = C;
			S2wrS2perm[5][G] = T;
			S2wrS2perm[5][C] = G;
			S2wrS2perm[5][T] = A;
			// (AC)(GT) . (CT):
			S2wrS2perm[6][A] = T;
			S2wrS2perm[6][G] = C;
			S2wrS2perm[6][C] = A;
			S2wrS2perm[6][T] = G;
			// (AC)(GT) . (AG)(CT):
			S2wrS2perm[7][A] = T;
			S2wrS2perm[7][G] = C;
			S2wrS2perm[7][C] = G;
			S2wrS2perm[7][T] = A;
		}
		return S2wrS2perm;
	}
	// Subgroup of S2wrS2: just the even permutations, i.e. generated by (AG)(CT) and (AC)(GT)
	private static int[][] evenPerm = null;
	private static int[][] even() {
		if (evenPerm == null) {
			evenPerm = new int[4][4];
			// Identity permutation:
			evenPerm[0][A] = A;
			evenPerm[0][G] = G;
			evenPerm[0][C] = C;
			evenPerm[0][T] = T;
			// (AG)(CT) permutation:
			evenPerm[1][A] = G;
			evenPerm[1][G] = A;
			evenPerm[1][C] = T;
			evenPerm[1][T] = C;
			// YR permutation, i.e. (AC)(GT)
			evenPerm[2][A] = C;
			evenPerm[2][G] = T;
			evenPerm[2][C] = A;
			evenPerm[2][T] = G;
			// Both permutations
			evenPerm[3][A] = T;
			evenPerm[3][G] = C;
			evenPerm[3][C] = G;
			evenPerm[3][T] = A;
		}
		return evenPerm;
	}
	// Subgroup of S2wrS2: generated by (AC)(GT) only (i.e. exchange purines and pyrmidines)
	private static int[][] YRperm = null;
	private static int[][] YR() {
		if (YRperm == null) {
			YRperm = new int[2][4];
			// Identity permutation:
			YRperm[0][A] = A;
			YRperm[0][G] = G;
			YRperm[0][C] = C;
			YRperm[0][T] = T;
			// YR permutation, i.e. (AC)(GT)
			YRperm[1][A] = C;
			YRperm[1][G] = T;
			YRperm[1][C] = A;
			YRperm[1][T] = G;
		}
		return YRperm;
	}
	// Subgroup maintaining purine/pyrimidine distinction, i.e. generated by (AG), (CT) 
	private static int[][] AG_CTperm = null;
	private static int[][] AG_CT() {
		if (AG_CTperm == null) {
			AG_CTperm = new int[4][4];
			// Identity permutation:
			AG_CTperm[0][A] = A;
			AG_CTperm[0][G] = G;
			AG_CTperm[0][C] = C;
			AG_CTperm[0][T] = T;
			// (AG) permutation:
			AG_CTperm[1][A] = G;
			AG_CTperm[1][G] = A;
			AG_CTperm[1][C] = C;
			AG_CTperm[1][T] = T;
			// (CT) permutation:
			AG_CTperm[2][A] = A;
			AG_CTperm[2][G] = G;
			AG_CTperm[2][C] = T;
			AG_CTperm[2][T] = C;
			// Both permutations
			AG_CTperm[3][A] = G;
			AG_CTperm[3][G] = A;
			AG_CTperm[3][C] = T;
			AG_CTperm[3][T] = C;
		}
		return AG_CTperm;
	}
	
	// Subgroup generated by (AG), (AC)(GT)) 
	private static int[][] AG_YRperm = null;
	private static int[][] AG_YR() {
		if (AG_YRperm == null) {
			AG_YRperm = new int[4][4];
			// Identity permutation:
			AG_YRperm[0][A] = A;
			AG_YRperm[0][G] = G;
			AG_YRperm[0][C] = C;
			AG_YRperm[0][T] = T;
			// (AG) permutation:
			AG_YRperm[1][A] = G;
			AG_YRperm[1][G] = A;
			AG_YRperm[1][C] = C;
			AG_YRperm[1][T] = T;
			// YR permutation, i.e. (AC)(GT)
			AG_YRperm[2][A] = C;
			AG_YRperm[2][G] = T;
			AG_YRperm[2][C] = A;
			AG_YRperm[2][T] = G;
			// Both permutations
			AG_YRperm[3][A] = T;
			AG_YRperm[3][G] = C;
			AG_YRperm[3][C] = A;
			AG_YRperm[3][T] = G;
		}
		return AG_YRperm;
	}
	
	// Subgroup swapping AG
	private static int[][] AGperm = null;
	private static int[][] AG() {
		if (AGperm == null) {
			AGperm = new int[2][4];
			// Identity permutation:
			AGperm[0][A] = A;
			AGperm[0][G] = G;
			AGperm[0][C] = C;
			AGperm[0][T] = T;
			// (AG) permutation:
			AGperm[1][A] = G;
			AGperm[1][G] = A;
			AGperm[1][C] = C;
			AGperm[1][T] = T;
		}
		return AGperm;
	}
	// Identity
	private static int[][] identityPerm = null;
	private static int[][] identity() {
		if (identityPerm == null) {
			identityPerm = new int[1][4];
			// Identity permutation:
			identityPerm[0][A] = A;
			identityPerm[0][G] = G;
			identityPerm[0][C] = C;
			identityPerm[0][T] = T;
		}
		return identityPerm;
	}
	private static double[][][] permute(int[][] perm, int[][] basePairs) {
		double[][][] orbit = new double[perm.length][N][N];
		double[][] exemplar = exemplar(basePairs);
		for (int mat=0; mat<perm.length; mat++) {
			for (int row : INDICES) {
				for (int col : INDICES) {
					orbit[mat][row][col] = exemplar[perm[mat][row]][perm[mat][col]];
				}
			}
		}
		return orbit;
	}
	
	private static double[][] exemplar(int[][] basePairs) {
		double[][] q = new double[N][N];
		double unit = 1./basePairs.length;
		for (int[] pair : basePairs) {
			q[pair[1]][pair[0]] += unit; // ensures normalization: sum off diagonal elements = 1.
		}
		// And implement the rows-sum-to-one condition.
		for (int row=0; row<N; row++) {
			double rowsum = 0;
			for (int col=0; col<N; col++) {
				rowsum += q[row][col];
			}
			q[row][row] = -rowsum;
		}
		return q;
	}
	
	// Adapted from code at http://stackoverflow.com/questions/80476/how-to-concatenate-two-arrays-in-java
	private static double[][][] concatAll(double[][][] first, double[][][]... rest) {
		int totalLength = first.length;
		for (double[][][] array : rest) {
			totalLength += array.length;
		}
		double[][][] result = Arrays.copyOf(first, totalLength);
		int offset = first.length;
		for (double[][][] array : rest) {
			System.arraycopy(array, 0, result, offset, array.length);
			offset += array.length;
		}
		return result;
	}
	
	
	private static void printOrbit(double[][][] f) {
		double minEntry = Double.MAX_VALUE;
		for (int row=0; row<N; row++) {
			for (int col=0; col<N; col++) {
				if (f[0][row][col] > 0) minEntry = pal.math.MathUtils.min(minEntry, f[0][row][col]);
			}
		}
		double scale = 1/minEntry;
		for (int col : INDICES) {
			for (int mat=0; mat<f.length; mat++) {
				System.out.printf("[%2.0f, %2.0f, %2.0f, %2.0f]  ", 
						f[mat][A][col]*scale, f[mat][G][col]*scale, f[mat][C][col]*scale, f[mat][T][col]*scale);
			}
			System.out.println();
		}
		System.out.println();
	}
	
	/**
	 * Prints out all the ray families, with indices in AGCT order, and transposed so that in display, 
	 * columns sum to zero.
	 * 
	 * Primarily for debugging purposes.
	 */
	public static void printOrbits() {
		System.out.println("Orbit 101: "); printOrbit(orbit101);
		System.out.println("Orbit 110: "); printOrbit(orbit110);
		System.out.println("Orbit 201a:"); printOrbit(orbit201a);
		System.out.println("Orbit 201b:"); printOrbit(orbit201b);
		System.out.println("Orbit 201c:"); printOrbit(orbit201c);
		System.out.println("Orbit 210: "); printOrbit(orbit210);
		System.out.println("Orbit 212: "); printOrbit(orbit212);
		System.out.println("Orbit 401c:"); printOrbit(orbit401c);
		System.out.println("Orbit 401d:"); printOrbit(orbit401d);
		System.out.println("Orbit 401e:"); printOrbit(orbit401e);
		System.out.println("Orbit 401f:"); printOrbit(orbit401f);
		System.out.println("Orbit 401a:"); printOrbit(orbit401a);
		System.out.println("Orbit 401b:"); printOrbit(orbit401b);
		System.out.println("Orbit 410a:"); printOrbit(orbit410a);
		System.out.println("Orbit 412a:"); printOrbit(orbit412a);
		System.out.println("Orbit 412b:"); printOrbit(orbit412b);
		System.out.println("Orbit 412c:"); printOrbit(orbit412c);
		System.out.println("Orbit 416: "); printOrbit(orbit416);
		System.out.println("Orbit 410b:"); printOrbit(orbit410b);
		System.out.println("Orbit 411a:"); printOrbit(orbit411a);
		System.out.println("Orbit 411b:"); printOrbit(orbit411b);		
		System.out.println("Orbit 412d:"); printOrbit(orbit412d);
		System.out.println("Orbit 412e:"); printOrbit(orbit412e);
		System.out.println("Orbit 412f:"); printOrbit(orbit412f);
		System.out.println("Orbit 414a:"); printOrbit(orbit414a);
		System.out.println("Orbit 414b:"); printOrbit(orbit414b);
		System.out.println("Orbit 414c:"); printOrbit(orbit414c);
		System.out.println("Orbit 432: "); printOrbit(orbit432);
		System.out.println("Orbit 801b:"); printOrbit(orbit801b);		
		System.out.println("Orbit 812a:"); printOrbit(orbit812a);		
		System.out.println("Orbit 812b:"); printOrbit(orbit812b);
		System.out.println("Orbit 811: "); printOrbit(orbit811);
	}
}
