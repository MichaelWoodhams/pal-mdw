package pal.math;
/**
 * A *very* simple class for some complex number calculations. Doesn't actually implement a complex
 * number object.
 * Created 2009/04/01
 * @author woodhams
 *
 */
public class Complex {
	
	public static double realPartComplexLog(double re, double im) {
		return 0.5*Math.log(re*re+im*im);
	}
	
	public static double imagPartComplexLog(double re, double im) {
		if (im != 0) {
			return 2*Math.atan(im/(Math.sqrt(re*re+im*im)+re));
		} else {
			if (re != 0) {
				return Math.PI;
			} else {
				return Double.NaN;
			}
		}
	}
	
	public static double realPartComplexExp(double re, double im) {
		return Math.exp(re)*Math.cos(im);
	}
	
	public static double imagPartComplexExp(double re, double im) {
		return Math.exp(re)*Math.sin(im);
	}
}
