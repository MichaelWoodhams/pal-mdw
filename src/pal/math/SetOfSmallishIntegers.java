package pal.math;

import java.util.Arrays;

import pal.misc.Utils;

/**
 * Represents a set of integers between 0 and maxvalue,
 * where maxvalue is moderately small. (Uses 'maxvalue' bits
 * of memory, or a little more.)
 * 
 * Is based of 'SetOfSmallIntegers', except it can 
 * hold sets of more than 64 elements, and it knows what its
 * max value is (i.e. the universal set is 0..maxvalue).
 * 
 * TODO: I'd like this and SetOfSmallIntegers to share
 * either a superclass or an interface, and to be able to have
 * the appropriate one chosen transparently at creation time,
 * but I haven't figured out a tidy way to do it yet.
 * 
 * If you know set will be smaller than size 64, and don't
 * need to know maxvalue, then SetOfSmallIntegers will be more efficient.
 * 
 * @author Michael
 *
 */
public class SetOfSmallishIntegers {
	private long[] bits;
	private int maxvalue;
	
	/**
	 * Creates an empty set
	 */
	public SetOfSmallishIntegers(int maxvalue) {
		int nWords = maxvalue/64+1; // NB: an integer division
		bits = new long[nWords];
		this.maxvalue=maxvalue;
	}

	/**
	 * Creates a set from the list of small integers
	 * @param list
	 */
	public SetOfSmallishIntegers(int maxvalue, int[] list) {
		this(maxvalue,list,false);
	}
	/**
	 * Creates a set from the list of small integers
	 * @param list
	 * @param enforceUnique If true, will throw an IllegalArgumentException if 'list' contains duplicates.
	 */
	public SetOfSmallishIntegers(int maxvalue, int[] list, boolean enforceUnique) {
		this(maxvalue);
		set(list,enforceUnique);
	}
	/**
	 * Creates a set from the binary (presence/absence) encoding in 'bits'.
	 * @param bits
	 */
	public SetOfSmallishIntegers(int maxvalue, long bits[]) {
		// TODO: check we don't set any bits beyond maxvalue
		this(maxvalue);
		// if bits longer than this.bits, will throw an error (as it should)
		// but better if I make that error more informative.
		for (int i=0; i<bits.length; i++) {
			this.bits[i] = bits[i];
		}
	}

	/**
	 * Creates a set from a boolean array
	 * @param isMember
	 */
	public SetOfSmallishIntegers(boolean[] isMember) {
		this(isMember.length-1);
		// Simple, but not as efficient as it could be
		for (int i=0; i<isMember.length; i++) {
			if (isMember[i]) {this.addToSet(i);}
		}
	}

	/**
	 *  copy constructor
	 * @param toCopy
	 */
	public SetOfSmallishIntegers(SetOfSmallishIntegers toCopy) {
		this.maxvalue = toCopy.maxvalue;
		this.bits = Arrays.copyOf(bits, bits.length);
	}
	
	/**
	 * 'set' in the sense of get/set methods.
	 * Make this set the same as the parameter.
	 * @param toCopy
	 */
	public void set(SetOfSmallishIntegers toCopy) {
		this.maxvalue = toCopy.maxvalue;
		this.bits = Arrays.copyOf(bits, bits.length);
	}
	
	/**
	 * 'set' in the sense of get/set methods.
	 * Make this set contain the elements in 'list'. (Elements of 'list' must be *small*: 0 to 63.)
	 * @param list
	 */
	public void set(int[] list) {
		set(list, false);
	}
	
	/**
	 * 'set' in the sense of get/set methods.
	 * Make this set contain the elements in 'list'. (Elements of 'list' must be *small*: 0 to 63.)
	 * @param list
	 * @param enforceUnique If true, will throw an error if 'list' contains duplicates.
	 */
	public void set(int[] list, boolean enforceUnique) {
		int n = list.length;
		this.setToEmpty();
		for (int i=0; i<n; i++) {
			if (list[i]<0 || list[i]>maxvalue) {
				throw new IllegalArgumentException("Set member number too large (maxvalue="+Integer.toString(maxvalue)+")");
			}
			int wordOffset = list[i]/64;
			int bitOffset = list[i]%64;
			long x = ((long)1)<<bitOffset;
			if (enforceUnique && (bits[wordOffset] & x) !=0) {
				throw new IllegalArgumentException("Set members must be unique");
			}
			bits[wordOffset] |= x;
		}
	}
	
	
	
	public int[] getList() {
		int count = 0;
		for (long copy : bits) {
			while (copy != 0) {
				if (copy%2 == 1) {
					count++;
				}
				copy >>>= 1;
			}
		}
		int[] list = new int[count];
		count = 0;
		int base = 0; // goes up by 64 for each word
		for (long copy : bits) {
			int i=0;
			while (copy != 0) {
				if (copy%2 == 1) {
					list[count++] = i+base;
				}
				copy >>>= 1;
				i++;
			}
			base += 64;
		}
		return list;
	}
	
	
	public void setToEmpty() {
		for (int i=0; i<bits.length; i++) bits[i] = 0;
	}
	
	public boolean isEmpty() {
		boolean empty = true;
		for (int i=0; i<bits.length; i++) empty = empty && (bits[i] == 0);
		return empty;
	}
	
	// Complicated form of 'equals' required to be able to use this as a hash key
	public boolean equals(Object that) {
		if (!(that instanceof SetOfSmallishIntegers)) return false;
		return equals((SetOfSmallishIntegers)that);
	}
	// Simpler version for when the compiler knows type of 'that'.
	public boolean equals(SetOfSmallishIntegers that) {
		if (this.maxvalue!=that.maxvalue) return false;
		for (int i=0; i<this.bits.length; i++) {
			if (this.bits[i]!=that.bits[i]) return false;
		}
		return true;
	}
	/* 
	 * @see java.lang.Object#hashCode()
	 * 
	 * I'm not sure that default wouldn't work - I'm worried that
	 * two 'bits' arrays, distinct but with same values, would
	 * cause different hashcode for default implementation.
	 */
	private static final int HASH_MULTIPLIER = 1234567890;
	public int hashCode() {
		int code = maxvalue*HASH_MULTIPLIER;
		for (int i=0; i<bits.length; i++) {
			// same result as (new Long(bits[i])).hashCode();
			int longHashCode = (int)(bits[i]^(bits[i]>>>32));
			code += 3*(i+1)*longHashCode; 
		}
		// Can't use 'code += bits.hashCode();' as it depends on the mem location of bits, rather than its contents.
		return code;
	}
	public SetOfSmallishIntegers union(SetOfSmallishIntegers that) {
		SetOfSmallishIntegers unionSet = new SetOfSmallishIntegers(this);
		unionSet.unionEquals(that);
		return unionSet;
	}
	public void unionEquals(SetOfSmallishIntegers that) {
		// TODO: Arguably this restriction is too harsh.
		if (this.maxvalue != that.maxvalue) 
			throw new IllegalArgumentException("Sets must have same sized superset to take union");
		for (int i=0; i<this.bits.length; i++) {
			this.bits[i] |= that.bits[i];
		}
	}
	public SetOfSmallishIntegers intersection(SetOfSmallishIntegers that) {
		SetOfSmallishIntegers intersectionSet = new SetOfSmallishIntegers(this);
		intersectionSet.intersectionEquals(that);
		return intersectionSet;
	}
	public void intersectionEquals(SetOfSmallishIntegers that) {
		// TODO: Arguably this restriction is too harsh.
		if (this.maxvalue != that.maxvalue) 
			throw new IllegalArgumentException("Sets must have same sized superset to take intersection");
		for (int i=0; i<this.bits.length; i++) {
			this.bits[i] &= that.bits[i];
		}
	}
	
	public boolean intersects(SetOfSmallishIntegers that) {
		// TODO: Arguably this restriction is too harsh.
		if (this.maxvalue != that.maxvalue) 
			throw new IllegalArgumentException("Sets must have same sized superset to take intersection");
		boolean intersects = false;
		for (int i=0; i<bits.length; i++) {
			intersects = intersects || ((this.bits[i] & that.bits[i])!=0);
		}
		return intersects;
	}
	public void addToSet(int newMember) {
		if (newMember<0 || newMember>maxvalue) {
			throw new IllegalArgumentException("Set members must all be between 0 and maxvalue="+Integer.toString(maxvalue));
		}
		int wordOffset = newMember/64;
		int bitOffset = newMember%64;
		long x = ((long)1)<<bitOffset;
		bits[wordOffset] |= x;
	}
	public boolean isMember(int value) {
		if (value<0 || value>maxvalue) {
			return false;
		}
		int wordOffset = value/64;
		int bitOffset = value%64;
		long x = ((long)1)<<bitOffset;
		
		return (bits[wordOffset]&x)!=0;
	}
	public String toString() {
		return "{"+Utils.toString(getList())+"}";
	}
	/*
	 * Complement the set 
	 */
	public void complement() {
		// all but last entry of 'bits' need to be XORed with 1111....111
		for (int i=0; i<bits.length-1; i++) {
			bits[i] = -(1+bits[i]);
		}
		int nRemainingBits = maxvalue%64;
		long bitmask = (((long)1)<<nRemainingBits)-1;
		bits[bits.length-1] ^= bitmask;
	}
	
	/**
	 * Unit test isn't nearly exhaustive.
	 */
	public static void unitTest() {
		boolean testFailed = false;
		SetOfSmallishIntegers a = new SetOfSmallishIntegers(10); // small set, one word
		SetOfSmallishIntegers b = new SetOfSmallishIntegers(127); // two words
		SetOfSmallishIntegers c = new SetOfSmallishIntegers(128); // three words

		testFailed = testFailed || a.bits.length != 1;
		testFailed = testFailed || b.bits.length != 2;
		testFailed = testFailed || c.bits.length != 3;

		int[] x = new int[]{12,0,3,6,12,128};
		int[] y = new int[]{1,32,33,63,64,128};
		try {
			a.set(x);
			testFailed = true; // should throw exception
		} catch (IllegalArgumentException e) {
		}
		try {
			b.set(x);
			testFailed = true; // should throw exception
		} catch (IllegalArgumentException e) {
		}
		try {
			c.set(x,true);
			testFailed = true; // should throw exception (duplicate value)
		} catch (IllegalArgumentException e) {
		}
		try {
			c.set(x);
		} catch (IllegalArgumentException e) {
			testFailed = true; // should not throw exception
		}
		
		SetOfSmallishIntegers d = new SetOfSmallishIntegers(128,y);
		testFailed = testFailed || !d.intersects(c);
		
		int[] yprimed = d.getList();
		testFailed = testFailed || (!Arrays.equals(y, yprimed));
		
		
		if (testFailed) {
			System.out.println("Unit test failed");
		} else {
			System.out.println("Unit test passed");
		}
	}
}
