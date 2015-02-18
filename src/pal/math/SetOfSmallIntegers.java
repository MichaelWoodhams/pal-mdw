package pal.math;

import pal.misc.Utils;

/**
 * Encode a set of integers (0-63) as a bitstring.
 * It is just a 'long'. It is only worth a class to make life easier later if we need to accommodate 
 * larger sets later on. 
 * 
 * TODO: two privates subclasses, one for up to size 64, one for larger sets. The correct subclass
 * will be transparently chosen and used.
 * 
 * TODO: Have an optional 'size of universal set' data member.
 * 
 * Functionality is not complete - mostly I'm adding only what I need.
 * @author woodhams
 *
 */

public class SetOfSmallIntegers {
	private long bits;
	public static final byte MAX_MEMBER=63; // Set contains numbers from 0 to MAX_MEMBER.
	
	/**
	 * Creates an empty set
	 */
	public SetOfSmallIntegers() {
		bits = 0;
	}
	/**
	 * Creates a set containing n members, 0..(n-1)
	 * @param n the number of members in the set
	 */
	public SetOfSmallIntegers(int n) {
		if (n>MAX_MEMBER+1) { throw new IllegalArgumentException("Set members must all be between 0 and MAX_MEMBER");}
		bits = 2^n -1;
	}
	/**
	 * Creates a set from the list of small integers
	 * @param list
	 */
	public SetOfSmallIntegers(int[] list) {
		set(list, false);
	}
	/**
	 * Creates a set from the list of small integers
	 * @param list
	 * @param enforceUnique If true, will throw an IllegalArgumentException if 'list' contains duplicates.
	 */
	public SetOfSmallIntegers(int[] list, boolean enforceUnique) {
		set(list,enforceUnique);
	}
	/**
	 * Creates a set from the binary (presence/absence) encoding in 'bits'.
	 * @param bits
	 */
	private SetOfSmallIntegers(long bits) {
		this.bits = bits;
	}
	
	/**
	 * 'set' in the sense of get/set methods.
	 * Make this set the same as the parameter.
	 * @param toCopy
	 */
	public void set(SetOfSmallIntegers toCopy) {
		this.bits = toCopy.bits;
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
		bits = 0;
		for (int i=0; i<n; i++) {
			if (list[i]<0 || list[i]>MAX_MEMBER) {
				throw new IllegalArgumentException("Set members must all be between 0 and MAX_MEMBER");
			}
			long x = ((long)1)<<list[i];
			if (enforceUnique && (bits & x) !=0) {
				throw new IllegalArgumentException("Set members must be unique");
			}
			bits = bits & x;
		}
	}
	
	
	/**
	 * Creates a set from a boolean array
	 * @param isMember
	 */
	public SetOfSmallIntegers(boolean[] isMember) {
		bits = 0;
		for (int i=0; i<isMember.length; i++) {
			if (isMember[i]) {this.addToSet(i);}
		}
	}
	
	public int[] getList() {
		int count = 0;
		long copy = bits;
		while (copy != 0) {
			if (copy%2 == 1) {
				count++;
			}
			copy >>>= 1;
		}
		int[] list = new int[count];
		count = 0;
		int i=0;
		copy = bits;
		while (copy > 0) {
			if (copy%2 == 1) {
				list[count++] = i;
			}
			copy >>>= 1;
			i++;
		}
		return list;
	}
	
	
	public void setToEmpty() {
		bits = 0;
	}
	
	public boolean isEmpty() {
		return (bits == 0);
	}
	
	// Complicated form of 'equals' required to be able to use this as a hash key
	public boolean equals(Object other) {
		return (other instanceof SetOfSmallIntegers && this.bits == ((SetOfSmallIntegers)other).bits);
	}
	// Simpler version for when the compiler knows type of 'other'.
	public boolean equals(SetOfSmallIntegers other) {
		return (this.bits == other.bits);
	}
	/* 
	 * @see java.lang.Object#hashCode()
	 * 
	 * Need to hash a long into an int, and furthermore, we want to avoid nearby
	 * 'long's mapping to nearby 'int's (for hashing efficiency.) 
	 * I do this by multiplying the long by a large number (on the order of Integer.MAX_VALUE)
	 * and truncating the result to an int.
	 * Note that the chosen number was not selected by some clever method, it is 
	 * off the top of my head.
	 */
	private static final long HASH_MULTIPLIER = 1234567890;
	public int hashCode() {
		//return (int)(bits*HASH_MULTIPLIER);
		int code = (int)(bits*HASH_MULTIPLIER);
		return code;
	}
	public SetOfSmallIntegers union(SetOfSmallIntegers other) {
		return new SetOfSmallIntegers(this.bits | other.bits);
	}
	public void unionEquals(SetOfSmallIntegers other) {
		this.bits = this.bits | other.bits;
	}
	public SetOfSmallIntegers intersection(SetOfSmallIntegers other) {
		return new SetOfSmallIntegers(this.bits & other.bits);
	}
	public void intersectionEquals(SetOfSmallIntegers other) {
		this.bits = this.bits & other.bits;
	}
	public boolean intersects(SetOfSmallIntegers other) {
		return ((this.bits & other.bits)!=0);
	}
	public void addToSet(int newMember) {
		if (newMember<0 || newMember>MAX_MEMBER) {
			throw new IllegalArgumentException("Set members must all be between 0 and MAX_MEMBER");
		}
		bits += ((long)1)<<newMember;
	}
	public boolean isMember(int value) {
		if (value<0 || value>MAX_MEMBER) {
			return false;
		}
		return ((((long)1)<<value & bits)!=0);
	}
	public String toString() {
		return "{"+Utils.toString(getList())+"}";
	}
	/*
	 * Complement of set, assuming there are n members of the universal set
	 * (i.e. 0 up to n-1.)
	 */
	public SetOfSmallIntegers complement(int n) {
		if (n<0 || n>MAX_MEMBER) {
			throw new IllegalArgumentException("Set members must all be between 0 and MAX_MEMBER");
		}
		long universalSet;
		// TODO: test whether special case is required
		if (n == MAX_MEMBER) {
			universalSet = -1;
		} else {
			universalSet = (((long)1)<<n)-1;
		}
		return new SetOfSmallIntegers(universalSet - bits);
	}
}
