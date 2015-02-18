package pal.algorithmics;

/**
 * <p>Title: UndoableAction</p>
 * <p>Description: A stateful, single thread object</p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 * @author Matthew Goode
 * @version 1.0
 */
import java.io.PrintWriter;
import java.util.*;

import pal.math.MathUtils;
import pal.math.MersenneTwisterFast;
public interface UndoableAction {
  /**
   * Perform an action
   * @param currentscore The current score before doing the action
   * @param desparationValue An indication by the processing machines of willingness to do more extreme actions. A value of 0 means not desparate at all, a value of 1 means very desparate
   * @return the current score after doing the action (or the input score if not successful)
   */
	
	// MDW: I've added crude statistics gathering. Could be upgraded: return a data structure instead of just printing
	// statistics, have a static flag to turn on or off all stats gathering within this class etc.
	// I've also added 'getDescription()' to aid in debugging output.
  public double doAction(double currentScore, double desparationValue);

  /**
   * Was the last action deterministic? That is, if it wasn't chosen and state is still as
   * before is it worth doing it again?
   * @return true if last action deterministic
   */
  public boolean isActionDeterministic();

  /**
   * Was the last call to doAction() succesful?
   * @return true if last action successful, false otherwise
   */
  public boolean isActionSuccessful();
  /**
   * Undo the last action (if it was successful)
   * Users of undoable actions should accept that sometimes it isn't possible.
   * If an undo was not possible the action should not change any state
   * @return true if undo was successful
   */
  public boolean undoAction();
  
  /**
   * Print statistics on how this action has performed.
   * Details are up to the class - usually this method does nothing.
   * @param out
   */
  public void printStats(PrintWriter out);
  
  /**
   * 
   * @return A string describing the action. If some sort of multi-action, describe which sub-action was performed last.
   */
  public String getDescription();

// -=-=-==-=--=-==--=-==--=-==---==-=--=-=-=-=-==-=-=-=-=-=-=--==-=-=--==-=--=-
// -=-=-= Utils -=-=-=-=-=-=-==--==-=--=-=-=-=-==-=-=-=-=-=-=--==-=-=--==-=--=-
// -=-=-==-=--=-==--=-==--=-==---==-=--=-=-=-=-==-=-=-=-=-=-=--==-=-=--==-=--=-

  public static final class Utils {
    /**
     * Create an action that selects uniformly from a number of sub actions
     * @param subActions
     * @return
     */
    public static final UndoableAction getSimpleUniformSelection(UndoableAction[] subActions) {
      return new Multi(subActions);
    }
    /**
     * Create an action that selects uniformly from a number of sub actions
     * @param subActions
     * @param actionProportions
     * @throws IllegalArgumentException if action array and proportion arrays are different lengths
     * @return
     */
    public static final UndoableAction getDistributedSelection(UndoableAction[] subActions, double[] actionProportions) {
    	return getDistributedSelection(subActions, actionProportions, 0);
    }
    /**
     * Create an action that selects uniformly from a number of sub actions.
     * Chances are you're better (or at least no worse) to use getShuffledSelection instead.
     * @param subActions
     * @param actionProportions
     * @param seed random number generator seed - if zero, will seed from system clock
     * @throws IllegalArgumentException if action array and proportion arrays are different lengths
     * @return
     */
    public static final UndoableAction getDistributedSelection(UndoableAction[] subActions, double[] actionProportions, long seed) {
      if(subActions.length>actionProportions.length) {
        throw new IllegalArgumentException("Actions and proportion array different lengths");
      }
      return new DistributedMulti(subActions, actionProportions, seed);
    }
    
    /**
     * Very similar to getDistributedSelection, except that the actions are shuffled rather than 
     * being randomly and independently chosen each time. This gives a more even coverage of the 
     * various sub-actions. Don't make the numbers in actionProportions very large - 'a few' is
     * probably optimal, more than a few hundred is wasteful and pointless (use getDistributedSelection instead.)
     * @param subActions
     * @param actionProportions
     * @param seed
     * @return
     */
    public static final UndoableAction getShuffledSelection(UndoableAction[] subActions, int[] actionProportions, long seed) {
        if(subActions.length>actionProportions.length) {
          throw new IllegalArgumentException("Actions and proportion array different lengths");
        }
        return new ProfiledMulti(subActions, actionProportions, seed);
      }

    /**
     * Just like getDistributedAction, except it maintains statistics on how each action is performing.
     * (These can be accessed by the printStats method.)
     * 
     */
    public static final UndoableAction getProfiledShuffledSelection(UndoableAction[] subActions, int[] actionProportions, long seed) {
      if(subActions.length>actionProportions.length) {
        throw new IllegalArgumentException("Actions and proportion array different lengths");
      }
      return new ProfiledMulti(subActions, actionProportions, seed);
    }

		/**
     * Create an action that combines multiple actions
     * @param subActions The actions that are do in turn.
     * @return An action that performs all the sub actions
     */
    public static final UndoableAction getCombined(UndoableAction[] subActions) {
      return new Combined(subActions);
    }
    /**
     * A simple tool for change actions when things get desparate
     * @param primaryAction The main action to do when things are going well
     * @param desparateAction The action to do when things get desparate. The desperation value for the desparate action will be scaled according to how much over the limit we are
     * @param desparationLimit The desparate value at which we start doing the desparate action
     * @param desparationInterval The time between desparate actions when we cross the cutoff (a value of one will mean do all the time after desparation value has crossed cutoff)
     */
    public static final UndoableAction getSimpleDesparation(UndoableAction primaryAction, UndoableAction desparateAction, double desparationLimit, int desparationInterval) {
      return new SimpleDesparation(primaryAction,desparateAction,desparationLimit,desparationInterval);
    }

    // -=-==-=--==--=-=-=-=-=-=-==--=-=
    private static class SimpleDesparation implements UndoableAction {
      private final UndoableAction primaryAction_;
      private final UndoableAction desparateAction_;
      private final double desparationLimit_;
      private final int desparationInterval_;
      private int currentDesparateCount_ = 0;
      private UndoableAction lastAction_ = null;
      /**
       * A simple tool for change actions when things get desparate
       * @param primaryAction The main action to do when things are going well
       * @param desparateAction The action to do when things get desparate. The desperation value for the desparate action will be scaled according to how much over the limit we are
       * @param desparationLimit The desparate value at which we start doing the desparate action
       * @param desparationInterval The time between desparate actions when we cross the cutoff (a value of one will mean do all the time after desparation value has crossed cutoff)
       */
      public SimpleDesparation(UndoableAction primaryAction, UndoableAction desparateAction, double desparationLimit, int desparationInterval) {
        this.primaryAction_ = primaryAction;
        this.desparateAction_ = desparateAction;
        this.desparationLimit_ = desparationLimit;
        this.desparationInterval_ = desparationInterval;
      }
      /**
       * @return false
       */
      public boolean isActionDeterministic() {
        return false;
      }
      public double doAction(double currentScore, double desparationValue) {
        if(desparationValue>=desparationLimit_) {
          currentDesparateCount_++;
          if(currentDesparateCount_==desparationInterval_) {
            currentDesparateCount_ = 0;
            lastAction_ = desparateAction_;
            desparationValue = (desparationLimit_-desparationValue)/(1-desparationLimit_);
          } else {
            lastAction_ = primaryAction_;
          }
        } else {
          lastAction_ = primaryAction_;
          currentDesparateCount_ = 0;
        }
        return lastAction_.doAction(currentScore,desparationValue);
      }
      public boolean isActionSuccessful() {
        if(lastAction_!=null) {
          return lastAction_.isActionSuccessful();
        }
        throw new RuntimeException("Assertion error : isActionSuccessful() called when no action has been done recently");
      }
      public boolean undoAction() {
        if(lastAction_!=null) {
          final boolean successful = lastAction_.undoAction();
          lastAction_ = null;
          return successful;
        } else {
          throw new RuntimeException("Assertion error : undoAction() called when no action has been done recently (or has already been undone)");
        }
      }
      public void printStats(PrintWriter out) {};
      public String getDescription() {
    	  return "SimpleDesparation(" + ((lastAction_ == primaryAction_) ? "primary=" : "desparation=") + lastAction_.getDescription() + ")"; 
      }
    } //End of class SimpleDesparation

    // -=-==-=--==--=-=-=-=-=-=-==--=-=
    private static class Multi implements UndoableAction {
      private final UndoableAction[] subActions_;
      private UndoableAction lastAction_ = null;
      private final MersenneTwisterFast random_;
      public Multi(UndoableAction[] subActions) {
        this.subActions_ = subActions;
        this.random_ = new MersenneTwisterFast();
      }
      public double doAction(double currentScore, double desparationValue) {
        lastAction_ = subActions_[random_.nextInt(subActions_.length)];
        return lastAction_.doAction(currentScore,desparationValue);
      }
      public boolean isActionSuccessful() {
        if(lastAction_!=null) {
          return lastAction_.isActionSuccessful();
        }
        throw new RuntimeException("Assertion error : isActionSuccessful() called when no action has been done recently");
      }
      /**
       * @return false
       */
      public boolean isActionDeterministic() {
        return false;
      }
      public boolean undoAction() {
        if(lastAction_!=null) {
          final boolean successful = lastAction_.undoAction();
          lastAction_ = null;
          return successful;
        } else {
          throw new RuntimeException("Assertion error : undoAction() called when no action has been done recently (or has already been undone)");
        }
      }
      public void printStats(PrintWriter out) {};
      public String getDescription() {
    	  return "Multi(" + lastAction_.getDescription() + ")";
      }
    } //End of class Multi
// -=-==-=--==--=-=-=-=-=-=-==--=-=
    private static class DistributedMulti implements UndoableAction {
      protected final UndoableAction[] subActions_;
      protected final double[] probabilities_;
      protected UndoableAction lastAction_ = null;
      protected final MersenneTwisterFast random_;

      public DistributedMulti(UndoableAction[] subActions, double[] proportions) {
    	  this(subActions,proportions,0);
      }
      public DistributedMulti(UndoableAction[] subActions, double[] proportions, long seed) {
        this.subActions_ = subActions;
        this.probabilities_ = new double[subActions.length];
        double total = 0;
        for(int i = 0 ; i < subActions.length ; i++) {
          total+=proportions[i];
        }
        for(int i = 0 ; i < subActions.length ; i++) {
          probabilities_[i] = proportions[i]/total;
        }
        if (seed == 0) {
        	this.random_ = new MersenneTwisterFast();
        } else {
        	this.random_ = new MersenneTwisterFast(seed);
        }
      }
      /**
       * @return false
       */
      public boolean isActionDeterministic() {
        return false;
      }
      public double doAction(double currentScore, double desparationValue) {
        double v = random_.nextDouble();
        double total = 0;
        int index = subActions_.length-1;
        for(int i = 0 ; i < subActions_.length ; i++) {
          total+=probabilities_[i];
          if(total>v) {
            index = i;
            break;
          }
        }
        lastAction_ = subActions_[index];
        return lastAction_.doAction(currentScore,desparationValue);
      }
      public boolean isActionSuccessful() {
        if(lastAction_!=null) {
          return lastAction_.isActionSuccessful();
        }
        throw new RuntimeException("Assertion error : isActionSuccessful() called when no action has been done recently");
      }
      public boolean undoAction() {
        if(lastAction_!=null) {
          boolean successful = lastAction_.undoAction();
          lastAction_ = null;
          return successful;
        } else {
          throw new RuntimeException("Assertion error : undoAction() called when no action has been done recently (or has already been undone)");
        }
      }
      public void printStats(PrintWriter out) {};
      public String getDescription() {
    	  return "DistributedMulti(" + lastAction_.getDescription() + ")";
      }
    } //End of class DistributedMulti
 // -=-==-=--==--=-=-=-=-=-=-==--=-=
    // Adds gathering of statistics to ShuffledMulti
    private static class ProfiledMulti extends ShuffledMulti {
 
    	// inner class to accumulate stats on results of actions
    	protected class ActionStats {
    		public double likelihoodChange;
    		public long timeTaken;
    		public boolean undone;
    		
    		public ActionStats(double lChange, long time) {
    			likelihoodChange = lChange;
    			timeTaken = time;
    			undone = false;
    		}
    	}
    	
    	Vector<ActionStats>[] stats_ = null;
    	ActionStats lastStats_ = null;
    	
      public ProfiledMulti(UndoableAction[] subActions, int[] proportions) {
    	  this(subActions,proportions,0);
      }
      public ProfiledMulti(UndoableAction[] subActions, int[] proportions, long seed) {
    	  super(subActions,proportions,seed);
    	  int n = subActions.length;
    	  stats_ = new Vector[n];
    	  for (int i=0; i<n; i++) {
    		  stats_[i] = new Vector<ActionStats>();
    	  }
      }
 
      public double doAction(double currentScore, double desparationValue) {
    	int index = nextActionIndex();
    	lastAction_ = subActions_[index];
        long startTime = System.currentTimeMillis();
        double result = lastAction_.doAction(currentScore,desparationValue);
        // On very first call, currentScore may be infinite, which of course would mess up our stats.
        if (!Double.isInfinite(currentScore)) {
        	lastStats_ = new ActionStats(result-currentScore,System.currentTimeMillis()-startTime);
        	stats_[index].add(lastStats_);
        }
        return result;
      }
      
      public boolean undoAction() {
    	  boolean success = super.undoAction();
    	  if (success) {
    		  lastStats_.undone = true;
    	  }
    	  return success;
      }
      public void printStats(PrintWriter out) {
    	  final double GOOD_THRESHOLD = 0.01;
    	  final double VERY_GOOD_THRESHOLD = 0.1;
    	  out.println("Statistics for ProfiledMultiAction:");
    	  for (int i=0; i<subActions_.length; i++) {
    		  int nAttempted = stats_[i].size();
    		  int nSucceeded = 0;
    		  for (int j=0; j<nAttempted; j++) {
    			  if (!stats_[i].elementAt(j).undone) {
    				  nSucceeded++;
    			  }
    		  }
    		  out.printf("Action %s succeeded %d times and was undone %d times\n", subActions_[i].getClass().getName(), nSucceeded, nAttempted-nSucceeded);
    		  long times[] = new long[nSucceeded];
    		  double delta[] = new double[nSucceeded];
    		  long sumTimes = 0;
    		  double sumDelta = 0;
    		  int index = 0;
    		  int countGood = 0;
    		  int countVeryGood = 0;
    		  for (int j=0; j<nAttempted; j++) {
    			  if (!stats_[i].elementAt(j).undone) {
    				  times[index]=stats_[i].elementAt(j).timeTaken;
    				  sumTimes += times[index];
    				  double change = stats_[i].elementAt(j).likelihoodChange;
    				  delta[index]=change;
    				  sumDelta += change;
    				  countGood += (change > GOOD_THRESHOLD && change <= VERY_GOOD_THRESHOLD) ? 1 : 0;
    				  countVeryGood += (change > VERY_GOOD_THRESHOLD) ? 1 : 0;
    				  index++;
    			  }
    		  }
    		  Arrays.sort(times);
    		  Arrays.sort(delta);
    		  if (nSucceeded > 0) {
    			  // some truncating of integers happens here
    			  out.printf("Times (ms): mean = %d, median = %d\n", sumTimes/nSucceeded, times[nSucceeded/2]);
    			  out.printf("Likelihood change: mean = %f, median = %f, %d good results, %d very good results of %d\n", 
    					  sumDelta/nSucceeded, delta[nSucceeded/2], countGood, countVeryGood, nSucceeded);
    		  }
    	  }
      };
      public String getDescription() {
    	  return "ProfiledMulti(" + lastAction_.getDescription() + ")";
      }
    } //End of class ProfiledMulti
 // -=-==-=--==--=-=-=-=-=-=-==--=-=
    // ShuffledMulti: has a list of actions  to perform. It shuffles the list, and each time an action is
    // requested, it takes the next one. Resuffles when through the list.
    // This is a somewhat more deterministic alternative to DistributedMulti: no action can (by chance) never
    // be executed.
    private static class ShuffledMulti implements UndoableAction {
      protected final UndoableAction[] subActions_;
      protected UndoableAction lastAction_ = null;
      protected final MersenneTwisterFast random_;
      private int[] shuffled_;
      private int shuffleIndex_;

      public ShuffledMulti(UndoableAction[] subActions, int[] proportions) {
    	  this(subActions,proportions,0);
      }
      public ShuffledMulti(UndoableAction[] subActions, int[] proportions, long seed) {
        this.subActions_ = subActions;
        int n = MathUtils.getTotal(proportions);
        shuffled_ = new int[n];
        int index = 0;
        for(int i = 0 ; i < subActions.length ; i++) {
        	for (int j=0; j<proportions[i]; j++) {
        		shuffled_[index++] = i;
        	}
        }
        if (seed == 0) {
        	this.random_ = new MersenneTwisterFast();
        } else {
        	this.random_ = new MersenneTwisterFast(seed);
        }
        reshuffle();
      }
        
      private void reshuffle() {
    	  shuffleIndex_ = 0;
    	  random_.shuffle(shuffled_);
      }
      protected int nextActionIndex() {
    	  if (shuffleIndex_ >= shuffled_.length) reshuffle();
    	  return shuffled_[shuffleIndex_++];
      }
      
      /**
       * @return false
       */
      public boolean isActionDeterministic() {
        return false;
      }
      public double doAction(double currentScore, double desparationValue) {
        lastAction_ = subActions_[nextActionIndex()];
        return lastAction_.doAction(currentScore,desparationValue);
      }
      public boolean isActionSuccessful() {
        if(lastAction_!=null) {
          return lastAction_.isActionSuccessful();
        }
        throw new RuntimeException("Assertion error : isActionSuccessful() called when no action has been done recently");
      }
      public boolean undoAction() {
        if(lastAction_!=null) {
          boolean successful = lastAction_.undoAction();
          lastAction_ = null;
          return successful;
        } else {
          throw new RuntimeException("Assertion error : undoAction() called when no action has been done recently (or has already been undone)");
        }
      }
      public void printStats(PrintWriter out) {};
      public String getDescription() {
    	  return "ShuffledMulti(" + lastAction_.getDescription() + ")";
      }
    } //End of class ShuffledMulti
 // -=-==-=--==--=-=-=-=-=-=-==--=-=
    private static class Combined implements UndoableAction {
      private final UndoableAction[] subActions_;
      private boolean deterministic_ = true;
      private boolean successful_ = false;
      public Combined(UndoableAction[] subActions) {
        this.subActions_ = subActions;
      }
      /**
       * @return false
       */
      public boolean isActionDeterministic() {  return deterministic_;     }
      public double doAction(double currentScore, double desparationValue) {
        boolean d = true;
				boolean s = true;
				for(int i = 0 ; i < subActions_.length ; i++) {
					UndoableAction a = subActions_[i];
					double score = a.doAction(currentScore, desparationValue);
					if(a.isActionSuccessful()) {
						s = true;
						currentScore = score;
				    d = d & a.isActionDeterministic();
					}
				}
				deterministic_ = d;
				successful_ = s;
				return currentScore;
      }
      public boolean isActionSuccessful() { return successful_; }
      public boolean undoAction() {
				boolean result = true;
				if(successful_) {
			    for(int i = subActions_.length -1 ; i >= 0 ; i--) {
					  UndoableAction a = subActions_[i];
					  if(a.isActionSuccessful()) {
					    result = result & a.undoAction();
					  }
				  }
					successful_ = false;
      		return result;
        } else {
          throw new RuntimeException("Assertion error : undoAction() called when not successful");
        }
		  }
      public void printStats(PrintWriter out) {};
      public String getDescription() {
    	  StringBuffer description = new StringBuffer("Combined(");
    	  for (UndoableAction subAction : subActions_) {
    		  description.append(subAction.getDescription()).append(',');
    	  }
    	  description.replace(description.length(), description.length(), ")"); // replace final ','
    	  return description.toString();
      }
    } //End of class Combined
  } //End of class Utils
}