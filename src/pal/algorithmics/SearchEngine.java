// SearchEngine.java
//
// (c) 1999-2003 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package pal.algorithmics;

import pal.math.MersenneTwisterFast;
import pal.util.*;

import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * A simplistic class (most of the work is done elsewhere) that handles basic search algorithms
 *
 * @version $Id: SearchEngine.java,v 1.2 2003/10/19 02:35:26 matt Exp $
 *
 * @author Matthew Goode
 */

public class SearchEngine {
  private final ProbabilityIterator.Factory probabilityIteratorFactory_;
  private double lastScore_ = Double.NEGATIVE_INFINITY;
  private MersenneTwisterFast random_;
  private static final SimpleDateFormat formatter = new SimpleDateFormat("hh:mm:ss a"); // time format, for debugging output
  private int printSubjectEveryNScoreChanges_ = 0; // Non-essential, does nothing if equal to zero, nasty hacked-in addition. 

  public SearchEngine(ProbabilityIterator.Factory probabilityIteratorFactory, MersenneTwisterFast random) {
	  this.probabilityIteratorFactory_ = probabilityIteratorFactory;
	  this.random_ = random;
  }
  public SearchEngine( ProbabilityIterator.Factory probabilityIteratorFactory) {
	  this(probabilityIteratorFactory, new MersenneTwisterFast());
  }
  
  public SearchEngine(ProbabilityIterator.Factory probabilityIteratorFactory, long seed) {
	  this(probabilityIteratorFactory, new MersenneTwisterFast(seed));
  }
  
  public double getLastScore() { return lastScore_; }
  
  public void setPrintSubjectEveryNScoreChanges(int n) { printSubjectEveryNScoreChanges_ = n; };
  
  public void run(AlgorithmCallback callback, final double initialScore, ObjectState subject, StoppingCriteria.Factory stoppingCriteria, Ranker ranker) {
	  run(callback,initialScore,subject,stoppingCriteria,ranker,null);
  }
  public void run(AlgorithmCallback callback, 
		  final double              initialScore, 
		  ObjectState               subject, 
		  StoppingCriteria.Factory  stoppingCriteria, 
		  Ranker                    ranker, 
		  PrintWriter               debug) 
  {
      double score = initialScore;
      StoppingCriteria stopper = stoppingCriteria.newInstance();
      ProbabilityIterator acceptanceProbability = probabilityIteratorFactory_.newInstance();
      int evaluationCount = 0;
      int scoreChangeCount = 0; // for printSubjectEveryNScoreChanges_
      final boolean maximising = subject.isMaximiseScore();
      while(!stopper.isTimeToStop()) {
        double newScore = subject.doAction(score,stopper.getRelativeStoppingRatio());
        evaluationCount++;
        double probability = acceptanceProbability.getNextProbability(score, newScore, maximising);
        if(ranker.isWorthAdding(newScore,maximising)) {
          ranker.add(subject.getStateReference(), newScore,maximising);
        }
        if ( (!maximising&&(newScore<=score)) || (maximising&&(newScore>=score)) || probability==1.0 || random_.nextDouble()<probability ) {
          score = newScore;
          scoreChangeCount++;
          if (debug != null) {
        	  debug.printf("Searcher: %s eval# %5d score %f (p=%5.3f) action=%s\n", 
        			  formatter.format(new Date()),evaluationCount,score,probability,subject.getActionDescription());
        	  
          }
          if (printSubjectEveryNScoreChanges_ !=0 && (scoreChangeCount % printSubjectEveryNScoreChanges_)==0) {
        	  if (debug != null) {
        		  debug.println(subject.toString());
        	  } else {
        		  System.out.printf("After %d score changes: %s\n",scoreChangeCount,subject.toString());
        	  }
          }
        } else {
          if(!subject.undoAction())  {
            //Undo was unsuccessful so we have to stick to what we have even if it was worse!
            score = newScore;
          }
        }
        stopper.newIteration(score,ranker.getBestScore(),!maximising,acceptanceProbability.isStablised(), callback);
        if(callback.isPleaseStop()) { break; }
      }
      lastScore_ = score;
      if (debug != null) { debug.printf("Searcher: %s eval# %5d score %f **EXIT**\n", formatter.format(new Date()),evaluationCount,score); }
      callback.updateStatus("Finished:"+score);
      //System.out.println("Evaluation Count:"+(evaluationCount+1));
    }

}

