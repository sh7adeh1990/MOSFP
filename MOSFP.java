//  MOSFP.java
// Add this code to jMetal Tool: http://jmetal.sourceforge.net/
/*
Shehadeh HA, Idna Idris MY, Ahmedy I, Ramli R, Mohamed Noor N. The Multi-Objective Optimization Algorithm Based on Sperm Fertilization Procedure (MOSFP) Method for Solving Wireless Sensor Networks Optimization Problems in Smart Grid Applications. Energies. 2018; 11(1):97. https://doi.org/10.3390/en11010097 


Shehadeh HA, Ldris MYI, Ahmedy I. Multi-Objective Optimization Algorithm Based on Sperm Fertilization Procedure (MOSFP). Symmetry. 2017; 9(10):241. https://doi.org/10.3390/sym9100241

*/

package jmetal.metaheuristics.MOSFP;

import jmetal.core.*;
import jmetal.operators.mutation.Mutation;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.NonDominatedSolutionList;
import jmetal.util.PseudoRandom;
import jmetal.util.archive.CrowdingArchive;
import jmetal.util.comparators.CrowdingDistanceComparator;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.comparators.EpsilonDominanceComparator;

import java.util.Comparator;

/**
 * This class representing an asynchronous version of MOSFP algorithm
 */
public class MOSFP extends Algorithm {
                 
  /**
   * Stores the number of sperms_ used
   */
  private int spermsSize_;
  
  /**
  * Stores the maximum size for the archive
  */
  private int archiveSize_;
  
  /**
  * Stores the maximum number of iteration_
  */
  private int maxIterations_;
  
  /**
  * Stores the current number of iteration_
  */
  private int iteration_;
  
  /**
  * Stores the perturbation used by the non-uniform mutation
  */
  private double perturbation_;
  
  /**
  * Stores the sperms
  */
  private SolutionSet sperms_;
  
  /**
   * Stores the best_ solutions founds so far for each sperms
   */
  private Solution[] best_;
  
  /**
  * Stores the winners_
  */
  private CrowdingArchive winners_ ;
  
  /**
  * Stores the epsilon-archive
  */
  private NonDominatedSolutionList eArchive_;
  
  /**
  * Stores the speed_ of each sperm
  */
  private double [][] speed_;  
  
  /**
  * Stores a comparator for checking dominance
  */
  private Comparator dominance_;
  
  /**
  * Stores a comparator for crowding checking
  */
  private Comparator crowdingDistanceComparator_;
  
  /**
   * Stores a <code>Distance</code> object
   */
  private Distance distance_;
  
  /**
  * Stores a operator for uniform mutations
  */
  private Operator uniformMutation_;
  
  /**
  * Stores a operator for non uniform mutations
  */ 
  private Operator nonUniformMutation_;
  
  /**
  * eta_ value
  */
  private double eta_ = 0.0075;

  /** 
  * Constructor
  * @param problem Problem to solve
  */    
  public MOSFP(Problem problem) {                
    super (problem) ;
  } // MOSFP
  
  /**
   * Initialize all parameter of the algorithm
   */
  public void initParams(){
    spermsSize_ = ((Integer)getInputParameter("swarmSize")).intValue();
    archiveSize_   = ((Integer)getInputParameter("archiveSize")).intValue();
    maxIterations_ = ((Integer)getInputParameter("maxIterations")).intValue();

    
    sperms_     = new SolutionSet(spermsSize_);        
    best_          = new Solution[spermsSize_];
    winners_       = new CrowdingArchive(archiveSize_,problem_.getNumberOfObjectives());
    eArchive_      = new NonDominatedSolutionList(new EpsilonDominanceComparator(eta_));
    
    uniformMutation_ = (Mutation)operators_.get("uniformMutation") ;
    nonUniformMutation_ = (Mutation)operators_.get("nonUniformMutation") ;
    
    // Create the dominator for equadless and dominance
    dominance_          = new DominanceComparator();    
    crowdingDistanceComparator_ = new CrowdingDistanceComparator();
    distance_           = new Distance();
    
    // Create the speed_ vector
    speed_ = new double[spermsSize_][problem_.getNumberOfVariables()];
  } // initParams
           
  
  /**
   * Update the spped of each sperm
   * @throws JMException 
   */
  private void computeSpeed() throws JMException{        
    double r1,r2,W,P0,P1,P2,T1,T2; 
    Variable[] bestGlobal;                                            
        
    for (int i = 0; i < spermsSize_; i++){
      Variable[] sperm     = sperms_.get(i).getDecisionVariables();
      Variable[] bestSperm = best_[i].getDecisionVariables();                        

      //Select a global best_ for calculate the speed of sperm i, bestGlobal
      Solution one, two;
      int pos1 = PseudoRandom.randInt(0,winners_.size()-1);
      int pos2 = PseudoRandom.randInt(0,winners_.size()-1);
      one = winners_.get(pos1);
      two = winners_.get(pos2);

      if (crowdingDistanceComparator_.compare(one,two) < 1)
        bestGlobal = one.getDecisionVariables();
      else
        bestGlobal = two.getDecisionVariables();
      //
            
      //Params for velocity equation
      r1 = PseudoRandom.randDouble();
      r2 = PseudoRandom.randDouble();
      W  = 0.55;
      P0 = PseudoRandom.randDouble(7.0,14.0);
      P1 = PseudoRandom.randDouble(7.0,14.0);
      P2 = PseudoRandom.randDouble(7.0,14.0);
      T1 = PseudoRandom.randDouble(35.1,38.5);
      T2 = PseudoRandom.randDouble(35.1,38.5);
      //

      for (int var = 0; var < sperm.length; var++){                                     
        //Computing the velocity of this sperm
        speed_[i][var] =  r1 * Math.log10(P0)* speed_[i][var] +
                   Math.log10(P1) * Math.log10(T1)* (bestSperm[var].getValue() - 
                              sperm[var].getValue()) +
                   Math.log10(P2) * Math.log10(T2) * (bestGlobal[var].getValue() - 
                              sperm[var].getValue());
      }
                
    }
  } // computeSpeed
     
  /**
   * Update the position of each sperm
   * @throws JMException 
   */
  private void computeNewPositions() throws JMException{
    for (int i = 0; i < spermsSize_; i++){
    	Variable[] sperm = sperms_.get(i).getDecisionVariables();
      //sperm.move(speed_[i]);
      for (int var = 0; var < sperm.length; var++){
        sperm[var].setValue(sperm[var].getValue()+ speed_[i][var]);
        if (sperm[var].getValue() < problem_.getLowerLimit(var)){
          sperm[var].setValue(problem_.getLowerLimit(var));                    
          speed_[i][var] = speed_[i][var] * -1.0;    
        }
        if (sperm[var].getValue() > problem_.getUpperLimit(var)){
          sperm[var].setValue(problem_.getUpperLimit(var));                    
          speed_[i][var] = speed_[i][var] * -1.0;    
        }                                             
      }
    }
  } // computeNewPositions
        
   
  /**
   * Apply a mutation operator to all sperms in the swarm
   * @throws JMException 
   */
  private void mosfpMutation(int actualIteration, int totalIterations) throws JMException{       
    //There are three groups of sperms_, the ones that are mutated with
    //a non-uniform mutation operator, the ones that are mutated with a 
    //uniform mutation and the one that no are mutated
    nonUniformMutation_.setParameter("currentIteration",actualIteration);
    //*/

    for (int i = 0; i < sperms_.size();i++)            
      if (i % 3 == 0) { //sperms_ mutated with a non-uniform mutation
        nonUniformMutation_.execute(sperms_.get(i));                                
      } else if (i % 3 == 1) { //sperms_ mutated with a uniform mutation operator
        uniformMutation_.execute(sperms_.get(i));                
      } else //sperms_ without mutation
          ;      
  } // mosfpMutation
   
    
  /**   
  * Runs of the MOSFP algorithm.
  * @return a <code>SolutionSet</code> that is a set of non dominated solutions
  * as a result of the algorithm execution  
   * @throws JMException 
  */  
  public SolutionSet execute() throws JMException, ClassNotFoundException {
    initParams();

    //->Step 1 (and 3) Create the initial population and evaluate
    for (int i = 0; i < spermsSize_; i++){
      Solution sperm = new Solution(problem_);
      problem_.evaluate(sperm);
      problem_.evaluateConstraints(sperm);
      sperms_.add(sperm);                   
    }
        
    //-> Step2. Initialize the speed_ of each sperm to 0
    for (int i = 0; i < spermsSize_; i++) {
      for (int j = 0; j < problem_.getNumberOfVariables(); j++) {
        speed_[i][j] = 0.0;
      }
    }
    
        
    // Step4 and 5   
    for (int i = 0; i < sperms_.size(); i++){
      Solution sperm = new Solution(sperms_.get(i));            
      if (winners_.add(sperm)){
        eArchive_.add(new Solution(sperm));
      }
    }
                
    //-> Step 6. Initialice the memory of each sperm
    for (int i = 0; i < sperms_.size(); i++){
      Solution sperm = new Solution(sperms_.get(i));           
      best_[i] = sperm;
    }
        
    //Crowding the winners_
    distance_.crowdingDistanceAssignment(winners_,problem_.getNumberOfObjectives());        

    //-> Step 7. Iterations ..        
    while (iteration_ < maxIterations_){
      //Compute the speed_        
      computeSpeed();
            
      //Compute the new positions for the sperms_            
      computeNewPositions();

      //Mutate the sperms_          
      mosfpMutation(iteration_,maxIterations_);                       
            
      //Evaluate the new sperms_ in new positions
      for (int i = 0; i < sperms_.size(); i++){
        Solution sperm = sperms_.get(i);
        problem_.evaluate(sperm);                
        problem_.evaluateConstraints(sperm);                
      }
            
      //Actualize the archive          
      for (int i = 0; i < sperms_.size(); i++){
        Solution sperm = new Solution(sperms_.get(i));                
        if (winners_.add(sperm)){
          eArchive_.add(new Solution(sperm));
        }                
      }
            
      //Actualize the memory of this sperm
      for (int i = 0; i < sperms_.size();i++){
        int flag = dominance_.compare(sperms_.get(i),best_[i]);
        if (flag != 1) { // the new sperm is best_ than the older remeber        
          Solution sperm = new Solution(sperms_.get(i));                    
          //this.best_.reemplace(i,sperm);
          best_[i] = sperm;
        }
      }       
            
      //Crowding the winners_
      distance_.crowdingDistanceAssignment(winners_,
                                              problem_.getNumberOfObjectives());            
      iteration_++;
    }
        
    return this.winners_;
    //return eArchive_;
  } // execute
    
  /** 
   * Gets the winners of the MOSFP algorithm
   */
  public SolutionSet getwinner(){
    return winners_;
  }  // getwinner 
} // MOSFP