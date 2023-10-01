This is an updated version of the original code V2.0 //  MOSFP.java
// Add this code to jMetal Tool: http://jmetal.sourceforge.net/
/*
Shehadeh HA, Idna Idris MY, Ahmedy I, Ramli R, Mohamed Noor N. The Multi-Objective Optimization Algorithm Based on Sperm Fertilization Procedure (MOSFP) Method for Solving Wireless Sensor Networks Optimization Problems in Smart Grid Applications. Energies. 2018; 11(1):97. https://doi.org/10.3390/en11010097 


Shehadeh HA, Ldris MYI, Ahmedy I. Multi-Objective Optimization Algorithm Based on Sperm Fertilization Procedure (MOSFP). Symmetry. 2017; 9(10):241. https://doi.org/10.3390/sym9100241

*/       



package jmetal.metaheuristics.MOSFP;


import jmetal.core.Algorithm;

import jmetal.core.Problem;

import jmetal.core.SolutionSet;

import jmetal.operators.mutation.Mutation;

import jmetal.operators.mutation.NonUniformMutation;

import jmetal.operators.mutation.UniformMutation;

import jmetal.problems.Kursawe;


import jmetal.problems.ProblemFactory;

import jmetal.qualityIndicator.QualityIndicator;

import jmetal.util.Configuration;

import jmetal.util.JMException;


import java.io.IOException;

import java.util.HashMap;

import java.util.logging.FileHandler;

import java.util.logging.Logger;

import jmetal.problems.DTLZ.DTLZ1;

import jmetal.problems.ZDT.ZDT1;

import jmetal.problems.ZDT.ZDT2;

import jmetal.problems.ZDT.ZDT3;

import jmetal.problems.ZDT.ZDT4;

import jmetal.problems.ZDT.ZDT6;

import jmetal.problems.singleObjective.Rastrigin;

import jmetal.problems.singleObjective.Sphere;

import jmetal.problems.singleObjective.Rosenbrock;



/**
 * Class for configuring and running the MOSFP algorithm
 */


public class MOSFP_main {
  


public static Logger      logger_ ;      // Logger object
  public static FileHandler fileHandler_ ; // FileHandler object



   public static void main(String [] args) throws JMException, IOException, ClassNotFoundException {
    Problem   problem   ;         // The problem to solve
    Algorithm algorithm ;         // The algorithm to use
    Mutation  uniformMutation ;
    Mutation nonUniformMutation ;
    
    QualityIndicator indicators ; // Object to get quality indicators
        
    HashMap  parameters ; // Operator parameters
    
    // Logger object and file to store log messages
    logger_      = Configuration.logger_ ;
    fileHandler_ = new FileHandler("MOSFP_main.log"); 
    logger_.addHandler(fileHandler_) ;

indicators = null ;
    if (args.length == 1) {
      Object [] params = {"Real"};
      problem = (new ProblemFactory()).getProblem(args[0],params);
    } // if
    else if (args.length == 2) {
      Object [] params = {"Real"};
      problem = (new ProblemFactory()).getProblem(args[0],params);
      indicators = new QualityIndicator(problem, args[1]) ;
    } // if
    else { // Default problem
      //problem = new Kursawe("Real", 3); 
      //problem = new Water("Real");
      // problem = new ZDT4("Real");
      //problem = new WFG1("Real");
      //problem = new DTLZ1("Real");
      //problem = new OKA2("Real") ;
      // problem = new ZDT1("Real");
    }
	}

algorithm = new MOSFP(problem) ;
    
    Integer maxIterations = 250 ;
    Double perturbationIndex = 0.5 ;
    Double mutationProbability = 1.0/problem.getNumberOfVariables() ;
    
    // Algorithm parameters
    algorithm.setInputParameter("swarmSize",20);
    algorithm.setInputParameter("archiveSize",20);
    algorithm.setInputParameter("maxIterations",maxIterations);


 parameters = new HashMap() ;
    parameters.put("probability", mutationProbability) ;
    parameters.put("perturbation", perturbationIndex) ;
    uniformMutation = new UniformMutation(parameters);
    
    parameters = new HashMap() ;
    parameters.put("probability", mutationProbability) ;
    parameters.put("perturbation", perturbationIndex) ;
    parameters.put("maxIterations", maxIterations) ;
    nonUniformMutation = new NonUniformMutation(parameters);

    // Add the operators to the algorithm
    algorithm.addOperator("uniformMutation",uniformMutation);
    algorithm.addOperator("nonUniformMutation",nonUniformMutation);

    // Execute the Algorithm 
    long initTime = System.currentTimeMillis();
    SolutionSet population = algorithm.execute();
    long estimatedTime = System.currentTimeMillis() - initTime;

// Print the results
    logger_.info("Total execution time: "+estimatedTime + "ms");
    logger_.info("Variables values have been writen to file VAR");
    population.printVariablesToFile("VAR");    
    logger_.info("Objectives values have been writen to file FUN");
    population.printObjectivesToFile("FUN");
  
    if (indicators != null) {
      logger_.info("Quality indicators") ;
      logger_.info("Hypervolume: " + indicators.getHypervolume(population)) ;
      logger_.info("GD         : " + indicators.getGD(population)) ;
      logger_.info("IGD        : " + indicators.getIGD(population)) ;
      logger_.info("Spread     : " + indicators.getSpread(population)) ;
      logger_.info("Epsilon    : " + indicators.getEpsilon(population)) ;  
    } // if
  }//main
} // MOSFP_main


