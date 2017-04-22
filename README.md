# BLEAQ2
Version uploaded on 22nd April 2017

Readme file for Bilevel Evolutionary Algorithm based on Quadratic Approximations - Version 2


Quick instruction for execution:
------------------------------------------------------------------------
Code the upper level optimization task in ulExternalProblem.m

Code the lower level optimization task in llExternalProblem.m

xu and xl are the upper and lower level decision vectors respectively

Provide the problem and algorithm parameters in externalProblem.m

Execute it as: externalProblem()

The code is written for maximization at both levels. Feasibility is given as g(x)<=0 and h(x)=0.

The code also allows handling of multiple objectives at upper level. Check the sample implementation in externalProblemMulti.m


BLEAQ2
------------------------------------------------------------------------
BLEAQ2 is the second version of a computationally efficient evolutionary algorithm for non-linear bilevel optimization problems. More information about the working of the algorithm can be found from the following paper.

Sinha, Ankur, Zhichao Lu, Kalyanmoy Deb, and Pekka Malo, "Bilevel Optimization based on Iterative Approximation of Mappings." arXiv preprint arXiv:1702.03394 (2017).

The current implementation is built on earlier versions and also supports solving bilevel optimization problems with multiple objectives at the upper level. Following papers provide further information about the algorithm.

BLEAQ Papers (First version)

Sinha, Ankur, Pekka Malo, and Kalyanmoy Deb. "Evolutionary algorithm for bilevel optimization using approximations of the lower level optimal solution mapping." European Journal of Operational Research 257.2 (2017): 395-411.

Sinha, Ankur, Pekka Malo, and Kalyanmoy Deb. "Efficient evolutionary algorithm for single-objective bilevel optimization." arXiv preprint arXiv:1303.3901 (2013).

Sinha, Ankur, Pekka Malo, and Kalyanmoy Deb. "An improved bilevel evolutionary algorithm based on Quadratic Approximations." In 2014 IEEE Congress on Evolutionary Computation, 2014.

m-BLEAQ Papers (For multiobjective bilevel problems)

Sinha, Ankur, Pekka Malo, Kalyanmoy Deb, Pekka Korhonen, and Jyrki Wallenius. "Solving bilevel multicriterion optimization problems with lower level decision uncertainty." IEEE Transactions on Evolutionary Computation 20, no. 2 (2016): 199-217.

Sinha, Ankur, Pekka Malo, and Kalyanmoy Deb. "Towards understanding bilevel multi-objective optimization with deterministic lower level decisions." International Conference on Evolutionary Multi-Criterion Optimization. Springer International Publishing, 2015.


Files in the package
------------------------------------------------------------------------
There are the following Matlab (.m) files in this package:

ulSearch.m: Performs search at the upper level.

llSearch.m: Performs search at the lower level.

quadApprox.m: Supporting file for quadratic creating quadratic approximations.

ulTestProblem.m: Source code for upper level SMD and TP suite.

llTestProblem.m: Source code for lower level SMD and TP suite.

ulExternalProblem.m: Source code for upper level optimization task for a user defined problem.

ulExternalProblem.m: Source code for lower level optimization task for a user defined problem

terminationCheck.m: Source code for the termination criteria used in the algorithm. This can be modified based on the user's requirements.

smd1.m - smd14.m: These files contain the problem and algorithm parameters for SMD-Suite. The implementation for the test problems can be found in the ulTestProblem.m and llTestProblem.m files.

tp1.m - tp10.m: These files contain the problem and algorithm parameters for TP-Suite. The implementation for the test problems can be found in the ulTestProblem.m and llTestProblem.m files.

externalProblem.m: It contains the problem and algorithm parameters for a user defined problem.

msmd1: This file contains the problem and algorithm parameters for a sample multiobjective bilevel optimization problem. The implementation for the test problem can be found in the ulTestProblem.m and llTestProblem.m files.

Following are other supporting files
------------------------------------------------------------------------
calculateCrowdingDistance.m

nonDominatedSorting.m

getMappings.m

getLowerLevelVariableFromMapping.m

getOptimalSolutionSMD.m


Executing a user-defined problem
------------------------------------------------------------------------
To execute a user-defined problem, code the upper level optimization task in ulExternalProblem.m and the lower level optimization task in llExternalProblem.m. The functions inside the files contain arguments xu and xl, which represent the upper level decision vector and lower level decision vector respectively.

Provide the problem and algorithm parameters for user-defined problem in externalProblem.m and call the following command to execute the user-defined bilevel optimization task.
externalProblem()

The results of the execution are printed on the screen as well as stored in 'externalProblem.mat'. A sample bilevel optimization problem is already coded as an external problem.

Problem parameters to be defined in externalProblem.m

ulDim: Number of dimensions at upper level.

llDim: Number of dimensions at lower level.

ulDimMax: Vector defining maximum values for upper level variables.

ulDimMin: Vector defining minimum values for upper level variables.

llDimMax: Vector defining maximum values for lower level variables.

llDimMin: Vector defining minimum values for lower level variables.

Algorithm parameters to be defined in externalProblem.m

ulPopSize: Upper level population size.

llPopSize: Lower level population size.

ulMaxGens: Upper level maximum generations.

llMaxGens: Lower level maximum generations.

ulStoppingCriteria: Stopping parameter at upper level. Smaller the value higher the accuracy.

llStoppingCriteria: Stopping parameter at lower level. Smaller the value higher the accuracy.

There are other parameters in the algorithm. However, most of them are either adaptive or not necessary to be adjusted.

Output of the execution

ulEliteFunctionValue: Upper level function value for the elite member.

llEliteFunctionValue: Lower level function value for the elite member.

ulEliteIndiv: Upper level elite member.

llEliteIndiv: Lower level elite member.

ulFunctionEvaluations: Upper level function evaluations required during the exection.

llFunctionEvaluations: Lower level function evaluations required during the exection.


Executing the SMD or TP suite
------------------------------------------------------------------------
To execute one of the SMD test problems (say SMD1), the following command needs to be called:

smd1()

This executes the SMD1 test problem with problem and algorithm parameters coded in smd1.m. It executes a 5 variable SMD1 problem with 2 upper level variables and 3 lower level variables. The results are printed on the screen as well as stored in 'smd1.mat'

To execute one of the problems in TP-Suite (say SMD1), the following command needs to be called:

tp1()

This executes tp1 with problem and algorithm parameters coded in tp1.m. The results are printed on the screen as well as stored in 'tp1.mat'


Contact
------------------------------------------------------------------------
In case you have any questions, comments, suggestions, or you want to report any bugs, you can send an email to Ankur Sinha (asinha@iima.ac.in)

Ankur Sinha, PhD

Indian Institute of Management

Ahmedabad, India

asinha@iima.ac.in
