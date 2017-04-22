function externalProblem()

problemName = 'externalProblem';             % Test problem name

ulPopSize=100;                    % Size of UL population
ulMaxGens=1000;                  % Maximum number of generations allowed at UL
ulDim=5;                         % Number of UL dimensions

llPopSize=100;                    % Size of LL population
llMaxGens=1000;                  % Maximum number of generations allowed at LL
llDim=5;                         % Number of LL dimensions

ulDimMin = [0 0 0 0 0];            % Minimum value accross UL dimensions
ulDimMax = [20 20 20 20 20];              % Maximum value accross UL dimensions

llDimMin = [0 0 0 0 0];      % Minimum value accross LL dimensions
llDimMax = [20 20 20 20 20];     % Maximum value accross LL dimensions

ulStoppingCriteria = 1e-5;
llStoppingCriteria = 1e-5;

[ulEliteFunctionValue, llEliteFunctionValue, ulEliteIndiv, llEliteIndiv, ulFunctionEvaluations, llFunctionEvaluations]=ulSearch(problemName, ulPopSize, ulMaxGens, ulDim, ulDimMin, ulDimMax, llPopSize, llMaxGens, llDim, llDimMin, llDimMax, ulStoppingCriteria, llStoppingCriteria)

save('externalProblem');


%Solution

%ulEliteFunctionValue =   1684.0966
%llEliteFunctionValue =  490.7658
%ulEliteIndiv =   12.0163    9.3334   15.6666   10.6250    9.3750
%llEliteIndiv =    8.8683    6.1317    9.8683    5.0000    2.6562