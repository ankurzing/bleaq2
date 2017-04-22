function [ulEliteFunctionValue, llEliteFunctionValue, ulEliteIndiv, llEliteIndiv, ulFunctionEvaluations, llFunctionEvaluations,llCalls,gen,ulDim,llDim,ulPopSize,llPopSize] = tp1()

problemName = 'tp1';             % Test problem name

ulPopSize=50;                    % Size of UL population
ulMaxGens=2000;                  % Maximum number of generations allowed at UL
ulDim=2;                         % Number of UL dimensions

llPopSize=50;                    % Size of LL population
llMaxGens=2000;                  % Maximum number of generations allowed at LL
llDim=2;                         % Number of LL dimensions

ulDimMin = [-30 -30];            % Minimum value accross UL dimensions
ulDimMax = [30 15];              % Maximum value accross UL dimensions

llDimMin = 0*ones(1,llDim);      % Minimum value accross LL dimensions
llDimMax = 10*ones(1,llDim);     % Maximum value accross LL dimensions

ulStoppingCriteria = 1e-4;
llStoppingCriteria = 1e-5;


[ulEliteFunctionValue, llEliteFunctionValue, ulEliteIndiv, llEliteIndiv, ulFunctionEvaluations, llFunctionEvaluations,llCalls,gen]=ulSearch(problemName, ulPopSize, ulMaxGens, ulDim, ulDimMin, ulDimMax, llPopSize, llMaxGens, llDim, llDimMin, llDimMax, ulStoppingCriteria, llStoppingCriteria);

save('tp1');
