function msmd1()

problemName = 'msmd1';           % SMD1 test problem modified to multi-objective problem

ulPopSize=40;                    % The size of the population
ulMaxGens=500;                   % The maximum number of generations allowed in a run
ulDim=2;                         % Number of dimensions

llPopSize=40;                    % The size of the population
llMaxGens=2000;                  % The maximum number of generations allowed in a run
llDim=3;                         % Number of dimensions

r = floor(ulDim/2);
p = ulDim - r;
q = llDim - r;

size_xu1 = p;
size_xu2 = r;
size_xl1 = q;
size_xl2 = r;

ulDimMin = -5*ones(1,ulDim);       % Minimum value accross dimensions
ulDimMax = 10*ones(1,ulDim);       % Maximum value accross dimensions

llDimMin = [-5*ones(1,q+r)];       % Minimum value accross dimensions
llDimMax = [10*ones(1,q+r)];       % Maximum value accross dimensions

ulStoppingCriteria = 1e-4;
llStoppingCriteria = 1e-5;

[ulEliteFunctionValue, llEliteFunctionValue, ulEliteIndiv, llEliteIndiv, ulFunctionEvaluations, llFunctionEvaluations]=ulSearch(problemName, ulPopSize, ulMaxGens, ulDim, ulDimMin, ulDimMax, llPopSize, llMaxGens, llDim, llDimMin, llDimMax, ulStoppingCriteria, llStoppingCriteria)

save('msmd1');