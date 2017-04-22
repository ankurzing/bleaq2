function [ulEliteFunctionValue, llEliteFunctionValue, ulEliteIndiv, llEliteIndiv, ulFunctionEvaluations, llFunctionEvaluations,llCalls,gen,ulDim,llDim,ulPopSize,llPopSize] = smd4(ulPopSize, llPopSize, ulMaxGens, llMaxGens, ulDim, llDim)

problemName = 'smd4';             % Test problem name

if nargin==0
    ulPopSize=20;                 % Size of UL population
    ulMaxGens=2000;               % Maximum number of generations allowed at UL
    ulDim=2;                      % Number of UL dimensions
    llPopSize=30;                 % Size of LL population
    llMaxGens=2000;               % Maximum number of generations allowed at LL
    llDim=3;                      % Number of LL dimensions
end

r = floor(ulDim/2);
p = ulDim - r;
q = llDim - r;

size_xu1 = p;
size_xu2 = r;
size_xl1 = q;
size_xl2 = r;

ulDimMin = [-5*ones(1,p) -1*ones(1,r)];       % Minimum bound accross UL dimensions
ulDimMax = [10*ones(1,p) 1*ones(1,r)];        % Maximum bound accross UL dimensions

eps = 0.00001;
llDimMin = [-5*ones(1,q) zeros(1,r)];         % Minimum bound accross LL dimensions
llDimMax = [10*ones(1,q) exp(1)*ones(1,r)];   % Maximum bound accross LL dimensions

ulStoppingCriteria = 1e-4;
llStoppingCriteria = 1e-5;

[ulEliteFunctionValue, llEliteFunctionValue, ulEliteIndiv, llEliteIndiv, ulFunctionEvaluations, llFunctionEvaluations,llCalls,gen]=ulSearch(problemName, ulPopSize, ulMaxGens, ulDim, ulDimMin, ulDimMax, llPopSize, llMaxGens, llDim, llDimMin, llDimMax, ulStoppingCriteria, llStoppingCriteria);

save('smd4');
