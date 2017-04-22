function [ulEliteFunctionValue, llEliteFunctionValue, ulEliteIndiv, llEliteIndiv, ulFunctionEvaluations, llFunctionEvaluations,llCalls,gen,ulDim,llDim,ulPopSize,llPopSize] = smd14(ulPopSize, llPopSize, ulMaxGens, llMaxGens, ulDim, llDim)

problemName = 'smd14';             % Test problem name

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
q = floor((llDim - r)/2 - eps);
s = ceil((llDim - r)/2 + eps);

size_xu1 = p;
size_xu2 = r;
size_xl1 = q+s;
size_xl2 = r;

ulDimMin = -5*ones(1,ulDim);      % Minimum bound accross dimensions
ulDimMax = 10*ones(1,ulDim);      % Maximum bound accross dimensions

llDimMin = -5*ones(1,llDim);      % Minimum bound accross dimensions
llDimMax = 10*ones(1,llDim);      % Maximum bound accross dimensions

ulStoppingCriteria = 1e-4;
llStoppingCriteria = 1e-5;

[ulEliteFunctionValue, llEliteFunctionValue, ulEliteIndiv, llEliteIndiv, ulFunctionEvaluations, llFunctionEvaluations,llCalls,gen]=ulSearch(problemName, ulPopSize, ulMaxGens, ulDim, ulDimMin, ulDimMax, llPopSize, llMaxGens, llDim, llDimMin, llDimMax, ulStoppingCriteria, llStoppingCriteria);

save('smd14');
