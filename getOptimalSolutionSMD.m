function [upperMember lowerMember upperFunctionValue lowerFunctionValue]=getOptimalSolutionSMD(upperDimensions,lowerDimensions,testProblemName)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function call here
    fhandle = str2func(testProblemName);
    [upperMember lowerMember] = fhandle(upperDimensions,lowerDimensions);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [lowerFunctionValue lowerEqualityConstrVals lowerInequalityConstrVals]=llTestProblem(lowerMember, testProblemName, upperMember);
    [upperFunctionValue upperEqualityConstrVals upperInequalityConstrVals]=ulTestProblem(upperMember, lowerMember, testProblemName);
    
function [upperMember lowerMember] = smd1(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = zeros(1,lowerDimensions);
    
function [upperMember lowerMember] = smd2(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = [zeros(1,q) ones(1,r)];
    
function [upperMember lowerMember] = smd3(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = zeros(1,lowerDimensions);
    
function [upperMember lowerMember] = smd4(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = zeros(1,lowerDimensions);
    
function [upperMember lowerMember] = smd5(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = [ones(1,q) zeros(1,r)];
    
function [upperMember lowerMember] = smd6(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = floor((lowerDimensions - r)/2 - eps);
    s = ceil((lowerDimensions - r)/2 + eps);

    upperMember = zeros(1,upperDimensions);
    lowerMember = zeros(1,lowerDimensions);

function [upperMember lowerMember] = smd7(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = [zeros(1,q) ones(1,r)];

function [upperMember lowerMember] = smd8(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = [ones(1,q) zeros(1,r)];

function [upperMember lowerMember] = smd9(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = zeros(1,lowerDimensions);

function [upperMember lowerMember] = smd10(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = 1/sqrt(p+r-1)*ones(1,upperDimensions);
    lowerMember = [1/sqrt(q-1)*ones(1,q) atan(1/sqrt(p+r-1)*ones(1,r))];

function [upperMember lowerMember] = smd11(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = [zeros(1,q) exp(-1/sqrt(r))*ones(1,r)];

function [upperMember lowerMember] = smd12(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = 1/sqrt(p+r-1)*ones(1,upperDimensions);
    lowerMember = [1/sqrt(q-1)*ones(1,q) atan(1/sqrt(p+r-1)-1/sqrt(r))*ones(1,r)];
    
function [upperMember lowerMember] = smd13(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = [ones(1,p) zeros(1,r)];
    lowerMember = [zeros(1,p) ones(1,r)];

function [upperMember lowerMember] = smd14(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = floor((lowerDimensions - r)/2 - eps);
    s = ceil((lowerDimensions - r)/2 + eps);

    upperMember = [ones(1,p) zeros(1,r)];
    lowerMember = zeros(1,lowerDimensions);