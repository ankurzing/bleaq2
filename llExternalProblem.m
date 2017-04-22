function [functionValue equalityConstrVals inequalityConstrVals] = llExternalProblem(xu, xl)

    %Upper level TP1 implemented

    x = xu;
    y = xl;

    functionValue =   20*y(1) + 40*y(2) + 50*y(3) - y(1)^2 - 2*y(2)^2 - 2*y(3)^2 - x(1)*y(2) - x(2)*y(3);
    functionValue = functionValue + 30*y(4) + 20*y(5) - 3*y(4)^2 - 2*y(5)^2 - x(5)*y(5);

    equalityConstrVals = [];
    inequalityConstrVals(1) = y(1) + y(2) - 15;
    inequalityConstrVals(2) = y(2) + y(3) - 16;
    inequalityConstrVals(3) = y(4) + y(5) - 10;



    