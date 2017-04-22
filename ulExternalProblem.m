function [functionValue equalityConstrVals inequalityConstrVals] = ulExternalProblem(xu, xl)

    %Upper level TP1 implemented

    x = xu;
    y = xl;

    functionValue =   60*x(1) + 80*x(2) + 70*x(3) - 2*x(1)^2 - 3*x(2)^2 - 2*x(3)^2 - x(1)*y(1) - 2*x(2)*y(3);
    functionValue = functionValue + 50*x(4) + 40*x(5) - 2*x(4)^2 - 2*x(5)^2 - x(4)*y(4);

    equalityConstrVals = [];
    inequalityConstrVals(1) = x(1) + x(2) - 22;
    inequalityConstrVals(2) = x(2) + x(3) - 25;
    inequalityConstrVals(3) = x(4) + x(5) - 20;

