function [StoppingCriteria, stoppingCondition] = terminationCheck(gen, tag, ulPop, llPop, ulDim, llDim, stoppingParameters)

    %If StoppingCriteria is returned 1 then the code terminates otherwise not
    StoppingCriteria = 0;
    stoppingCondition = [];
    
    [accuracyStoppingCriteria, accuracyStoppingCondition] = accuracyTerminationCheck(stoppingParameters);
    if accuracyStoppingCriteria == 1
        StoppingCriteria = 1;
        stoppingCondition = accuracyStoppingCondition;
    end
    
    [improvementStoppingCriteria, improvementStoppingCondition] = improvementTerminationCheck(gen, tag, ulPop, llPop, ulDim, llDim, stoppingParameters);
    if improvementStoppingCriteria==1
        StoppingCriteria = 1;
        stoppingCondition = improvementStoppingCondition;
    end
    [varianceStoppingCriteria, varianceStoppingCondition] = varianceTerminationCheck(gen, tag, ulPop, llPop, ulDim, llDim, stoppingParameters);
    if varianceStoppingCriteria==1
        StoppingCriteria = 1;
        stoppingCondition = varianceStoppingCondition;
    end
    
function [StoppingCriteria, stoppingCondition] = accuracyTerminationCheck(stoppingParameters)
    
    %desiredAccuracy = stoppingParameters.ulEpsilonStopping;
    desiredAccuracy = 1e-2;
    testProblemName = stoppingParameters.testProblemName;
    eliteFunctionValue = stoppingParameters.eliteFunctionValueAtGen(end);
    llEliteFunctionValue = stoppingParameters.llEliteFunctionValue;
    
    StoppingCriteria = 0;
    stoppingCondition = [];
    
    if strcmp(testProblemName(1:2),'tp')
        ulBestKnownFunctionValue(1) = 225.0000;
        llBestKnownFunctionValue(1) = 100.0000;

        ulBestKnownFunctionValue(2) = 0.0000;
        llBestKnownFunctionValue(2) = 100.0000;

        ulBestKnownFunctionValue(3) = -18.6787;
        llBestKnownFunctionValue(3) = -1.0156;

        ulBestKnownFunctionValue(4) = -29.2;
        llBestKnownFunctionValue(4) = 3.2;

        ulBestKnownFunctionValue(5) = -3.6;
        llBestKnownFunctionValue(5) = -2.0;

        ulBestKnownFunctionValue(6) = -1.2091;
        llBestKnownFunctionValue(6) = 7.6145;

        ulBestKnownFunctionValue(7) = -1.96;
        llBestKnownFunctionValue(7) = 1.96;

        ulBestKnownFunctionValue(8) = 0;
        llBestKnownFunctionValue(8) = 100;
        
        ulBestKnownFunctionValue(9) = 0;
        llBestKnownFunctionValue(9) = 1.0;

        ulBestKnownFunctionValue(10) = 0;
        llBestKnownFunctionValue(10) = 1.0;

        i = str2num(testProblemName(3:end));

        if abs(-ulBestKnownFunctionValue(i)-eliteFunctionValue)<desiredAccuracy %&& abs(-llBestKnownFunctionValue(i)-llEliteFunctionValue)<0.01
            StoppingCriteria = 1;
            stoppingCondition = 'Accuracy based'
        end

    elseif strcmp(testProblemName(1:3),'smd')
        
        i = str2num(testProblemName(4:end));
        
        eliteIndiv = stoppingParameters.eliteIndiv;
        llEliteIndiv = stoppingParameters.llEliteIndiv;
        
        [ulBestKnownMember llBestKnownMember ulBestKnownFunctionValue llBestKnownFunctionValue]=getOptimalSolutionSMD(length(eliteIndiv),length(llEliteIndiv),testProblemName);

        if abs(ulBestKnownFunctionValue-eliteFunctionValue) < desiredAccuracy 
            if abs(llBestKnownFunctionValue-llEliteFunctionValue) < desiredAccuracy
                StoppingCriteria = 1;
                stoppingCondition = 'Accuracy based';
            end
        end

    end
   
function [StoppingCriteria, stoppingCondition] = varianceTerminationCheck(gen, tag, ulPop, llPop, ulDim, llDim, stoppingParameters)

    ulEpsilonStopping = stoppingParameters.ulEpsilonStopping;
    llEpsilonStopping = stoppingParameters.llEpsilonStopping;
    alphaStoppingInitial = stoppingParameters.alphaStoppingInitial;
    eliteConstrViolation = stoppingParameters.eliteConstrViolation;
    
    StoppingCriteria = 0;
    stoppingCondition = [];
    
    if sum(tag.ulPop==1) <= 1
        alphaStopping = Inf;
    else
        alphaStopping = sum(var([ulPop(tag.ulPop==1,:) llPop(tag.ulPop==1,:)]))/(ulDim+llDim);
        alphaStopping = alphaStopping/alphaStoppingInitial;
    end
    if alphaStopping<ulEpsilonStopping && eliteConstrViolation<=0
        StoppingCriteria = 1;
        stoppingCondition = 'Variance based';
    end
    
    
function [StoppingCriteria, stoppingCondition] = improvementTerminationCheck(gen, tag, ulPop, llPop, ulDim, llDim, stoppingParameters)

    ulEpsilonStopping = stoppingParameters.ulEpsilonStopping;
    llEpsilonStopping = stoppingParameters.llEpsilonStopping;
    eliteFunctionValueAtGen = stoppingParameters.eliteFunctionValueAtGen;
    eliteConstrViolation = stoppingParameters.eliteConstrViolation;
    
    StoppingCriteria = 0;
    stoppingCondition = [];
    
    %Stops if normalized improvement is less than ulEpsilonStopping in ulImprovementGenDiff generations
    ulEpsilonStoppingImprovement = ulEpsilonStopping;
    ulImprovementGenDiff = 200;

    if (gen>ulImprovementGenDiff)
        if abs(eliteFunctionValueAtGen(gen)-eliteFunctionValueAtGen(gen-ulImprovementGenDiff)) == 0
            betaStopping = 0;
        elseif abs(eliteFunctionValueAtGen(gen)+eliteFunctionValueAtGen(1)) == 0
            betaStopping = abs(eliteFunctionValueAtGen(gen)-eliteFunctionValueAtGen(gen-ulImprovementGenDiff));
        else
            betaStopping = abs(eliteFunctionValueAtGen(gen)-eliteFunctionValueAtGen(gen-ulImprovementGenDiff))/abs(eliteFunctionValueAtGen(gen)+eliteFunctionValueAtGen(1));
        end
    else
        betaStopping = Inf;
    end
    if betaStopping<ulEpsilonStoppingImprovement && eliteConstrViolation<=0
        StoppingCriteria = 1;
        stoppingCondition = 'Improvement based';
    end