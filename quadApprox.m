function approxmodel = quadApprox(y, X)

% modification - 2016/08/19
% eliminates input "model"
% "linear", "pure-quadratic", "quadratic" model is chosen based on bic
% number

%enter y as a column vector 
%enter X as a matrix where each column represents a variable
%model is 'quadratic' or 'purequadratic' or 'interaction' or 'linear'
%general model is y = f(x+deltax) = f(x) + grad(x)*deltax +
%1/2*deltax'*sqmatrix*deltax
%Approximation y = f(x) = constant + linear'*x + x'*sqmatrix*x

%Note:
%In case of linear model sqmatrix is zero
%In case of purequadratic model off diagonal terms of sqmatrix are zero
%In case of interaction model the diagonal terms of sqmatrix are zero

%Important:
%Backslash operator is used for regression
%If the matrix in the later part of the code happens to be rank deficient
%then the backslash operator issues a warning and produces the basic least-
%squares solution (coefficients as 0) out of the infinitely many solutions.

dimensions = size(X,2);
datasetSize = size(X,1);
if (datasetSize<dimensions)
    error('Cannot compute model as the datasetSize is smaller than dimensions.');
end

try
    if sum(isnan(y))>=1 || sum(sum(isnan(X)))>=1
        error('Cannot proceed as data contains NaN elements in it');
    end
catch
    display('hi');
end

modelConsidered = {'linear','purequadratic','quadratic'};
numModel = length(modelConsidered);
XX = cell(1,numModel); output = cell(1,numModel);
mse = zeros(1,numModel); bic = zeros(1,numModel);

for i = 1:numModel
    model = modelConsidered{i};
    XX{i} = x2fx(X, model);
    output{i}.beta = XX{i}\y;
    mse(i) = real(sum((y-XX{i}*output{i}.beta).^2)/length(y));
    bic(i) = size(XX{i},2)*log(size(XX{i},1))/size(XX{i},1) + 2*log(mse(i));      
end

[~,index] = min(bic);

modelSelected = modelConsidered(index);
constant = output{index}.beta(1);
linear = output{index}.beta(2:2+dimensions-1);
sqmatrix = zeros(dimensions,dimensions);

if strcmp(modelSelected,'purequadratic')
    diagonal = output{i}.beta(end-dimensions+1:end);
    for i=1:dimensions %#ok<ALIGN>
        sqmatrix(i,i) = diagonal(i);
	end
end

if strcmp(modelSelected,'quadratic')
    cross = output{i}.beta(2+dimensions:end-dimensions);
	diagonal = output{i}.beta(end-dimensions+1:end);
	k=0;
	for i=1:dimensions
        sqmatrix(i,i) = diagonal(i);
        for j=i+1:dimensions
            k=k+1;
            sqmatrix(i,j) = cross(k)/2;
            sqmatrix(j,i) = cross(k)/2;
        end
	end
end

stdXX = std(XX{index});
stdy = std(y);
stdXX(stdXX==0) = -realmin;
stdy(stdy==0) = -realmin;
XXNorm = (XX{index}-ones(size(XX{index},1),1)*mean(XX{index}))./(ones(size(XX{index},1),1)*stdXX);
yNorm = (y-mean(y))./stdy;
betaNorm = XXNorm\yNorm;
mseNorm = sum((yNorm-XXNorm*betaNorm).^2)/length(yNorm);

approxmodel = struct('model',modelSelected,'constant',real(constant),...
                'linear',real(linear),'sqmatrix',real(sqmatrix),...
                'mse',mse(index),'mseNorm',mseNorm,'bic',bic(index));

if isnan(approxmodel.mse)
    error('The code has ended up with NAN values. One of the reasons could be duplicate rows in the input matrix that makes the system under-defined and a solution cannot be found using least-squares.')
end