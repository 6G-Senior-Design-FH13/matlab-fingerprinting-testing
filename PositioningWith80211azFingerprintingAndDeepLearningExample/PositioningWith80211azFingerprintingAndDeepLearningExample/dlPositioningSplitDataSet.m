function [training, validation] = dlPositioningSplitDataSet(data, labels, valfrac)
%dlPositioningSplitDataSet Randomly sample to create training and validation data
%   dlPositioningSplitDataSet(DATA, LABELS, VALFRAC) segments TRAINING data
%   for a Convolutional Neural Network model to fit and unseen VALIDATION
%   data to evaluate model performance. The features, DATA, is a 4D double
%   array and the LABELS are a vector of either position coordinates or
%   categorical location labels. The portion of the data assigned to
%   VALIDATION is determined by value of VALFRAC (e.g. 0.2 - results in
%   80/20 % split). If categorical, number of categories in TRAINING and
%   VALIDATION must be equal.

%   Copyright 2020 The MathWorks, Inc. 

% Determine number of validation samples.
numVal = floor(valfrac*(size(data,4)));
numClass = 7; % based on locations areas set during map creation
 
% Get logical index of validation set within training data
validx = randsampleidx(size(data,4),numVal);

% Filter training data and validation data
training.X = data(:,:,:,~validx);
validation.X = data(:,:,:,validx);

% Filter positioning training and validation labels
training.Y.regression = labels.position(:,~validx);
validation.Y.regression = labels.position(:,validx);

% Filter localization training and validation labels
training.Y.classification = labels.class(:,~validx);
validation.Y.classification = labels.class(:,validx);

% Check all classes are used in the validation and data set
if (numel(unique(training.Y.classification)) ~= numClass) || (numel(unique(validation.Y.classification)) ~= numClass)
    error("The resultant validation and training sets do not contain at least one value per class (expected 7 classes). Try increasing the number of features provided.")
end

end

function x = randsampleidx(n,k)
    % Select random samples until we have the required number of unique ones
    x = false(1,n);
    sumx = 0;
    while sumx < k
        x(randi(n,1,k-sumx)) = true; % sample w/replacement
        sumx = nnz(x); % count how many unique elements so far
    end
end


