function [xvalR2, xvalsqErr, yhatLOOCV, coefStruct , lsFitStruct , permStruct ] = ... 
    nc_FitAndEvaluateModels(y, x, model, crossvalidate, bootIter, params , permIter)
% Calculate R2 using leave one out cross validation for a variety of models
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% J. Faskowitz edit;
% using the script initially provided here: 
% https://github.com/jyeatman/lifespan/blob/master/nc_FitAndEvaluateModels.m
% editing it for the project
%
% By the time data reaches this function, it should be conditioned so that
% it does not have NaNs or so that it only looks at a specific age subset
% that one would want 
% 
% adding permutation functionality
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% [xvalR2, xvalsqErr, yhatLOOCV, coefStruct , lsFitStruct , permStruct] = 
%    nc_FitAndEvaluateModels(y, x, model, crossvalidate, bootIter, params , permIter)
%
% This function will fit and evaluate a number of different types of models 
% using cross validation
%
% The following model classes have been implimented:
% 'linear'          - Linear model
% 'quadratic'       - Second order polynomial
% 'piecewise'       - Piecewise linear model consisting of a line with a
%                     slope joined to a flat line by a hinge
% 'piecewise2'      - Piecewise linear model consisting of 2 lines with
%                     independent slopes each connected by a flat line in 
%                     the middle with 2 independent hinges
% 'piecewisenoflat' - Same as piecewise but the second line has an
%                     independent slope (rather than being flat)
% 'piecewise2noflat'- Same as piecewise2 but the middle line also has a
%                     slope
% 'exponent'        - y = a*x^n + c
% 'lowess'          - local regression model. Requires kendrick kays code.
%                     See https://github.com/knk
% 'poisson'         - A Poisson curve
%
% Inputs:
%
% y             - vector of y values
% x             - vector of x values
% model         - string denoting model type (e.g., 'linear')
% crossvalidate - estimate R2 using cross validation (logical)
% bootIter      - number of bootstrap iterations (scalar >=1)
%
% Outputs:
%
% R2    - Cross validated R2 (model accuracy)
% sqErr - squared error
% yhat  - model prediction
% coef  - model coefficients
% permPval - permutation pval of xval R2 
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications.

%% Argument checking
if ~exist('crossvalidate','var') || isempty(crossvalidate)
    crossvalidate = true;
end
if ~exist('bootIter','var') || isempty(bootIter)
    bootstrap = false;
else
    bootstrap = true;
end

% If a lowess model was requested then parse the bandwidth parameter
if ~isempty(strfind(model,'lowess'))
    if length(model)>6
        params = str2double(model(7:end));
        model = 'lowess';
    end
end

% %% Echo to disp to make sure things look good
% disp(strcat('xval?: ', num2str(crossvalidate == 1)));
% disp(strcat('bootIter?:', num2str(bootIter)));
% disp(strcat('permIter?:', num2str(permIter)));
% good to go

%% Setup vars

% Set up structure for coeficients
coefStruct = struct('full', [], ...
    'xval', [], ...
    'boot', [], ...
    'name',[],...
    'x',[],...
    'y',[]);

lsFitStruct = struct( 'regressStat', [], ...
    'R2', [],...
    'SSE', [] );

permStruct = struct ( 'permDist', [], ...
    'permMat', [], ...
    'permPvalR2', [] );

% save points in this struct
coefStruct.x = x; 
coefStruct.y = y;
% And other outputs
xvalR2 = ''; 
xvalsqErr = ''; 
yhatLOOCV = '';

if crossvalidate == 1
    
    % preallocate so matlab does not complain
    
    % indLOOCV, rows = n -1 , cols = n
    indLOOCV = zeros( (length(x)-1) , length(x)) ;
    leftoutLOOCV = zeros( 1 , length(x) ) ;
    yhatLOOCV = zeros( 1 , length(x) ) ;
    
    % First generate cross validation indices
    % for leave one out x-val
    for idx = 1:length(x)
        indLOOCV(:,idx) = horzcat(1:(idx-1), (idx+1):length(x))';
        leftoutLOOCV(idx) = idx;
    end
    
end

% if xval and permIter not empty..;
if exist('permIter','var') && ~isempty(permIter)
       
    % check to see if we have a good number
    if ~isnumeric(permIter)
        permIter = 5000 ;
    end
  
else
    permIter = 0 ;
end

%% Fit model

% case statement for different model fits
switch(model)
    
    case {'linear' 'lin'} 
    %% linear
        
        coefStruct.name = model;
        % Make regressor matrix by concatenating x with a column of ones
        X = horzcat(x, ones(length(x),1));
        
        % Fit the full model without leaving out data
        % FULL FIT
        [coefStruct.full, ~,~,~, lsFitStruct.regressStat] = regress(y, X);     

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % make lsFit stats
        yhatFull = coefStruct.full(1) .* X(:,1) + coefStruct.full(2) ;
        residr = norm(y - yhatFull) ;
        SSE = residr.^2; % Error sum of squares.
        TSS = norm(y-mean(y))^2;     % Total sum of squares.
        R2 = 1 - SSE/TSS; 
        
        lsFitStruct.R2 = R2 ;
        lsFitStruct.SSE = SSE ;
        lsFitStruct.TSS = TSS ;
        lsFitStruct.RMSE = sqrt(SSE/(size(x,1))) ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        % Bootstrap
        if bootstrap == 1
            coefStruct.boot = bootstrp(bootIter, @regress, y, X);
        end
        
        % Cross validate
        if crossvalidate == 1
            for idx = 1:size(indLOOCV,2)
                % Fit the model to the data leaving one out
                coefStruct.xval(:,idx) = regress(y(indLOOCV(:,idx)), X(indLOOCV(:,idx),:));
                % Predict the y value for the one that was left out
                yhatLOOCV(idx)   = coefStruct.xval(1,idx) .* X(leftoutLOOCV(idx),1) + coefStruct.xval(2,idx);
            end
        end % xval
        
        % permutation test 
        if permIter > 0

            % preallocate array 
            permMat = zeros(size(coefStruct.full,1),permIter);
            permDist = zeros( [ permIter 1 ] ) ;
            indMat = zeros( [ size(x,1) permIter ] ) ;             

            % get randperm indicies
            for idx=1:permIter
                indMat(:,idx) = randperm(size(x,1))' ;
            end

            % loop for perms
            for idx=1:permIter
                
                % lets get an R2 at each permuation 
                
                % run regression with the permuted X mat and full y
                permMat(:,idx) = regress(y, X(indMat(:,idx),:)) ;
                
                % get the predicted Y
                % from linear equation
                % we should use unshuffled X values here?! NO
                yhatPerm = permMat(1,idx) .* X(indMat(:,idx),1) + permMat(2,idx);

                % residuals 
                residr = norm(y - yhatPerm) ;
                SSE = residr.^2; % Error sum of squares.
                TSS = norm(y-mean(y))^2;     % Total sum of squares.
                % R2 at this permutation
                permDist(idx) = 1 - SSE/TSS; 
                
            end    

            % record pvalR2 from permuation test
            % add 1 to numerator & demoninatior: https://www.ncbi.nlm.nih.gov/pubmed/21044043
            permStruct.permPvalR2 = (sum(permDist > lsFitStruct.R2) + 1) / (permIter + 1) ;
            permStruct.permDist = permDist ; 
            % record for future interests
            permStruct.permMat = permMat ;
            
        end % perm 
        
    case {'quadratic' 'quad'}
    %% quadratic
        
        coefStruct.name = model;
        
        % Make regressor matrix by concatenating x^2, x and a constant
        X = horzcat(x.^2, x, ones(length(x),1));
        
        % Fit the full model without leaving out data
        % FULL FIT
        [coefStruct.full, ~,~,~, lsFitStruct.regressStat] = regress(y, X);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % make lsFit stats
        yhatFull = coefStruct.full(1) .* X(:,1) ...
            + coefStruct.full(2) .* X(:,2) + coefStruct.full(3);
        residr = norm(y - yhatFull) ;
        SSE = residr.^2; % Error sum of squares.
        TSS = norm(y-mean(y))^2;     % Total sum of squares.
        R2 = 1 - SSE/TSS; 
            
        lsFitStruct.R2 = R2 ;
        lsFitStruct.SSE = SSE ;
        lsFitStruct.TSS = TSS ;
        lsFitStruct.RMSE = sqrt(SSE/(size(x,1))) ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Bootstrap
        if bootstrap == 1
            coefStruct.boot = bootstrp(bootIter, @regress, y, X);
        end
        
        % Cross validate
        if crossvalidate == 1
            for idx = 1:size(indLOOCV,2)
                % Fit the model to the data leaving one out
                coefStruct.xval(:,idx) = regress(y(indLOOCV(:,idx)), X(indLOOCV(:,idx),:));
                % Predict the y value for the one that was left out
                yhatLOOCV(idx)   = coefStruct.xval(1,idx) .* X(leftoutLOOCV(idx),1) ...
                    + coefStruct.xval(2,idx) .* X(leftoutLOOCV(idx),2) + coefStruct.xval(3,idx);
            end
        end
        
        % permutation test 
        if permIter > 0

            % preallocate array 
            permMat = zeros(size(coefStruct.full,1),permIter);
            permDist = zeros( [ permIter 1 ] ) ;
            indMat = zeros( [ size(x,1) permIter ] ) ;             

            % get randperm indicies
            for idx=1:permIter
                indMat(:,idx) = randperm(size(x,1))' ;
            end

            % loop for perms
            for idx=1:permIter
                
                % lets get an R2 at each permuation 
                
                % run regression with the permuted X mat and full y
                permMat(:,idx) = regress(y, X(indMat(:,idx),:)) ;
                
                % get the predicted Y
                % from equation
                yhatPerm = permMat(1,idx) .* X(indMat(:,idx),1) ...
                    + permMat(2,idx) .* X(indMat(:,idx),2) ...
                    + permMat(3,idx);

                % residuals 
                residr = norm(y - yhatPerm) ;
                SSE = residr.^2; % Error sum of squares.
                TSS = norm(y-mean(y))^2;     % Total sum of squares.
                permDist(idx) = 1 - SSE/TSS; 
                
            end 
            
            % record pvalR2 from permuation test
            permStruct.permPvalR2 = (sum(permDist > lsFitStruct.R2) + 1) / (permIter + 1) ;
            permStruct.permDist = permDist ; 
            permStruct.permMat = permMat ;
            
        end % perm 

%% piecewise
    case {'piecewise'}
        coefStruct.name = model;
        % We will try cutpoints between age 12 and 30
        c = 12:2:30;
        
        % Fit the full model without leaving out data
        coefStruct.full = piecewiseFit(x, y, c);
        
        % Bootstrap
        if bootstrap == 1
            coefStruct.boot = bootstrp(bootIter,@(x,y) piecewiseFit(x,y,c),x,y);
        end
        
        if crossvalidate == 1
            for idx = 1:size(indLOOCV,2)
                % Fit the piecewise model on the subset of the data
                coefStruct.xval(:,idx) = piecewiseFit(x(indLOOCV(:,idx)), y(indLOOCV(:,idx)), c);
                % Evaluate the model at the left out point
                yhatLOOCV(idx) = piecewiseEval(coefStruct.xval(:,idx)', x(leftoutLOOCV(idx)));
            end
        end
    case {'piecewise2'}
        coefStruct.name = model;
        % We will try cutpoints between age 12 and 30 and 62 and 80
        %c = [14:1:30; 58:1:74]';
        c = [14:2:28; 60:2:74]';
        % Fit the full model without leaving out data
        coefStruct.full = piecewiseFit(x, y, c);
        
        % Bootstrap
        if bootstrap == 1
            coefStruct.boot = bootstrp(bootIter,@(x,y) piecewiseFit(x,y,c),x,y);
        end
        
        if crossvalidate == 1
            for idx = 1:size(indLOOCV,2)
                % Fit the piecewise model on the subset of the data
                coefStruct.xval(:,idx) = piecewiseFit(x(indLOOCV(:,idx)), y(indLOOCV(:,idx)), c);
                % Evaluate the model at the left out point
                yhatLOOCV(idx) = piecewiseEval(coefStruct.xval(:,idx)', x(leftoutLOOCV(idx)));
            end
        end
    case {'piecewisenoflat'}
        coefStruct.name = model;
        % We will try cutpoints between age 12 and 30
        c = 12:2:30;
        
        % Fit the full model without leaving out data
        coefStruct.full = piecewiseFit(x, y, c,'fit');
        
        % Bootstrap
        if bootstrap == 1
            coefStruct.boot = bootstrp(bootIter,@(x,y) piecewiseFit(x,y,c,'fit'),x,y);
        end
        
        if crossvalidate == 1
            for idx = 1:size(indLOOCV,2)
                % Fit the piecewise model on the subset of the data
                coefStruct.xval(:,idx) = piecewiseFit(x(indLOOCV(:,idx)), y(indLOOCV(:,idx)), c,'fit');
                % Evaluate the model at the left out point
                yhatLOOCV(idx) = piecewiseEval(coefStruct.xval(:,idx)', x(leftoutLOOCV(idx)));
            end
        end
    case {'piecewise2noflat'}
        coefStruct.name = model;
        % We will try cutpoints between age 12 and 30 and 62 and 80
        %c = [14:1:30; 58:1:74]';
        c = [14:2:28; 60:2:74]';
        % Fit the full model without leaving out data
        coefStruct.full = piecewiseFit(x, y, c,'fit');
        
        % Bootstrap
        if bootstrap == 1
            coefStruct.boot = bootstrp(bootIter,@(x,y) piecewiseFit(x,y,c,'fit'),x,y);
        end
        
        if crossvalidate == 1
            for idx = 1:size(indLOOCV,2)
                % Fit the piecewise model on the subset of the data
                coefStruct.xval(:,idx) = piecewiseFit(x(indLOOCV(:,idx)), y(indLOOCV(:,idx)), c,'fit');
                % Evaluate the model at the left out point
                yhatLOOCV(idx) = piecewiseEval(coefStruct.xval(:,idx)', x(leftoutLOOCV(idx)));
            end
        end
        
    case {'exponent' 'exponential' 'exp'}
    %% exp    
        
        coefStruct.name = model;
        
        % Set options for nonlinear fitting
        options = optimset('Display','off');
        % Write out the function for the exponential
        expfun = @(p,xval) p(1).*xval.^p(2) + p(3);
        
        % Do a least squares fit of a line to seed the nonlin fit
        L = polyfit(x, y, 1);
        % Fit the full model without leaving out data
        % FULL FIT
        coefStruct.full = lsqcurvefit(expfun, [L(1) 1 L(2)], x, y, [], [], options);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % make lsFit stats
        yhatFull = feval(expfun, coefStruct.full(:), x);
        residr = norm(y - yhatFull) ;
        SSE = residr.^2; % Error sum of squares.
        TSS = norm(y-mean(y))^2;     % Total sum of squares.
        R2 = 1 - SSE/TSS; 
            
        lsFitStruct.R2 = R2 ;
        lsFitStruct.SSE = SSE ;
        lsFitStruct.TSS = TSS ;
        lsFitStruct.RMSE = sqrt(SSE/(size(x,1))) ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Bootstrap
        if bootstrap == 1
            coefStruct.boot = bootstrp(bootIter,@(x,y) lsqcurvefit(expfun, [L(1) 1 L(2)], x, y, [], [], options), x, y);
        end
        
        % Cross validate
        if crossvalidate == 1
            for idx = 1:size(indLOOCV,2)
                % Do a least squares fit of a line to seed the nonlin fit
                L = polyfit(x(indLOOCV(:,idx)), y(indLOOCV(:,idx)), 1);
                
                % Fit the exponential seeding with a second order polynomial
                coefStruct.xval(:,idx) = lsqcurvefit(expfun, [L(1) 1 L(2)], x(indLOOCV(:,idx)), y(indLOOCV(:,idx)), [], [], options);
                % Evaluate the exponential to get the prediction for the
                % leftout point
                yhatLOOCV(idx) = feval(expfun, coefStruct.xval(:,idx), x(leftoutLOOCV(idx)));
            end
        end
        
        % permutation test 
        if permIter > 0

            % preallocate array 
            % size along the second dim here
            permMat = zeros(size(coefStruct.full,2),permIter);
            permDist = zeros( [ permIter 1 ] ) ;
            indMat = zeros( [ size(x,1) permIter ] ) ;             

            % get randperm indicies
            for idx=1:permIter
                indMat(:,idx) = randperm(size(x,1))' ;
            end

            % loop for perms
            for idx=1:permIter
                
                % lets get an R2 at each permuation 
                
                % run regression with the permuted X mat and full y
                L = polyfit(x(indMat(:,idx)), y, 1);
                permMat(:,idx) = lsqcurvefit(expfun, [L(1) 1 L(2)], x(indMat(:,idx)), y, [], [], options) ;
                
                % get the predicted Y
                % from equation
                yhatPerm = feval(expfun, permMat(:,idx), x(indMat(:,idx)));

                % residuals 
                residr = norm(y - yhatPerm) ;
                SSE = residr.^2; % Error sum of squares.
                TSS = norm(y-mean(y))^2;     % Total sum of squares.
                permDist(idx) = 1 - SSE/TSS; 
                
            end 
            
            % record pvalR2 from permuation test
            permStruct.permPvalR2 = (sum(permDist > lsFitStruct.R2) + 1) / (permIter + 1) ;
            permStruct.permDist = permDist ; 
            permStruct.permMat = permMat ;
            
        end % perm 
        
    case {'lowess' 'loes' 'loess'}
    %% local    
        
        coefStruct.name = model;
        
        % Transpose x and y into row vectors if necesary
        if size(x,1) > size(x,2)
            x = x';
        end
        if size(y,1) > size(y,2)
            y = y';
        end
        % kernal width
        if ~exist('params','var') || isempty(params)
            k = 15;
        else
            k = params;
        end
        
        % For the coeficients we save first the point at which the model is
        % evaluated and second the predicted value at that point.
        coefStruct.full(:,1) = min(x):0.1:max(x);
        % FULL FIT        
        coefStruct.full(:,2) = localregression(x, y, coefStruct.full(:,1) , 3, [], k);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % get the yhat for local regression

        % first get xVals and indicies of intersection
        % this will only work if x0 and x vars intersect
        [xFound,ind] = intersect(coefStruct.full(:,1),x);
        yhatFull = zeros([size(x,1) 1]);
        for idx=1:size(x,1)
            yhatFull(idx) = coefStruct.full(ind(x(idx) == xFound),2);           
        end
            
        % make lsFit stats
        residr = norm(y - yhatFull) ;
        SSE = residr.^2; % Error sum of squares.
        TSS = norm(y-mean(y))^2;     % Total sum of squares.
        R2 = 1 - SSE/TSS; 
            
        % this is kind of dumb...but I'll calculate anyways
        lsFitStruct.R2 = R2 ;
        lsFitStruct.SSE = SSE ;
        lsFitStruct.TSS = TSS ;
        lsFitStruct.RMSE = sqrt(SSE/(size(x,1))) ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Bootstrap
        if bootstrap == 1
            coefStruct.boot = bootstrp(bootIter,@(x,y) localregression(x, y, coefStruct.full(:,1) , 3, [], k), x, y);
        end
        
        % Crossvalidate
        if crossvalidate == 1
            % Predict y from local regression
            for idx = 1:size(indLOOCV,2)
                yhatLOOCV(idx) = localregression(x(indLOOCV(:,idx)), y(indLOOCV(:,idx)), x(leftoutLOOCV(idx)) ,[],[],k);
            end
        end
        
        % not sure if good idea to permute this...
        % r-squared pseudo too...
        
    case{'poisson'}
    %% poisson    
        
        coefStruct.name = model;
        
        % FULL FIT
        coefStruct.full = fitPoissonCurve(x,y);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % make lsFit stats
        yhatFull = evalPoissonCurve(coefStruct.full(:), x);
        residr = norm(y - yhatFull) ;
        SSE = residr.^2; % Error sum of squares.
        TSS = norm(y-mean(y))^2;     % Total sum of squares.
        R2 = 1 - SSE/TSS; 
            
        lsFitStruct.R2 = R2 ;
        lsFitStruct.SSE = SSE ;
        lsFitStruct.TSS = TSS ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Bootstrap
        if bootstrap == 1
            coefStruct.boot = bootstrp(bootIter, @fitPoissonCurve, x, y);
        end
        
        % Cross validate
        if crossvalidate == 1
            for idx = 1:size(indLOOCV,2)
                % Fit the model to the data leaving one out
                coefStruct.xval(:,idx) = fitPoissonCurve(x(indLOOCV(:,idx)), y(indLOOCV(:,idx),:));
                % Predict the y value for the one that was left out
                yhatLOOCV(idx) = evalPoissonCurve(coefStruct.xval(:,idx), x(leftoutLOOCV(idx)));
            end
        end
           
        % permutation test 
        if permIter > 0

            % preallocate array 
            % size along the second dim here
            permMat = zeros(size(coefStruct.full,2),permIter);
            permDist = zeros( [ permIter 1 ] ) ;
            indMat = zeros( [ size(x,1) permIter ] ) ;             

            % get randperm indicies
            for idx=1:permIter
                indMat(:,idx) = randperm(size(x,1))' ;
            end

            % loop for perms
            for idx=1:permIter
                
                % lets get an R2 at each permuation  
                permMat(:,idx) = fitPoissonCurve(x(indMat(:,idx)), y);
                
                % get the predicted Y
                % from equation
                yhatPerm = evalPoissonCurve(permMat(:,idx), x(indMat(:,idx)));

                % residuals 
                residr = norm(y - yhatPerm) ;
                SSE = residr.^2; % Error sum of squares.
                TSS = norm(y-mean(y))^2;     % Total sum of squares.
                permDist(idx) = 1 - SSE/TSS; 
                
            end 
            
            % record pvalR2 from permuation test
            permStruct.permPvalR2 = (sum(permDist > lsFitStruct.R2) + 1) / (permIter + 1) ;
            permStruct.permDist = permDist ; 
            permStruct.permMat = permMat ;
            
        end % perm   
        
    otherwise
        error('%s is not an implimented model class',model)
end

%% X-val
% calculate xval squarred error and coef of determination

if crossvalidate == 1
    % match the dimensions of y and yhat
    if size(yhatLOOCV,2) ~= size(y,2)
        yhatLOOCV = yhatLOOCV';
    end
    
    % Calculate the squared error for the model prediction
    xvalsqErr = (y - yhatLOOCV).^2;
    
    % Calculate crossvalidated R2
    % comparing(predicted y vals) & (actual y vals)
    % how much variance in actual y vals can be explained by predicted yhat
    xvalR2 = calccod(yhatLOOCV,y);
    
else
    xvalR2 = nan;
end

return