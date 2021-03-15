function [uniques,numUnique] = count_unique(x,option)
%COUNT_UNIQUE  Determines unique values, and counts occurrences
%   [uniques,numUnique] = count_unique(x)
%
%   This function determines unique values of an array, and also counts the
%   number of instances of those values.
%
%   This uses the MATLAB builtin function accumarray, and is faster than
%   MATLAB's unique function for intermediate to large sizes of arrays for integer values.  
%   Unlike 'unique' it cannot be used to determine if rows are unique or 
%   operate on cell arrays.
%
%   If float values are passed, it uses MATLAB's logic builtin unique function to
%   determine unique values, and then to count instances.
%
%   Descriptions of Input Variables:
%   x:  Input vector or matrix, N-D.  Must be a type acceptable to
%       accumarray, numeric, logical, char, scalar, or cell array of
%       strings.
%   option: Acceptable values currently only 'float'.  If 'float' is
%           specified, the input x vector will be treated as containing
%           decimal values, regardless of whether it is a float array type.
%
%   Descriptions of Output Variables:
%   uniques:    sorted unique values
%   numUnique:  number of instances of each unique value
%
%   Example(s):
%   >> [uniques] = count_unique(largeArray);
%   >> [uniques,numUnique] = count_unique(largeArray);
%
%   See also: unique, accumarray

% Author: Anthony Kendall
% Contact: anthony [dot] kendall [at] gmail [dot] com
% Created: 2009-03-17

testFloat = false;
if nargin == 2 && strcmpi(option,'float')
    testFloat = true;
end

nOut = nargout;
if testFloat
    if nOut < 2
        [uniques] = float_cell_unique(x,nOut);
    else
        [uniques,numUnique] = float_cell_unique(x,nOut);
    end
else
    try %this will fail if the array is float or cell
        if nOut < 2
            [uniques] = int_log_unique(x,nOut);
        else
            [uniques,numUnique] = int_log_unique(x,nOut);
        end
    catch %default to standard approach
        if nOut < 2
            [uniques] = float_cell_unique(x,nOut);
        else
            [uniques,numUnique] = float_cell_unique(x,nOut);
        end
    end
end

end

function [uniques,numUnique] = int_log_unique(x,nOut)
%First, determine the offset for negative values
minVal = min(x(:));

%Check to see if accumarray is appropriate for this function
maxIndex = max(x(:)) - minVal + 1;
if maxIndex / numel(x) > 1000
    error('Accumarray is inefficient for arrays when ind values are >> than the number of elements')
end

%Now, offset to get the index
index = x(:) - minVal + 1;

%Count the occurrences of each index value
numUnique = accumarray(index,1);

%Get the values which occur more than once
uniqueInd = (1:length(numUnique))';
uniques = uniqueInd(numUnique>0) + minVal - 1;

if nOut == 2
    %Trim the numUnique array
    numUnique = numUnique(numUnique>0);
end
end 

function [uniques,numUnique] = float_cell_unique(x,nOut)

if ~iscell(x)
    %First, sort the input vector
    x = sort(x(:));
    numelX = numel(x);
    
    %Check to see if the array type needs to be converted to double
    currClass = class(x);
    isdouble = strcmp(currClass,'double');
    
    if ~isdouble
        x = double(x);
    end
    
    %Check to see if there are any NaNs or Infs, sort returns these either at
    %the beginning or end of an array
    if isnan(x(1)) || isinf(x(1)) || isnan(x(numelX)) || isinf(x(numelX))
        %Check to see if the array contains nans or infs
        xnan = isnan(x);
        xinf = isinf(x);
        testRep = xnan | xinf;
        
        %Remove all of these from the array
        x = x(~testRep);
    end
    
    %Determine break locations of unique values
    uniqueLocs = [true;diff(x) ~= 0];
else
    isdouble = true; %just to avoid conversion on finish
    
    %Sort the rows of the cell array
    x = sort(x(:));
    
    %Determine unique location values
    uniqueLocs = [true;~strcmp(x(1:end-1),x(2:end)) ~= 0] ;
end

%Determine the unique values
uniques = x(uniqueLocs);

if ~isdouble
    x = feval(currClass,x);
end

%Count the number of duplicate values
if nOut == 2
    numUnique = diff([find(uniqueLocs);length(x)+1]);
end
end
