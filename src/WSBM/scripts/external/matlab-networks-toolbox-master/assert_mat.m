function [ ] = assert_mat( A, B )
%ASSERT_MAT Overloading of assert statement for arbitrary-dimensioned data 
%   This method does not throw an error, it only outputs an error to the
%   screen. It is intended for handling 'assert' in the testAllFunctions
%   method not supported on MATLAB.

%   @input A, an arbitrarily-dimensioned data 
%   @input B, (optional) another arbitarily-dimensioned data, the same size as A.
%   @output NONE (note: usually there would be a boolean response, but the
%   testAllFunctions script doesn't ';' the output so this is cleaner). 


if(~exist('B', 'var'))
    ret = all(A(:));
else
    if(iscell(A))
        ret = all(cellfun(@all, cellfun(@eq, A, B, 'UniformOutput', 0)));
    else
        ret = all(A(:) == B(:));
    end
    
    
end

if(~ret)
    fprintf('Error\n');
end

end

