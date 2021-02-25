function [ ret, implemented_str ] = getTemporalNetworkMeasure( dyn_net, measure )
%GETTEMPORALNETWORKMEASURE Compute a network measure over a sequence of adjacency matrices
%   This is a meta/calling function for usability. Given a sequence of adjacency matrices [NxNxT], compute the network measure/feature at each time step.

%   This method is optimized for speed, so the time 'for' loop is on the
%   interior and replicated several times. This saves on a strcmp each
%   iteration. A 'helper' output is written to show currently implemented methods.

%   @input dyn_net, a [NxNxT] sequence of adjacency matrices over T time steps.
%   @input measure, a string specifying the network measure to compute
%   @output ret, a T-length vector of network measures
%   @output implemented_str [optional], a cell array of implemented method strings, for ease of method selection (e.g. to loop over all features)

%   IB, Update: Testing on R2013a, parfor now only runs a small constant
%   slower serially vs. regular 'for.' removed regular for 'if' blocks

%   Last updated: fixed return from clique number. 

%   IB, last updated: 3/24/14

%% keep track of implemented methods for helper output
implemented_str = sort({'max component', 'num components', 'average degree', 'average path', 'clustering coeff', 'density', 'radius', 'diameter', 'transitivity', 'clique number'});


if(~exist('dyn_net', 'var') || isempty(dyn_net))
    ret = [];
else
    
    [~,~,s] = size(dyn_net);
    
    %% preallocate
    ret = zeros(s, 1);
        
    %% logic for measures
    if(exist('measure', 'var')) %if measure given
        if(strcmpi(measure, 'max component'))
            parfor(i=1:s, matlabpool('size'))
                [~, ret(i)] = connectedComponentMeasures( dyn_net(:,:, i) );
            end
        elseif(strcmpi(measure, 'num components'))
            parfor(i=1:s, matlabpool('size'))
                [ret(i), ~] = connectedComponentMeasures(dyn_net(:,:, i));
            end
        elseif(strcmpi(measure, 'clique number'))
            parfor(i=1:s, matlabpool('size'))
                ret(i) = cliqueNumber(dyn_net(:,:, i));
            end
        elseif(strcmpi(measure, 'average degree'))
            parfor(i=1:s, matlabpool('size'))
                ret(i) = averageDegree(dyn_net(:,:, i));
            end
        elseif(strcmpi(measure, 'average path'))
            parfor(i=1:s, matlabpool('size'))
                [ret(i), ~] = averagePathLength(dyn_net(:,:, i));
            end
        elseif(strcmpi(measure, 'clustering coeff'))
            parfor(i=1:s, matlabpool('size'))
                [~, ret(i), ~] = clustCoeff(dyn_net(:,:,i));
            end
        elseif(strcmpi(measure, 'transitivity'))
            parfor(i=1:s, matlabpool('size'))
                [ret(i), ~, ~] = clustCoeff(dyn_net(:,:,i));
            end
        elseif(strcmpi(measure, 'density'))
            parfor(i=1:s, matlabpool('size'))
                ret(i) = edgeDensity(dyn_net(:,:,i));
            end
        elseif(strcmpi(measure, 'radius'))
            parfor(i=1:s, matlabpool('size'))
                [ret(i), ~] = radiusAndDiameter(dyn_net(:,:,i));
            end
        elseif(strcmpi(measure, 'diameter'))
            parfor(i=1:s, matlabpool('size'))
                [~, ret(i)] = radiusAndDiameter(dyn_net(:,:,i));
            end
        else
            fprintf(['String ' measure ' not valid, implemented methods:\n' strjoin(implemented_str, '\n') '\n']);
        end
    else %if measure not given
        fprintf(['No network measure given, implemented methods:\n' strjoin(implemented_str, '\n') '\n']);
    end
end

end

