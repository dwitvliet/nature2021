function G = fcn_group_average(A,dist,hemiid)
%FCN_GROUP_AVERAGE      build group-average matrix
%
%   G = FCN_GROUP_AVERAGE(A,DIST) takes as input a WEIGHTED "stack" of
%   connectivity matrices, A, with dimensions [N x N x SUBJECT] where N is
%   number of nodes and SUBJECT is the number of matrices in the stack. The
%   matrices MUST be weighted, and ideally with continuous weights (e.g.
%   fractional anisotropy rather than streamline count). The second input
%   is DIST, which is a pairwise distance matrix (i.e. DIST(i,j) is the
%   Euclidean distance between nodes i and j). The final input is an [N x
%   1] vector, HEMIID, which labels nodes as in the left (1) or right (2)
%   hemisphere.
%
%   This function estimates the average edge length distribution and builds
%   a group-averaged connectivity matrix that approximates this
%   distribution with density equal to the mean density across subjects.
%
%   The algorithm works as follows:
%   1. Estimate the cumulative edge length distrbiution.
%   2. Divide the distribution into M length bins, one for each edge that
%      will be added to the group-average matrix.
%   3. Within each bin, select the most edge that is most consistently
%      expressed across subjects, breaking ties according to average edge
%      weight (which is why A must be weighted).
%
%   The algorithm works separately on within/between hemisphere links.
%
%   Inputs:
%               A,      stack of connectivity matrices
%               DIST,   pairwise distance matrix
%
%   Outputs:    G,      binary group-average matrix
%
%   Richard Betzel, Indiana University, 2015
%   Updated by Andrea A-K, 2016

[n,~,nsub] = size(A);
C = sum(A > 0,3);
W = sum(A,3)./C;
Grp = zeros(n,n,2);
for j = 1:2
    if j == 1
        d = +(hemiid == 1)*(hemiid' == 2);
        d = d | d';
    else
        d = +(hemiid == 1)*(hemiid' == 1) | +(hemiid == 2)*(hemiid' == 2);
        d = d | d';
    end
    D = nonzeros(bsxfun(@times,(A > 0),dist.*triu(d)));
    M = round(length(D))/nsub;
    dist_hemi = dist.*d;
    [x,y] = ecdf(D);
    x = round(x.*M);
    G = zeros(n);
    for i = 1:M
        ind = (x >= (i - 1)) & (x < i);
        yy = y(ind);
        if ~isempty(yy)
            mask = dist_hemi >= min(yy) & dist_hemi <= max(yy);
            [u,v] = find(triu(mask,1));
            indx = (v - 1)*n + u;
            c = C(indx);
            w = W(indx);
            zz = sum(c == max(c));
            if zz == 1
                [~,indmax] = max(c);
                G(indx(indmax)) = 1;
            else
                aa = find(c == max(c));
                ww = w(aa);
                [~,indsort] = sort(ww,'descend');
                G(indx(aa(indsort(1)))) = 1;
            end
        end
    end
    Grp(:,:,j) = G;
end
G = sum(Grp,3); G = G + G';