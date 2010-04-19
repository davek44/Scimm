function [minSegs,minNcut] = split_optimal(W,V)
%
%   minSegs = split_optimal(W,V)
%
%   Purpose
%   =======
%   Compute the normalized cut with the smallest normalized cut value.

%
%    if( V(1)<0 )
%        V=-V;
%    end

    % I did this, they used as input variable
    splits = 10;

    % here they set the interval size
    n=length(V);
    sortedV=sort(V);
    cut_int = (sortedV(n)-sortedV(1))/splits;
   
    minNcut=inf;
    splitV=inf;
    minSegs=sparse(n,1);
    cut=sortedV(1);

    for i=1:splits-1
        cut=cut+cut_int;

        % here they partition with cut
        segs=sparse(n,1);
        g=find(cut<V);
	segs(g)=1;

        % here they compute the cut goodness
        ncut=norm_cut(W,0,segs);
        
        % here they compare to past good cuts
        g1=find(cut>=V);
        if( ncut < minNcut )
            minNcut=ncut;
            minSegs=segs;
            splitV=cut;
        end
    end
    
