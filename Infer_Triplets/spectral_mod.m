function newData = spectral_mod(M,kks)
    % returns a data matrix M projected onto the first 'kks' PCs. If 'kks'
    % is not given, M is randomized, and will return the number of
    % components with eigenvalue greater than the largest randomized
    % eigenvalue
    [coef,~,~,~,oexp,~] = pca(M');
    size(coef)
    if nargin==2
        ncomp = kks;
    else
        rm = row_randomized_matrix(M);
        [~,~,~,~,rexp,~] = pca(rm');
        ncomp = sum(oexp>rexp(1));
    end
    newData = M'*coef(:,1:ncomp);
end


function rm = row_randomized_matrix(M)
    [mrow,mcol] = size(M);
    rowIndex = repmat((1:mrow)',[1 mcol]);
    [~,rcolind] = sort(rand(mrow,mcol),2);
    newLinInd = sub2ind([mrow,mcol],rowIndex,rcolind);
    rm = M(newLinInd);
end
