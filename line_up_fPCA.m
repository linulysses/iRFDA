function [phi,lam] = line_up_fPCA(Z)

% Z: coef-process data, d-m-n matrix, not m-by-n-by-d array
% phi: d-by-m-by-K matrix
% lam: vector
    [d,m,n] = size(Z);
    K = min(m,n);
    
    X = zeros(m*d,n);
    for j = 1:d
        i = (j-1)*m+1;
        k = j*m;
        X(i:k,:) = squeeze(Z(j,:,:));
    end
    
    XX = X * X' / (n*m);
    [U,S] = svd(XX);
    lam = diag(S);
    U = U * sqrt(m);
    
    phi = zeros(d,m,K);
    
    for k = 1:K
        pp = U(:,k);
        for j = 1:d
            idx = ((j-1)*m+1):(j*m);
            phi(j,:,k) = pp(idx);
        end
    end
end