function [varphi,lambda,totalvar,score,mu] = mFPCA(X,maxK,centered)

% X: vector-valued process, d-m-n matrix, d dimensional, m observations, n
% sample curves

% phi: d-by-m-by-K matrix, eigenfunctions
% lam: 1-by-K, eigenvalues
% score: n-by-K matrix

    if nargin == 2
        centered = true;
    end

    [d,m,n] = size(X);
    K = min(m,n) * d;
    
    if ~centered
        mu = mean(X,3);
        X = X - repmat(mu,1,1,n);
    else
        mu = zeros(d,m);
    end
    
    phi = cell(1,d);
    lam = cell(1,d);
    Xi = cell(1,d);
    for j = 1:d
        [a,b,c] = uFPCA(squeeze(X(j,:,:)));
        phi{j} = a;
        lam{j} = b;
        Xi{j} = c;
    end
    
    Z = zeros(K,K);
    
    for j = 1:d
        for k = 1:d
            bridx = ((j-1)*m+1):(j*m);
            bcidx = ((k-1)*m+1):(k*m);
            Z(bridx,bcidx) = Xi{j}' * Xi{k} / n;
        end
    end
    
    [C,T] = svd(Z);
    lambda = diag(T);
    totalvar = sum(lambda);
    
    if nargin == 1
        maxK = K;
    end
    
    lambda = lambda(1:maxK);
    varphi = zeros(d,m,maxK);
    
    for k = 1:maxK
        cm = C(:,k);
        for j = 1:d
            idx = ((j-1)*m+1):(j*m);
            varphi(j,:,k) = sum(repmat(cm(idx)',m,1).*phi{j},2);
        end
    end
    
    if nargout >= 4
        score = zeros(n,maxK);
        for k = 1:maxK
            cm = C(:,k);
            for j = d
                idx = ((j-1)*m+1):(j*m);
                score(:,k) = score(:,k) + sum(repmat(cm(idx)',n,1).*Xi{j},2);
            end
        end
    end
    
    %% helper functions
    
    % univariate fpca
    % X: m-by-n matrix, already centered
    % phi: m-by-K matrix, eigenfunctions
    % lam: 1-by-K vector, eigenvalues
    % Xi:  n-by-K matrix, fpc scores
    function [phi,lam,Xi] = uFPCA(X) 
        [mm,nn] = size(X);
        %KK = min(mm,nn);
        XX = X * X' / (nn*mm);
        [U,S] = svd(XX);
        lam = diag(S); 
        phi = U * sqrt(mm);
        Xi = X' * phi / mm;

    end
    
end