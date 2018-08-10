function [rslt] = ExRFPCA(X,mfd,varargin)


    parm = inputParser;
    parm.KeepUnmatched = true;
    addRequired(parm,'dat');
    addOptional(parm,'K',10);
    addOptional(parm,'mu',[]);
    addOptional(parm,'newX',[]);
    addOptional(parm,'valX',[]);
    
    parse(parm,X,mfd,varargin{:});
    pm = parm.Results;
    
    K = pm.K;
    [d,m,n] = size(X);  
        
    if isempty(pm.mu)
        mu_hat = mfd.intrinsic_mean(X);
    else
        mu_hat = pm.mu;
    end
    
    V = log_map(mu_hat,X,mfd);
    [phi,lam,totalvar] = mFPCA(V,K);
    Xi = zeros(n,K);
    XX = cell(1,K);
    VV = cell(1,K);
    
    if ~isempty(pm.newX)
        newn = size(pm.newX,3);
        XiNew = zeros(newn,K);
        XXNew = cell(1,K);
        VVNew = cell(1,K);
        VNew = log_map(mu_hat,pm.newX,mfd);
    end
    
    if ~isempty(pm.valX)
        valn = size(pm.valX,3);
        XiVal = zeros(valn,K);      
        XXVal = cell(1,K);
        VVVal = cell(1,K);
        VVal = log_map(mu_hat,pm.valX,mfd);

    
    end
    
    for k = 1:K
        VV{k} = zeros(d,m,n);
        phik = phi(:,:,k); % d-by-m
        for i = 1:n
            Xi(i,k) = project(V(:,:,i),phik);
            if k == 1
                VV{k}(:,:,i) = Xi(i,k) * phik;
            else
                VV{k}(:,:,i) = VV{k-1}(:,:,i) + Xi(i,k)*phik;
            end
        end 
        XX{k} = exp_map(mu_hat,VV{k},mfd);
        
        if ~isempty(pm.newX)
            VVNew{k} = zeros(d,m,newn);
            for i = 1:newn
                XiNew(i,k) = project(VNew(:,:,i),phik);
                if k == 1
                    VVNew{k}(:,:,i) = XiNew(i,k) * phik;
                else
                    VVNew{k}(:,:,i) = VVNew{k-1}(:,:,i) + XiNew(i,k)*phik;
                end
            end
            XXNew{k} = exp_map(mu_hat,VVNew{k},mfd);
        end
        
        if ~isempty(pm.valX)
            VVVal{k} = zeros(d,m,valn);
            for i = 1:valn
                XiVal(i,k) = project(VVal(:,:,i),phik);
                if k == 1
                    VVVal{k}(:,:,i) = XiVal(i,k) * phik;
                else
                    VVVal{k}(:,:,i) = VVVal{k-1}(:,:,i) + XiVal(i,k)*phik;
                end
            end
            XXVal{k} = exp_map(mu_hat,VVVal{k},mfd);
        end
    end
    
    rslt.totalvar = totalvar;
    rslt.mu = mu_hat;
    rslt.XX = XX;
    rslt.VV = VV;
    rslt.K = K;
    rslt.lam = lam;
    rslt.phi = phi;
    rslt.phiV = phi;
    rslt.Xi = Xi;
    rslt.Name = 'ExRFPCA';
    
    if ~isempty(pm.newX)
        rslt.XXNew = XXNew;
        rslt.VVNew = VVNew;
        rslt.XiNew = XiNew;
    end
    
    if ~isempty(pm.valX)
        rslt.XXVal = XXVal;
        rslt.VVVal = VVVal;
        rslt.XiVal = XiVal;
    end
    
    %% helper function

    
    function xi = project(Y,phik)
        xi = mean(sum(Y.*phik,1));
    end
    
    
    % X: D-m-n
    % V: D-m-n
    % mu: D-m
    function V = log_map(mu,X,mfd)
        V = mfd.Log(mu,X);
    end

    % X: D-m-n
    % V: D-m-n
    % mu: D-m
    function X = exp_map(mu,V,mfd)
        X = mfd.Exp(mu,V);
    end
end