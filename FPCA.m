function [rslt] = FPCA(dat,mfd,varargin)


    % dat.X:  m-n-3

    parm = inputParser;
    parm.KeepUnmatched = true;
    addRequired(parm,'dat');
    addOptional(parm,'K',10);
    addOptional(parm,'newX',[]);
    
    parse(parm,dat,varargin{:});
    pm = parm.Results;
    
    K = pm.K;
    maxK = pm.K;
    X = dat.X;
    [d,m,n] = size(X);  
        
    mu_hat = mean(X,3);
    
    V = X - repmat(mu_hat,1,1,n);
    
    [phi,lam] = line_up_fPCA(V);
    phi = phi(:,:,1:maxK);
    lam = lam(1:maxK);
    
    Xi = zeros(n,K);
    
    if ~isempty(pm.newX)
        newn = size(pm.newX,3);
        XiNew = zeros(newn,K);
        VVNew = cell(1,K);
        VNew = pm.newX - repmat(mu_hat,1,1,newn);
        XXNew = cell(1,K);
    end
    
    VV = cell(1,K);
    XX = cell(1,K);
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
        XX{k} = repmat(mu_hat,1,1,n) + VV{k};
        
        if ~isempty(pm.newX)
            VVNew{k} = zeros(d,m,newn);
            for i = 1:newn
                XiNew(i,k) = project(VNew(:,:,i),phik);
                if k == 1
                    VVNew{k}(:,:,i) = XiNew(i,k) * phik;
                else
                    VVNew{k}(:,:,i) = VVNew{k-1}(:,:,i) + XiNew(i,k) * phik;
                end
            end
            XXNew{k} = repmat(mu_hat,1,1,newn) + VVNew{k};
        end
  
    end
    
    if ~isempty(mfd)
        rslt.XX = project_to_mfd(XX,mfd,n);
    else
        rslt.XX = XX;
    end
    
    if ~isempty(pm.newX)
        if ~isempty(mfd)
            rslt.XXNew = project_to_mfd(XXNew,mfd,newn);
        else
            rslt.XXNew = XXNew;
        end
        rslt.XiNew = XiNew;
    end
    
    rslt.VV = VV;
    rslt.K = K;
    rslt.lam = lam;
    rslt.phi = phi;
    rslt.phiV = phi;
    rslt.Xi = Xi;
    rslt.Name = 'EuFPCA';
    rslt.mu = mu_hat;

    
    rslt.VVNew = VVNew;
    
    if ~isempty(mfd)
        for j = 1:m
            mu_hat(:,j) = mfd.project(mu_hat(:,j));
        end
    end
    rslt.mu = mu_hat;
    
    %% helper function

    
    function xi = project(Y,phik)
        xi = mean(sum(Y.*phik,1));
    end

    function [XX] = project_to_mfd(VV,mfd,n)
        XX = cell(1,K);
        for kk = 1:K
            XX{kk} = zeros(d,m,n);
            for jj = 1:m
                for ii = 1:n
                    u = mfd.project(VV{kk}(:,jj,ii));
                    XX{kk}(:,jj,ii) = u;
                 
                end
            end
        end
    end
end