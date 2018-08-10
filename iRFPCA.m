function [rslt] = iRFPCA(X,mfd,varargin)


    % X:  D-m-n

    parm = inputParser;
    parm.KeepUnmatched = true;
    addRequired(parm,'X');
    addRequired(parm,'mfd');
    addOptional(parm,'newX',[]); % X from test set
    addOptional(parm,'valX',[]); % X from validation 
    addOptional(parm,'K',10);
    
    parse(parm,X,mfd,varargin{:});
    pm = parm.Results;
    
    d = mfd.d;
    K = pm.K;
    [~,m,n] = size(X);

    mu_hat = mfd.intrinsic_mean(X);
    V = log_map(mu_hat,X,mfd);
    Z = mfd.coef_process(mu_hat,V); % d-m-n
    [phi,lam,totalvar] = mFPCA(Z,K);

    Xi = zeros(n,K);
    ZZ = cell(1,K);
    XX = cell(1,K);
    VV = cell(1,K);
    
    if ~isempty(pm.newX)
        newn = size(pm.newX,3);
        XiNew = zeros(newn,K);
        ZZNew = cell(1,K);
        XXNew = cell(1,K);
        VVNew = cell(1,K);
        VNew = log_map(mu_hat,pm.newX,mfd);
        ZNew = mfd.coef_process(mu_hat,VNew); % d-m-n
    
    end
    
    if ~isempty(pm.valX)
        valn = size(pm.valX,3);
        XiVal = zeros(valn,K);
        ZZVal = cell(1,K);
        XXVal = cell(1,K);
        VVVal = cell(1,K);
        VVal = log_map(mu_hat,pm.valX,mfd);
        ZVal = mfd.coef_process(mu_hat,VVal); % d-m-n
    
    end
    
    for k = 1:K
        ZZ{k} = zeros(d,m,n);
        phik = phi(:,:,k); % d-by-m
        for i = 1:n
            Xi(i,k) = project(Z(:,:,i),phik);
            if k == 1
                ZZ{k}(:,:,i) = Xi(i,k) * phik;
            else
                ZZ{k}(:,:,i) =ZZ{k-1}(:,:,i) + Xi(i,k)*phik;
            end
        end 
        VV{k} = to_log_process(mu_hat,ZZ{k},mfd);
        XX{k} = exp_map(mu_hat,VV{k},mfd);
        
        if ~isempty(pm.newX)
            ZZNew{k} = zeros(d,m,newn);
            for i = 1:newn
                XiNew(i,k) = project(ZNew(:,:,i),phik);
                if k == 1
                    ZZNew{k}(:,:,i) = XiNew(i,k) * phik;
                else
                    ZZNew{k}(:,:,i) = ZZNew{k-1}(:,:,i) + XiNew(i,k)*phik;
                end
            end
            VVNew{k} = to_log_process(mu_hat,ZZNew{k},mfd);
            XXNew{k} = exp_map(mu_hat,VVNew{k},mfd);
        end
        
        if ~isempty(pm.valX)
            ZZVal{k} = zeros(d,m,valn);
            for i = 1:valn
                XiVal(i,k) = project(ZVal(:,:,i),phik);
                if k == 1
                    ZZVal{k}(:,:,i) = XiVal(i,k) * phik;
                else
                    ZZVal{k}(:,:,i) = ZZVal{k-1}(:,:,i) + XiVal(i,k)*phik;
                end
            end
            VVVal{k} = to_log_process(mu_hat,ZZVal{k},mfd);
            XXVal{k} = exp_map(mu_hat,VVVal{k},mfd);
        end
        
    end
    
    rslt.totalvar = totalvar;
    rslt.mu = mu_hat;
    rslt.ZZ = ZZ;
    rslt.XX = XX;
    rslt.VV = VV;
    rslt.K = K;
    rslt.lam = lam;
    rslt.phi = phi;
    rslt.phiV = mfd.coef_to_log(mu_hat,phi);
    rslt.Xi = Xi;
    rslt.Name = 'iRFPCA';
    
    if ~isempty(pm.newX)
        rslt.ZZNew = ZZNew;
        rslt.XXNew = XXNew;
        rslt.VVNew = VVNew;
        rslt.XiNew = XiNew;
    end
    
    if ~isempty(pm.valX)
        rslt.ZZVal = ZZVal;
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

    % mu: 3-m, mean curve
    % Z: d-m-n, coefficient process
    % V: 3-m-n
    function V = to_log_process(mu,Z,mfd)
        [frame] = mfd.orthonormal_frame(mu); % 3-m-d
        [dd,mm,nn] = size(Z);
        V = zeros(mfd.D,mm,nn);
        for u = 1:dd
            for jj = 1:mm
                e = frame(:,jj,u); 
                V(:,jj,:) = squeeze(V(:,jj,:)) + repmat(squeeze(Z(u,jj,:))',mfd.D,1) ...
                                               .* repmat(e,1,nn);
            end
        end
    end
end