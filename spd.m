% SPD manifold with Log-Euclidean metric
% matD: the dimension of the matrices

% parameterization: see the manuscript.
% each matrix is vectorized row by row

% NOTE: all public interface must represent a matrix by D-1, and all
% interfaces must accept both D-1 or matD-matD format

function mfd = spd(matD,varargin)

    parm = inputParser;
    parm.KeepUnmatched = true;
    addRequired(parm,'matD');

    addOptional(parm,'WCache',[]);
    addOptional(parm,'DlogCache',[]);
    addOptional(parm,'DexpCache',[]);
    
    parse(parm,matD,varargin{:});
    pm = parm.Results;

    %% private interface
    D = matD*matD;
    d = matD+matD*(matD-1)/2;
    [BM,BV] = basis_sym();
    
    % cache for speedup
    cacheSize = 5;
    if ~isempty(pm.WCache)
        WCache = pm.WCache;
    else
        WCache = Cache(cacheSize);
    end
    if ~isempty(pm.DlogCache)
        DlogCache = pm.DlogCache;
    else
        DlogCache = Cache(cacheSize);
    end
    if ~isempty(pm.DexpCache)
        DexpCache = pm.DexpCache;
    else
        DexpCache = Cache(cacheSize);
    end
    
    % private functions
    mfd.issym = @issym;
    mfd.isspd = @isspd;
    mfd.count_matrix = @count_matrix;
    mfd.frobenius = @frobenius;
    mfd.vec_to_mat = @vec_to_mat;
    mfd.mat_to_vec = @mat_to_vec;
    mfd.coef = @coef;
    mfd.coef_to_sym = @coef_to_sym;
    mfd.getW = @getW;
    mfd.left_translate = @left_translate;
    mfd.inv = @myInv;
    mfd.metric_workhorse = @metric_workhorse;
    mfd.BM = BM;
    mfd.BV = BV;
    mfd.WCache = WCache;
    mfd.DlogCache = DlogCache;
    mfd.DexpCache = DexpCache;
    mfd.Dexp = @Dexp;
    mfd.Dlog = @Dlog;
    mfd.matrix_exp = @matrix_exp;
    mfd.matrix_log = @matrix_log;
    
    %% public interface

    mfd.name = 'SPD';
    mfd.tau = Inf;
    
    % metric tensor on tangent space
    % x and y are of same dimension, D-by-n, columns for vectors
    mfd.metric = @metric; %@(x,y) x(1,:).*y(1,:)+x(2,:).*y(2,:)-x(3,:).*y(3,:);
    mfd.metricZ = @metricZ;
    
    % norm on tangent space induced by the metric
    % x: 3-by-n matrix, each column is a vector
    
    mfd.norm = @mynorm;
    
    % the extrinsic L2 norm in 3D Euclidean space
    mfd.l2norm = @frobenius;
    
    mfd.d = d; % intrinsic dim
    mfd.D = D; % ambient dim
       
    % geodesic: geodesic starting at p with speed h, norm(h) = 1
    mfd.geodesic = @geodesic;
    
    % Exponential and Log map
    mfd.Exp = @Exp;
    mfd.Log = @Log;
    
    mfd.sample_point = @sample_point;
    
    mfd.sample_tangent = @sample_tangent;
    
    mfd.orthonormal_frame = @orthonormal_frame;
    
    mfd.coef_process = @coef_process;
    
    mfd.coef_to_log = @coef_to_log; % coef process to log process wrt oframe
    
    mfd.dist = @(P,Q) distance(P,Q);
    
    mfd.intrinsic_mean = @intrinsic_mean;
    
    mfd.project = @project;
    
    mfd.vfinprod_V = @vfinprod_V;
    mfd.vfinprod_Z = @vfinprod_Z;
    
    mfd.vfnorm_V = @vfnorm_V;
    mfd.vfnorm_Z = @vfnorm_Z;
    
    mfd.parallel_transport = @parallel_transport;
    
    %% helper function
    
    
    
    % inverse matrices given by X
    % X: D-m or matD-matD-m
    % Y: matD-matD-m
    function [Y] = myInv(X)
        m = count_matrix(X);
        X = vec_to_mat(X);
        Y = zeros(matD,matD,m);
        for i = 1:m
            Y(:,:,i) = inv(X(:,:,i));
        end
    end
    
    % get inner product matrix for the metric at X
    % X: points on spd, D-m or matD-matD-m
    % W: d-d-m matrix
    function [W] = getW(X)
        m = count_matrix(X);
        X = vec_to_mat(X);
        W = zeros(d,d,m);
        for i = 1:m
            for j = 1:d
                for k = j:d
                   W(j,k,i) = metric_workhorse(BM(:,:,j),BM(:,:,k),X(:,:,i));
                   W(k,j,i) = W(j,k,i); % W is symmetric
                end
            end
        end
    end

    
    % X,Y: D-m matrix or matD-matD-m matrix, spds
    % R  : 1-m matrix
    function [R] = distance(X,Y)
        R = frobenius(matrix_log(X)-matrix_log(Y));
    end
    
    % metric of tangent vectors wrt an oframe
    % Z1,Z2: coef process, same dimension, d-by-m
    % M    : 1-by-m matrix
    function [M] = metricZ(Z1,Z2)
        M = sum(Z1.*Z2,1);
    end

    % X,Y: D-by-m matrix, or matD-matD-m matrix, if only X, return norm
    % nm: 1-m matrix
    function [nm] = frobenius(X,Y)
        if nargin == 1
            Y = X;
        end
        if size(X,1) == D
            nm = sum(X.*Y,1);
        else
            nm = squeeze(sum(sum(X.*Y,1),2))';
        end
        
        if nargin == 1
            nm = sqrt(nm);
        end
    end

    % the workhorse of metric function
    % U,V: D-by-m matrices; or matD-matD-m matrices
    % S  : D-by-n matric or matD-matD-n, n=1 or n=m; the base point
    % W  : inner product matrices, d-d-n matrix
    % M  : 1-by-m matric
    function [M] = metric_workhorse(U,V,S)
        if nargin <= 2
            S = [];
        end
        
        
        if isempty(S) % S is treated as identity
            if size(U,1) == D
                M = sum(U.*V,1);
            else
                M = squeeze(sum(sum(U.*V,1),2))';
            end
        else
            U = vec_to_mat(U);
            V = vec_to_mat(V);
            S = vec_to_mat(S);

            m = count_matrix(U);
            same_S = (count_matrix(S)==1); % for speedup
                
            if same_S
                DSlog = Dlog_func(S);
            end
            M = zeros(1,m);
            for i = 1:m
                if ~same_S
                    DSlog = Dlog_func(S(:,:,i));
                end
                M(i) = frobenius(DSlog(U(:,:,i)),DSlog(V(:,:,i)));
            end
            
        end
    end
    

    
    % U,V: D-by-m matrices; or matD-matD-m matrices
    % S  : D-by-n matric or matD-matD-n, n=1 or n=m; the base point
    % W  : inner product matrices, d-d-n matrix
    % M  : 1-by-m matric
    function [M,W] = metric(U,V,S,W)
        if nargin <= 2
            S = [];
        end
        
        if nargin <= 3
            W = [];
        end
        
        
        if isempty(S) % S is treated as identity
            if size(U,1) == D
                M = sum(U.*V,1);
            else
                M = squeeze(sum(sum(U.*V,1),2))';
            end
        else
            
            if isempty(W) % need to compute inner product matrix
                
                W = WCache.get(S);
                if isempty(W)
                    W = getW(S);
                    WCache.add(S,W);
                end
            end
            
            U = vec_to_mat(U);
            V = vec_to_mat(V);
            S = vec_to_mat(S);

            m = count_matrix(U);
            same_S = (count_matrix(S)==1); % for speedup
                
            
            
            M = zeros(1,m);
            for i = 1:m
                cU = coef(U(:,:,i));
                cV = coef(V(:,:,i));
                if ~same_S
                    M(i) = cU' * W(:,:,i) * cV;
                else
                    M(i) = cU' * W * cV;
                end
            end

        end
    end

    % X: D-m or matD-matD-m
    % S,W: same as metric
    % N  : 1-m matrix
    function [N] = mynorm(X,S,W)
        if nargin <= 1
            S = [];
        end
        if nargin <= 2
            W = [];
        end
        
        N = sqrt(metric(X,X,S,W));
    end
    
    % v: D-by-n matrix, one spd per column
    % mat: matD-matD-n matrix
    function [mat] = vec_to_mat(v)
        
        if size(v,1) == matD
            mat = v;
            return;
        end
        
        [~,n] = size(v);
        mat = zeros(matD,matD,n);
        for k = 1:matD
            idx = (matD*(k-1)+1):(matD*k);
            mat(k,:,:) = v(idx,:);
        end
    end

    % mat: matD-matD-n matrix, 
    % v:   D-by-n matrix, one spd per column
    function [v] = mat_to_vec(mat)
        if size(mat,1) == D
            v = mat;
            return;
        end
        n = size(mat,3);
        v = zeros(D,n);
        for k = 1:matD
            idx = (matD*(k-1)+1):(matD*k);
            v(idx,:) = mat(k,:,:);
        end
    end
    
    % V: D-by-n matrix, one spd per column; or matD-matD-n matrix
    % X: matD-matD-n matrix
    % Y: D-by-n matrix
    function [Y,X] = matrix_exp(V)
        if size(V,1) == D
            V = vec_to_mat(V);
        end
        n = size(V,3);
        X = zeros(matD,matD,n);
        for i = 1:n
            X(:,:,i) = expm(V(:,:,i));
        end
        Y = mat_to_vec(X);
    end


    % V: D-by-n matrix, one spd per column; or matD-matD-n matrix
    % X: matD-matD-n matrix
    % Y: D-by-n matrix
    function [Y,X] = matrix_log(V)
        if size(V,1) == D
            V = vec_to_mat(V);
        end
        n = size(V,3);
        X = zeros(matD,matD,n);
        for i = 1:n
            X(:,:,i) = logm(V(:,:,i));
        end
        Y = mat_to_vec(X);
    end
    

    % Riemannian Exponential map, p and v of the same dimension, one
    % tangent vector per point
    % p: D-m or matD-matD-m spds,  
    % v: D-by-m or matD-matD-m symmetric matrices, corresponding tangent
    % vectors on spd
    % E: D-m matrix
    function [E] = Exp_single(p,v)
        
        F = get_Dlog(p);
        
        NP = count_matrix(p);
        if NP == 1
            m = count_matrix(v);
            p = replicate(p,m);
        else
            m = NP;
        end
        
        E = zeros(D,m);
        for i = 1:m
            if NP == 1
                f = F{1};
            else
                f = F{i};
            end
            S = p(:,i);
            L = v(:,i);
            E(:,i) = matrix_exp(matrix_log(S)+mat_to_vec(f(L)));
        end
    end

    % Riemannian Exponential map, p and v of the same dimension, multiple
    % tangent vectors per point, n > 1
    % p: D-m or matD-matD-m spds,  
    % v: D-by-m-n or matD-matD-m symmetric matrices, corresponding tangent
    % vectors on spd
    % F: the Dlog functions, one per point
    % E: D-m-n matrix
    function [E] = Exp_multiple(p,v,m,n,F)
        
        E = zeros(D,m,n);
        for j = 1:m
            f = F{j};
            S = p(:,j);
            L = v(:,j,:);
            logS = matrix_log(S);
            E(:,j,:) = matrix_exp(repmat(logS,1,n)+mat_to_vec(f(L)));
        end
    end


    % Riemannian Exponential map, allow n > 1 tangent vectors per point
    % p: D-m or matD-matD-m spds,  
    % v: D-by-m-n or matD-matD-m-n symmetric matrices, corresponding tangent
    % vectors on spd, n could be 1 or >1
    % E: D-m-n matrix
    function [E] = Exp(p,v)
        
        if size(v,1) == D
            n = size(v,3);
        else
            n = size(v,4);
        end
        
        if n == 1
            E = Exp_single(p,v);
        else

            F = get_Dlog(p);
            m = count_matrix(p);
            E = Exp_multiple(p,v,m,n,F);
        end
        
    end

    % the Log of q's at p
    % pp: D-by-k or matD-matD-k matrix, k = 1 or k = m
    % qq: D-by-m-n or matD-matD-m-n matrix
    % L : D-m-n matrix if n>1, or D-m matrix
    function [L] = Log(pp,qq)
        
        if size(qq,1) == D
            n = size(qq,3);
            m = size(qq,2);
        else
            n = size(qq,4);
            m = size(qq,3);
        end
        
        if n == 1
            L = Log_single(pp,qq);
        else
            L = Log_multiple(pp,qq,m,n);
        end
    end

    
    % the Log of q's at p, allow two cases:
    %     1) m = 1, so pp is a single point, and qq could be m points
    %     2) m > 1, then qq and pp are of the same dimension
    % pp: D-by-m or matD-matD-m matrix, m = 1 or m = n
    % qq: D-by-n or matD-matD-n matrix
    % L : D-n matrix
    function [L] = Log_single(pp,qq)
        
        %F = get_Dexp(pp);
        
        NP = count_matrix(pp);
        if NP == 1
            n = count_matrix(qq);
            pp = replicate(pp,n);
        else
            n = NP;
        end
       
        %pp = vec_to_mat(pp);
        %qq = vec_to_mat(qq);
        logS = matrix_log(pp);
        logQ = matrix_log(qq);
        F = get_Dexp(logS);
        L = zeros(matD,matD,n);
        for i = 1:n
            if NP == 1
                f = F{1};
            else
                f = F{i};
            end
            L(:,:,i) = f(logQ(:,i)-logS(:,i));
        end
        L = mat_to_vec(L);
    end

    % the Log of q's at p, require: n > 1
    % pp: D-by-m-n or matD-matD-m-n matrix, m = 1 or m = n
    % qq: D-by-m-n or matD-matD-m matrix
    % L : D-m-n matrix
    function [L] = Log_multiple(pp,qq,m,n)
        
        logS = matrix_log(pp);
        
        F = get_Dexp(logS);
        L = zeros(D,m,n);
        for j = 1:m
            f = F{j};
            if size(qq,1) == D
                logQ = matrix_log(squeeze(qq(:,j,:)));
            else
                logQ = matrix_log(squeeze(qq(:,:,j,:)));
            end
            L(:,j,:) = mat_to_vec(f(logQ-repmat(logS(:,j),1,n)));
        end
    end

    % X: D-n or matD-matD-n matrix
    % n: the number of matrices given by X
    function [n] = count_matrix(X)
        if size(X,1) == D
            n = size(X,2);
        else
            n = size(X,3);
        end
    end

    % replicate a matrix n times
    % x: a matrix given by D-m or matD-matD-m, m=1 or m=n
    % n: the number of times to be replcated, do nonthing if m=n
    % M: D-n or matD-matD-n matrices
    function [M] = replicate(x,n)
        n0 = count_matrix(x);
        if n0 == n
            M = x;
        else
            if n0 > 1;  error(' #x must be 1 or n ');  end;
            if size(x,1) == D
                M = repmat(x,1,n);
            else
                M = repmat(x,1,1,n);
            end
        end
    end

    % basis of m-by-m symmetric matrices
    % NOTE: the basis might not be orthonormal!!!!!
    % BM: m-m-n matrix, n = m+m(m-1)/2, the dim(Sym(m))
    % BV: (m*m)-n matrix, the vectorized version of B1
    function [BM,BV] = basis_sym()
        m = matD;
        Dsym = m + m*(m-1)/2;
        BM = zeros(m,m,Dsym);
        k = 0;
        for i = 1:m
            for j = i:m
                k = k + 1;
                BM(i,j,k) = 1;
                BM(j,i,k) = 1;
            end
        end
        BV = vec_to_mat(BM);
    end

    % coefficient of V wrt BM
    % V: D-m matrix or matD-matD-m
    % C: d-m matrix
    function [C] = coef(V)
        V = vec_to_mat(V);
        m = size(V,3);
        C = zeros(d,m);
        for n = 1:m
            k = 0;
            for i = 1:matD
                for j = i:matD
                    k = k + 1;
                    C(k,n) = V(i,j,n);
                end
            end
        end
    end

    % C: d-by-m matrix
    % V: matD-matD-m matrix
    function [V] = coef_to_sym(C)
        m = size(C,2);
        V = zeros(matD,matD,m);
        for n = 1:m
            k = 0;
            for i = 1:matD
                for j = i:matD
                    k = k + 1;
                    V(i,j,n) = C(k,n);
                    V(j,i,n) = C(k,n);
                end
            end
        end
    end

    % the differential of matrix_log at V
    % X: D-1 or matD-matD matrix, spd
    % M: d-by-d matrix representation of D(exp) wrt basis B
    function [M] = Dlog(X)
        logX = matrix_log(X);
        A = Dexp(logX);
        M = inv(A);
    end

    % X: D-m or matD-matD-m
    % F: cell array of function handles of length m
    function [F] = get_Dlog(X)
        F = DlogCache.get(X);
        if isempty(F)
            m = count_matrix(X);
            F = cell(1,m);
            X = mat_to_vec(X);
            for j = 1:m
                F{j} = Dlog_func(X(:,j));
            end
            DlogCache.add(X,F);
        end
    end
    
    % X: D-1 or matD-matD
    function [f] = Dlog_func(X)
        M = Dlog(X);
        f = @(S) coef_to_sym(M*coef(S));
    end

    % X: D-m or matD-matD-m
    % F: cell array of function handles of length m
    function [F] = get_Dexp(X)
        F = DexpCache.get(X);
        if isempty(F)
            m = count_matrix(X);
            F = cell(1,m);
            X = mat_to_vec(X);
            for j = 1:m
                F{j} = Dexp_func(X(:,j));
            end
            DexpCache.add(X,F);
        end
    end

    
    % X: D-1 or matD-matD
    function [f] = Dexp_func(X)
        M = Dexp(X);
        f = @(S) coef_to_sym(M*coef(S));
    end

    % the differential of matrix_exp at V
    % V: D-1 or matD-matD matrix, symmetric
    % M: d-by-d matrix representation of D(exp) wrt basis B
    function [M] = Dexp(V)
        V = vec_to_mat(V);
        V = (V+V') /2; % make sure symmetric
        [U,S] = eig(V);
        K = 25; % make it larger to have a better approximation
        M = zeros(d,d);
        for j = 1:d
            R = zeros(matD,matD);
            phij = BM(:,:,j);
            kfact = 1;
            for k = 1:K
                kfact = kfact * k;
                W = zeros(matD,matD);
                Q = U'*phij*U;
                for ell = 0:(k-1)
                    diagS = diag(S);
                    W = W + diag(diagS.^(k-ell-1)) * Q * diag(diagS.^(ell));
                end
                R = R + W / (kfact);
            end
            R = U * R * U';
            R = (R + R') / 2;
            M(:,j) = coef(R);
        end
        
    end

    % left translation of Y by X
    % X: D-m matrix or matD-matD-m matrix, if m=1, treated as replicated n,
    % if m > 1, then it must be m = n;
    % Y: D-n matrix or matD-matD-n matrix
    % L: D-n matrix
    function [L] = left_translate(X,Y)
        n = count_matrix(Y);
        X = replicate(X,n);
        L = matrix_exp(matrix_log(X) + matrix_log(Y));
    end


    % p: spds, count_matrix(p)=1 or count_matrix(p)=count_matrix(h)
    % h: same dimension, D-m or matD-matD-m
    % t  : 1-by-m or scalar, treated as repmat(t,1,m)
    % gd : D-m matrix
    function [gd] = geodesic(t,p,h)
        
        %n = count_matrix(h);
        if length(t) == 1
            tU = t*h;
        else
            m = length(t);
            if count_matrix(h) == 1
                h = replicate(h,m);
            end
            tU = repmat(t,D,1) .* mat_to_vec(h);
        end
        tU = vec_to_mat(tU);
           
        p = vec_to_mat(p);
        p = replicate(p,m);    
        gd = left_translate(p,matrix_exp(tU));
        gd = mat_to_vec(gd);
    end

    % turn a frame into an orthonormal one
    % frame: D-d matrix or matD-matD-d matrix
    % p: a spd, D-1 or matD-matD
    % W: inner product matrix, d-d
    % G: matD-matD-d matrix
    function [G] = gram_schmidt(frame,p,W)
        F = vec_to_mat(frame);
        G = zeros(matD,matD,d);
        for j = 1:d
            g = F(:,:,j);
            if j > 1
                M = metric(replicate(g,j-1),G(:,:,1:(j-1)),p,W);
                for k = 1:j-1
                    g = g - M(k)*G(:,:,k);
                end
            end
            g = g / mynorm(g,p,W);
            G(:,:,j) = g;
        end
    end

    % p: D-by-m or matD-matD-m matrix
    % frame: D-m-d matrix
    function [frame,mframe] = orthonormal_frame(p)
        m = count_matrix(p);
        W = WCache.get(p);
        mframe = zeros(matD,matD,m,d);
        frame = zeros(D,m,d);
        
        if isempty(W)
            W = getW(p);
            WCache.add(p,W);
        end
        for j = 1:m
            mframe(:,:,j,:) = gram_schmidt(BM,p(:,j),W(:,:,j));
            frame(:,j,:) = mat_to_vec(squeeze(mframe(:,:,j,:)));
        end
    end

    % X: D-n
    function [X] = sample_point(n)

        %low = 0.1;
        %high = 10;
        lam = 10.^(2*rand(matD,n)-1);
        
        X = zeros(matD,matD,n);
        for i = 1:n
            U = orth(rand(matD,matD));
            S = diag(lam(:,i));
            V = U*S*U';
            X(:,:,i) = (V+V')/2;
        end
        
        X = mat_to_vec(X);
    end

    % p: a spd, D-m or matD-matD
    % v: D-by-m-n
    % z: the coef of v, d-m-n matrix
    % if m == 1, it v and z are squeezed
    function [v,z] = sample_tangent(p,n,distribution)
        
        if nargin == 2
            distribution = 'uniform';
        end
        
        m = count_matrix(p);
        
        if m == 1
        
            if strcmp(distribution,'uniform')
                z = 2*rand(d,n)-1; %[-1,1]
            else
                z = randn(d,n);
            end

            [frame] = orthonormal_frame(p);
            vframe = squeeze(mat_to_vec(frame)); % D-d
            v = vframe * z;
            v = mat_to_vec(v);
        else
            if strcmp(distribution,'uniform')
                z = 2*rand(d,m,n)-1; %[-1,1]
            else
                z = randn(d,m,n);
            end

            [vframe] = orthonormal_frame(p);
            %vframe = squeeze(mat_to_vec(frame)); % D-d
            v = zeros(D,m,n);
            for i = 1:n
                for j = 1:m
                    v(:,j,i) = mat_to_vec(squeeze(vframe(:,j,:)) * z(:,j,i));
                end
            end
        end
    end

    % Z: d-m-n, the coefficient process Z(t) with respect to orth frame
    % V: D-m-n, or matD-matD-m-n, tangent vectors on spd
    % mu: D-by-m matrix or matD-matD-matrix
    function [Z] = coef_process(mu,V,oframe)
        if nargin == 2
            [oframe] = orthonormal_frame(mu); 
        end
        if size(V,1) == D
            n = size(V,3);
        else
            n = size(V,4);
        end
        
        %oframe = mat_to_vec(oframe);
        m = count_matrix(mu);
        W = WCache.get(mu);
        if isempty(W)
            W = getW(mu);
            WCache.add(mu,W);
        end
        
        mu = mat_to_vec(mu);
        
        Z = zeros(d,m,n);
        for k = 1:m
            if size(V,1) == D
                vv = squeeze(V(:,k,:)); %D-by-n
            else
                vv = squeeze(V(:,:,k,:));
            end
            for j = 1:d
                Z(j,k,:) = metric(replicate(oframe(:,k,j),n),vv,mu(:,k),W(:,:,k));
            end
        end
    end

    % X: D-m-n or matD-matD-m-n matrix, the m dimension across design
    % points t1, t2, ..., tm, n is the sample size #curves
    % M: D-m matrix
    function [M] = intrinsic_mean(X)
        if size(X,1) == D
            m = size(X,2);
            n = size(X,3);
        else
            m = size(X,3);
            n = size(X,4);
        end
        
        if n > 1
            M = zeros(D,m);
            for j = 1:m
                if size(X,1) == D
                    Y = squeeze(X(:,j,:));
                    Y = vec_to_mat(Y);
                else
                    Y = squeeze(X(:,:,j,:));
                end
                M(:,j) = matrix_exp(mean(mat_to_vec(matrix_log(Y)),2));
            end
        else
            M = matrix_exp(mean(mat_to_vec(matrix_log(X)),2));
        end
        
        M = mat_to_vec(M);
    end




    % "project" a symmetric matrix onto spd
    % X: D-m matrix or matD-matD-m matrix
    % Y: D-m matrix
    % g: 1-m matrix
    function [Y,g] = project(X)
        m = count_matrix(X);
        X = vec_to_mat(X);
        Y = zeros(matD,matD,m);

        for i = 1:m
            [U,S] = eig(X(:,:,i));
            S(S<0) = 1e-8;
            Y(:,:,i) = U*S*U';
        end
        Y = mat_to_vec(Y);
        g = sqrt(frobenius(X-Y));
    end

    % mu: D-by-m matrix
    % Z: d-by-m-by-n matrix
    % V: D-by-m-by-n matrix
    function V = coef_to_log(mu,Z)
        [frame] = orthonormal_frame(mu);
        [~,m,n] = size(Z);
        V = zeros(D,m,n);
        if n == 1
            for j = 1:d
                V = V + repmat(Z(j,:),D,1) .* frame(:,:,j);
            end
        else
            for j = 1:d
                 V = V + repmat(Z(j,:,:),D,1,1) .* repmat(frame(:,:,j),1,1,n);
            end
        end
    end

    % U  : same dimension, D-m-n or matD-matD-m-n
    % S  : base point, D-m or matD-matD-m
    % W  : d-d-m matrix, inner product space
    % r  : 1-n matrix
    function [r] = vfnorm_V(U,S,W)
        if nargin <= 1
            S = [];
        end
        if nargin  <= 2
            W = [];
        end
        r = sqrt(vfinprod_V(U,U,S,W));
    end

    % U,V: same dimension, D-m-n or matD-matD-m-n
    % S  : base point, D-m or matD-matD-m
    % W  : d-d-m matrix, inner product space
    % r  : 1-n
    function [r] = vfinprod_V(U,V,S,W)
        if nargin <= 1
            S = [];
        end
        if nargin  <= 2
            W = [];
        end
        
        if size(U,1) == D
            type = 1;
            n = size(U,3);
        else
            type = 2;
            n = size(U,4);
        end
        r = zeros(1,n);
        
        if nargin <= 2;  S = []; end
        if nargin <= 3;  W = []; end
        
        if n == 1
            [A,~] = metric(U,V,S,W);
            r = mean(A);
        else
        
            for i = 1:n
                if type == 1
                    [A,W] = metric(U(:,:,i),V(:,:,i),S,W);
                else
                    [A,W] = metric(U(:,:,:,i),V(:,:,:,i),S,W);
                end
                r(i) = mean(A);
            end
        end
    end

    % W,Z: d-m-n matrix
    % r  : 1-n
    function [r] = vfnorm_Z(W)
        r = vfinprod_Z(W,W);
    end

    % W,Z: d-m-n matrix
    % r  : 1-n
    function [r] = vfinprod_Z(Q,Z)
        n = size(Z,3);
        if n == 1
            r = mean(sum(Q.*Z,1));
        else
            r = squeeze(mean(sum(Q.*Z,1),2));
        end
    end

    % parallel transport U from P to Q
    % P,Q,U: same dim, D-m or matD-matD-m matrix
    % V    : D-m matrix
    function [V] = parallel_transport(P,Q,U)
        m = count_matrix(U);
        P = vec_to_mat(P);
        Q = vec_to_mat(Q);
        U = vec_to_mat(U);
        V = zeros(matD,matD,m);
        for i = 1:m
            if size(P,3) == 1
                p = P;
                q = Q;
            else
                p = P(:,:,i);
                q = Q(:,:,i);
            end
            a = matrix_exp(matrix_log(q)-matrix_log(p)); %p \ q;
            b = p;
            f = Dlog_func(b);
            g = Dexp_func(matrix_log(a)+matrix_log(b));
            V(:,:,i) = g(f(U(:,:,i)));
        end
        V = mat_to_vec(V);
    end

    % test whether given matrices are spds
    % X: D-m or matD-matD-m
    % b: 1-m vector
    function b = isspd(X)
        X = vec_to_mat(X);
        m = count_matrix(X);
        b = zeros(1,m);
        for i = 1:m
            lam = eig(X(:,:,i));
            b(i) = all(lam > 0);
        end
    end

    % test whether given matrices are symmetric
    % X: D-m or matD-matD-m
    % b: 1-m vector
    function b = issym(X)
        X = vec_to_mat(X);
        m = count_matrix(X);
        b = zeros(1,m);
        for i = 1:m
            b(i) = (frobenius(X(:,:,i)-X(:,:,i)') < 1e-12);
        end
    end
end