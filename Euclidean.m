function mfd = Euclidean(d)

    D = d;

    % metric tensor on tangent space
    % x and y are of same dimension, 3-by-n, columns for vectors
    mfd.metric = @metric; %@(x,y) x(1,:).*y(1,:)+x(2,:).*y(2,:)-x(3,:).*y(3,:);
    
    % norm on tangent space induced by the metric
    % x: 3-by-n matrix, each column is a vector
    mfd.norm = @(x,~) sqrt(metric(x,x));
    
    % the extrinsic L2 norm in 3D Euclidean space
    mfd.l2norm = mfd.norm;

    % normal coordinate frame
    % theta and phi are vector/matrix of the same dimension
    mfd.rtheta = @(theta,phi) [0;1];
    mfd.rphi = @(theta,phi) [1;0];
       
    % convertion of Cartesian and spherical coordinates
    mfd.co_convert = @(t) t;
    
    mfd.d = d; % intrinsic dim
    mfd.D = d; % ambient dim
       
    % geodesic: geodesic starting at p with speed h, norm(h) = 1
    mfd.geodesic = @geodesic;
    
    % Exponential and Log map
    mfd.Exp = @Exp;
    mfd.Log = @Log;
    
    mfd.sample_point = @(n) 10*rand(mfd.D,n);
    
    mfd.sample_tangent = @sample_tangent;
    
    mfd.orthonormal_frame = @orthonormal_frame;
    
    mfd.coef_process = @(~,t) t;
    
    mfd.dist = @(P,Q) mfd.norm(P-Q);
    
    mfd.intrinsic_mean = @intrinsic_mean;
    
    mfd.coef_to_log = @coef_to_log;
    
    mfd.project = @(t) t;
    
    mfd.vfinprod_V = @vfinprod_V;
    mfd.vfinprod_Z = @vfinprod_Z;
    
    mfd.vfnorm_V = @vfnorm_V;
    mfd.vfnorm_Z = @vfnorm_Z;
    
    mfd.parallel_transport = @parallel_transport;
    
    function [met] = metric(U,V)
        met = squeeze(sum(U.*V,1));
    end
    
    % p: 3-by-1 or 3-by-m matrix, points on H2, 
    % v: a 3-by-m matrix, corresponding tangent vectors on H2
    function [emap] = Exp(p,v)
        n = size(v,3);
        if n == 1
            if size(p,2) == 1
                p = repmat(p,1,size(v,2));
            end
            emap = p + v;
        else
            emap = repmat(p,1,1,n) + v;
        end
    end

    % the Log of q's at p
    % pp: 3-by-1 vector or 3-by-m
    % qq: 3-by-m-n matrix
    % hh: 3-by-m-n matrix
    function [hh] = Log(pp,qq)
        n = size(qq,3);
        if n == 1
            if size(pp,2) == 1
                pp = repmat(p,1,size(v,2));
            end
            hh = qq - pp;
        else
            hh = qq - repmat(pp,1,1,n);
        end
    end

    function [gd] = geodesic(t,p,h)
        % p and h are of same dimension, 3-by-n
        % t is 1-by-n matrix or a scalar, treated as repmat(t,1,n)
        % the result is a 3-by-n matrix
        n = size(p,2);
        if n == 1
            n = size(t,2);
            p = repmat(p,1,n);
            h = repmat(h,1,n);
        end
        
        if length(t) == 1
            t = repmat(t,1,n);
        end
        gd = p + repmat(t,mfd.D,1).* h;
    end

    % p: 3-by-n matrix, each column is a point on H2
    % e1,e2:  3-by-n matrix, orthnormal basis at each p(:,i)
    function [frame] = orthonormal_frame(p)
        n = size(p,2);
        frame = zeros(d,n,d);
        for j = 1:d
            frame(j,:,j) = ones(1,n);
        end
    end

    % v: 3-by-n
    function v = sample_tangent(~,n)
        v = rand(mfd.D,n);
    end

    

    % X: 3-m-n
    function [Cmu,Smu] = intrinsic_mean(X)
        Cmu = squeeze(mean(X,3));
        Smu = Cmu;
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

function [r] = vfnorm_V(U,S,W)
        r = sqrt(vfinprod_V(U,U,S,W));
    end

    % U,V: same dimension, D-m-n 
    % r  : 1-n
    function [r] = vfinprod_V(U,V,~,~)
        n = size(U,3);
        r = zeros(1,n);

        if n == 1
            [A] = metric(U,V);
            r = mean(A);
        else
        
            for i = 1:n
                [A] = metric(U(:,:,i),V(:,:,i));
                r(i) = mean(A,1);
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

    function [R] = parallel_transport(~,~,U)
        R = U;
    end

    % project a point in R3 to H2
    % x: 3-by-1, any point in R3 (not necessary in H2);
    % y: the projected point on H2
    % g: the R2 distance between x and y
    function [y,g] = project(x)
        warning('off','all');
        if h2.metric(x,x) == -1
            y = x;
            g = 0;
        else
            f = @(u) (u(1)-x(1))^2+(u(2)-x(2))^2 + (sqrt(u(1)^2+u(2)^2+1)-x(3))^2;
            options = optimoptions('fminunc','Algorithm','trust-region','Display','off');
            x0 = [1;1;sqrt(3)];
            [y,g] = fminunc(f,x0,options);
        end
    end

end