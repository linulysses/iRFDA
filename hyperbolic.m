function mfd = hyperbolic()

    %% private interface
    
    D = 3;
    d = 2;
    
    % normal coordinate frame
    % theta and phi are vector/matrix of the same dimension
    mfd.rtheta = @(theta,phi) [cosh(theta).*cos(phi); ...
              cosh(theta).*sin(phi); ...
              sinh(theta)];
    mfd.rphi = @(theta,phi) [-sinh(theta).*sin(phi);...
             sinh(theta).*cos(phi);...
             zeros(1,length(theta))];

    %% public interface
    
    mfd.name = 'H2';
    mfd.tau = Inf;

    % metric tensor on tangent space
    % x and y are of same dimension, 3-by-n, columns for vectors
    mfd.metric = @metric; %@(x,y) x(1,:).*y(1,:)+x(2,:).*y(2,:)-x(3,:).*y(3,:);
    
    % norm on tangent space induced by the metric
    % x: 3-by-n matrix, each column is a vector
    mfd.norm = @mynorm;
    
    % the extrinsic L2 norm in 3D Euclidean space
    mfd.l2norm = @(x) sqrt(dot(x,x));
    
    % the normal vector field in 3D Euclidean space 
    % theta and phi are vector/matrix of the same dimension
    mfd.l2nml = @(theta,phi) [-sinh(theta).*cos(phi);...
           -sinh(theta).*sin(phi);...
           cosh(theta)];
       
    % convertion of Cartesian and spherical coordinates
    mfd.co_convert = @co_convert_H2;
    
    mfd.d = 2; % intrinsic dim
    mfd.D = 3; % ambient dim
       
    % geodesic: geodesic starting at p with speed h, norm(h) = 1
    mfd.geodesic = @geodesic;
    
    % Exponential and Log map
    mfd.Exp = @Exp;
    mfd.Log = @Log;
    
    mfd.sample_point = @sample_H2;
    
    mfd.sample_tangent = @sample_tangent;
    
    mfd.orthonormal_frame = @orthonormal_frame;
    
    mfd.coef_process = @coef_process;
    
    mfd.coef_to_log = @ coef_to_log; % coef process to log process wrt oframe
    
    mfd.dist = @(P,Q) acosh(-metric(P,Q));
    
    mfd.intrinsic_mean = @intrinsic_mean;
    
    mfd.project = @project;
    
    mfd.vfinprod_V = @vfinprod_V;
    mfd.vfnorm_V = @(U,V,~) sqrt(vfinprod_V(U,V));
    mfd.vfinprod_Z = @vfinprod_Z;
    mfd.vfnorm_Z = @(U,V) sqrt(vfinprod_Z(U,V));
    
    mfd.parallel_transport = @parallel_transport_H2;
    
    %% helper function
    
    function [R] = mynorm(x,~)
        R = sqrt(lorentz(x,x));
    end
    
    
    function [met] = metric(U,V,~)
        met = squeeze(U(1,:,:).*V(1,:,:)+U(2,:,:).*V(2,:,:)-U(3,:,:).*V(3,:,:));
  
    end
    
    % p: 3-by-1 (when n==1) or 3-by-m matrix, points on H2, 
    % v: a 3-m-n matrix, corresponding tangent vectors on H2
    function [E] = Exp(p,v)
        
        n = size(v,3);
        if n == 1
            E = Exp_single(p,v);
        else
            m = size(v,2);
            E = Exp_multiple(p,v,m,n);
        end
    end

    % p: 3-by-1 or 3-by-n matrix, points on H2, 
    % v: a 3-by-n matrix, corresponding tangent vectors on H2
    function [emap] = Exp_single(p,v)
        NP = size(p,2);
        if NP == 1
            t = mfd.norm(v);
            z = (t == 0);
            n = size(v,2);
            emap = zeros(3,n);
            emap(:,z) = repmat(p,1,sum(z));
            if sum(~z) > 0
                emap(:,~z) = geodesic(t(1,~z),repmat(p,1,n-sum(z)),v(:,~z)./repmat(t(1,~z),3,1));
            end
        else
            t = mfd.norm(v);
            z = (t == 0);
            n = size(v,2);
            emap = zeros(3,n);
            emap(:,z) = p(:,z); 
            if sum(~z) > 0
                emap(:,~z) = geodesic(t(1,~z),p(:,~z),v(:,~z)./repmat(t(1,~z),3,1));
            end
        end
    end

    % p: 3-m, points on H2, 
    % v: 3-m-n, corresponding tangent vectors on H2, n > 1
    % E: 3-m-n matrix
    function [E] = Exp_multiple(p,v,m,n)
        E = zeros(mfd.D,m,n);
        for i = 1:n
            E(:,:,i) = Exp_single(p,v(:,:,i));
        end
    end

    % the Log of q's at p
    % pp: 3-by-1 vector or 3-by-n
    % qq: 3-by-n matrix
    % hh: 3-by-n matrix
    function [hh] = Log_single(pp,qq)
        % find h and tq such that q = cosh(t)*p + sinh(t)*h
        [Cp,Sp] = mfd.co_convert(pp);
        l2nml_p = mfd.l2nml(Sp(1,:),Sp(2,:));
        
        eps = 1e-8;
        
        NP = size(pp,2);
        
        n = size(qq,2);
        hh = zeros(3,n);
        for i = 1:n
            
            if NP == 1
                p = Cp;
                l2np = l2nml_p;
            else
                p = Cp(:,i);
                l2np = l2nml_p(:,i);
            end
            
            q = qq(:,i);
            
            w = cross(p,q);
            h = cross(l2np,w);
            h = h / sqrt(mfd.metric(h,h));
            np = p;
            if abs(mfd.metric(h,np)) > eps || abs(dot(h,w)) > eps
                warning('h is not in tangent space');
            end

            px = p(1); hx = h(1); qx = q(1); 
            [k1,k2] = quadroot(px+hx,-2*qx,px-hx);
            if k2 < 0 && k1 < 0
                warning('the solution is negative');
            end

            if k2 > 0;   t2 = log(k2);  else  t2 = 0;  end
            if k1 > 0;   t1 = log(k1);  else  t1 = 0;  end

            
            delta1 = mfd.l2norm(geodesic(t1,p,h)-q);
            delta2 = mfd.l2norm(geodesic(t2,p,h)-q);
            
            if delta1 < delta2
                tq = t1;
            else
                tq = t2;
            end

            hh(:,i) = tq * h;
        end

    end


    % the Log of q's at p
    % pp: 3-by-m vector or 3-by-1 (only when n==1)
    % qq: 3-by-m-n matrix
    % hh: 3-by-m-n matrix
    function [hh] = Log(pp,qq)
        n = size(qq,3);
        if n == 1
            hh = Log_single(pp,qq);
        else
            m = size(qq,2);
            hh = zeros(D,m,n);
            for i = 1:n
                hh(:,:,i) = Log_single(pp,qq(:,:,i));
            end
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

        gd = repmat(cosh(t),3,1).*p + repmat(sinh(t),3,1).*h;
    end

    % p: 3-by-n matrix, each column is a point on H2
    % e1,e2:  3-by-n matrix, orthnormal basis at each p(:,i)
    function [frame] = orthonormal_frame(p)
       
        n = size(p,2);
        e1 = zeros(3,n);
        e2 = zeros(3,n);
        for i = 1:n
            if abs(p(3,i)-1) < 1e-12
                e1(:,i) = [1;0;0];
                e2(:,i) = [0;1;0];
            else
                [~,s] = mfd.co_convert(p(:,i));
                e1(:,i) = mfd.rtheta(s(1),s(2));
                e1(:,i) = e1(:,i) ./ mfd.norm(e1(:,i));
                e2(:,i) = mfd.rphi(s(1),s(2));
                e2(:,i) = e2(:,i) ./ mfd.norm(e2(:,i));
            end  
        end
        frame = zeros(3,n,2);
        frame(:,:,1) = e1;
        frame(:,:,2) = e2;
    end

    % v: 3-by-n
    function v = sample_tangent(p,n,distribution)
        
        if nargin == 2
            distribution = 'uniform';
        end
        
        [frame] = orthonormal_frame(p);
        e1 = frame(:,:,1);
        e2 = frame(:,:,2);
        if mod(n,2) == 0
            nn = n/2;
            if strcmp(distribution,'uniform')
                coef = 2*rand(2,nn)-1;
            else
                coef = randn(2,nn);
            end
            v1 = repmat(coef(1,:),3,1).*repmat(e1,1,nn) ...
                + repmat(coef(2,:),3,1).*repmat(e2,1,nn);
            v2 = repmat(coef(1,:),3,1).*repmat(-e1,1,nn) + ...
                 repmat(coef(2,:),3,1).*repmat(-e2,1,nn);
            v = [v1 v2];
        else
            nn = n;
            if strcmp(distribution,'uniform')
                coef = 2*rand(2,nn)-1;
            else
                coef = randn(2,nn);
            end
            v = repmat(coef(1,:),3,1).*repmat(e1,1,nn) ...
                + repmat(coef(2,:),3,1).*repmat(e2,1,nn);
        end
    end

    % Z: 2-m-n, the coefficient process Z(t) with respect to frame [e1,e2]
    % V: 3-m-n, tangent vectors on H2
    % mu: 3-by-m matrix
    function [Z] = coef_process(mu,V)
        [frame] = orthonormal_frame(mu); %3-by-m-by-d
        m = size(mu,2);
        n = size(V,3);
        d = 2;
        Z = zeros(d,m,n);
        for k = 1:m
            vv = squeeze(V(:,k,:)); %3-by-n
            for j = 1:d
                Z(j,k,:) = mfd.metric(repmat(frame(:,k,j),1,n),vv);
            end
        end
    end

    % X: 3-m-n
    function [Cmu,Smu] = intrinsic_mean(X)
        [~,m,~] = size(X);
        Cmu = zeros(3,m);
        Smu = zeros(2,m);
        for j = 1:m
            P = squeeze(X(:,j,:));
            [c,s] = frechet_mean_H2(P);
            Cmu(:,j) = c;
            Smu(:,j) = s;
        end
    end




    % project a point in R3 to H2
    % x: 3-by-1, any point in R3 (not necessary in H2);
    % y: the projected point on H2
    % g: the R2 distance between x and y
    function [y,g] = project(x)
        warning('off','all');
        if mfd.metric(x,x) == -1
            y = x;
        else
            if x(1) == 0 && x(2) == 0 % x is on the z-axis
                if x(3) <= 2
                    y = [0;0;1];
                else
                    y(3) = x(3)/2;
                    y(1) = 0;
                    y(2) = sqrt(y(3)^2-1);
                end
                y = y(:);
            elseif x(1) == 0 || x(2) == 0
                if x(1) == 0
                    a = x(2);  b = x(3);
                else
                    a = x(1);  b = x(3);
                end
                rt = roots([4 4*a 4+a^2-b^2 -4*a a^2]);
                fval = a^2 + (1-b)^2;
                s = [0 1];
                for j = 1:length(rt)
                    if isreal(rt(j))
                        tmp = (rt(j)-a)^2+(sqrt(rt(j)^2+1)-b)^2;
                        if tmp < fval
                            fval = tmp;
                            s = [rt(j) sqrt(rt(j)^2+1)];
                        end
                    end
                end
                if x(1) == 0
                    y(1) = 0;
                    y(2:3) = s';
                else
                    y(2) = 0;
                    y(1) = s(1);
                    y(3) = s(3);
                end
                y = y(:);
            else
                if x(3)^2 == 1
                    s = zeros(3,3);
                    s(1,:) = [0 x(1)/3 x(1)];
                    s(2,:) = [0 x(2)/3 x(2)];
                    s(3,:) = sqrt(s(1,:).^2 + s(2,:).^2 + 1);
                    fval = Inf;
                    y = [0;0;1];
                    for j = 1:3
                        tmp = sum((s(:,j)-x).^2);
                        if tmp < fval
                            fval = tmp;
                            y = s(:,j);
                        end
                    end    
                else
                    a(1) = 4*(x(1)^2+x(2)^2);
                    a(2) = -x(1)*a(1);
                    a(3) = 4*x(1)^2+x(1)^2*(x(1)^2+x(2)^2)-x(1)^2*x(3)^2;
                    a(4) = -4*x(1)^3;
                    a(5) = x(1)^4;
                    rt = roots(a);
                    y = [0;0;1];
                    fval = sum((y-x).^2);
                    for j = 1:length(rt)
                        if isreal(rt(j))
                            s = zeros(3,1);
                            s(1) = rt(j);
                            s(2) = x(2)/x(1)*s(1);
                            s(3) = sqrt(s(1)^2+s(2)^2+1);
                            tmp = sum((s-x).^2);
                            if tmp < fval
                                fval = tmp;
                                y = s;
                            end
                        end
                    end
                end
            end
        end
        g = sqrt(sum((x-y).^2));
    end

    % mu: D-by-m matrix
    % Z: d-by-m-by-n matrix
    % V: D-by-m-by-n matrix
    function V = coef_to_log(mu,Z)
        [frame] = orthonormal_frame(mu);
        n = size(Z,3);
        if n == 1
            D = size(mu,1);
            e1 = frame(:,:,1);
            e2 = frame(:,:,2);
            V = repmat(Z(1,:),D,1).*e1 + repmat(Z(2,:),D,1).*e2;
        else
            D = size(mu,1);
           
            e1 = frame(:,:,1);
            e2 = frame(:,:,2);
            V = repmat(Z(1,:,:),D,1,1).* repmat(e1,1,1,n) + ...
                repmat(Z(2,:,:),D,1,1).* repmat(e2,1,1,n);
        end
    end


    function [r] = vfinprod_V(U,V,~)
        n = size(U,3);
        if n == 1
            r = mean(mfd.metric(U,V));
        else
            r = mean(mfd,metric(U,V),1);
        end
    end

    function [r] = vfinprod_Z(W,Z,~)
        n = size(Z,3);
        if n == 1
            r = mean(sum(W.*Z,1));
        else
            r = squeeze(mean(sum(W.*Z,1),2));
        end
    end
end