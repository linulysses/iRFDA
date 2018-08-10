function [Cmu,Smu,fval] = frechet_mean_H2(points,varargin)

% points:  3-by-n

% parameterization
%      canonical:  (sinh(theta)*cos(phi), sinh(theta)*sin(phi), cosh(theta)
%      R2       :  (sinh(u)*cosh(v), sinh(v), cosh(u)*cosh(v))
parm = inputParser;
parm.KeepUnmatched = true;
addRequired(parm,'points');
addOptional(parm,'parameterization','R2'); % or R2 or canonical
addOptional(parm,'random_start',10);
addOptional(parm,'verbose',false);

parse(parm,points,varargin{:});

sumsquare = @(x) sum(x.^2);

warning('off','all');

[Cp,Sp] = co_convert_H2(points);

if strcmp(parm.Results.parameterization,'canonical')
    
    % traditional parametrization
    
    options = optimoptions('fmincon','Display','off');

    f = @(x) sumsquare(acosh(-lorentz([sinh(x(1))*cos(x(2)); ...
        sinh(x(1))*sin(x(2)); ...
        cosh(x(1))],Cp)));
    x0 =  mean(Sp,2);
    [mu,fval] = fmincon(f,x0,[],[],[],[],[0,0],[Inf,2*pi],[],options);
    
    theta_min = min(Sp(1,:));
    theta_max = max(Sp(2,:));
    
    for k = 1:parm.Results.random_start
        x0 = [theta_min+rand(1,1)*(theta_max-theta_min); 2*rand(1,1)*pi];
        [mu_tmp,fval_tmp] = fmincon(f,x0,[],[],[],[],[0,0],[Inf,2*pi],[],options);
        if fval_tmp < fval
            fval = fval_tmp;
            mu = mu_tmp;
        end
    end
    
    Smu = mu;
    Cmu = co_convert_H2(mu);
else
    
    Tp = C_to_T(Cp);
    x0 = mean(Tp,2);
    
    grad = false;
    
    if ~grad
        
        f = @(x) sumsquare(acosh(-lorentz([sinh(x(1))*cosh(x(2)); ...
            sinh(x(2)); ...
            cosh(x(1))*cosh(x(2))],Cp)));
        options = optimoptions('fminunc','Algorithm','trust-region','Display','off');
        [mu,fval] = fminunc(f,x0,options);
    else
        f = @(x) objfun(x(1),x(2),Cp);
        options = optimoptions('fminunc','Algorithm','trust-region',...
            'Display','off','GradObj','on');
        [mu,fval] = fminunc(f,x0,options);
    end
    
    u_max = max(Tp(1,:));
    u_min = min(Tp(1,:));
    v_max = max(Tp(2,:));
    v_min = min(Tp(2,:));
    
    for k = 1:parm.Results.random_start
        x0 = [u_min+rand(1,1)*(u_max-u_min);
            v_min+rand(1,1)*(v_max-v_min)];
        [mu_tmp,fval_tmp] = fminunc(f,x0,options);
        if fval_tmp < fval
            fval = fval_tmp;
            mu = mu_tmp;
        end
    end
    
    [~,mu] = co_convert_H2(T_to_C(mu));
    Smu = mu;
    Cmu = co_convert_H2(mu);
end

warning('on','all')

%% helper functions
    function [f,g] = objfun(u,v,Cp)
        
        a = -lorentz([sinh(u)*cosh(v); ...
            sinh(v); ...
            cosh(u)*cosh(v)],Cp);
        
        f = sumsquare(acosh(a));
        
        au = -lorentz(Cp,[cosh(u)*cosh(v); ...
            0; ...
            sinh(u)*cosh(v)]);
        av = -lorentz(Cp,[sinh(u)*sinh(v); ...
            cosh(v); ...
            cosh(u)*sinh(v);]);
        
        pu = sumsquare(2*acosh(a)./sqrt(a.^2-1).*au);
        pv = sumsquare(2*acosh(a)./sqrt(a.^2-1).*av);
        g = [pu;pv];
    end

    function [Tp] = C_to_T(Cp)
        v = asinh(Cp(2,:)); 
        u = asinh(Cp(1,:)./sqrt(Cp(2,:).^2+1));
        Tp = [u;v];
    end

    function [Cp] = T_to_C(Tp)
        u = Tp(1,:);
        v = Tp(2,:);
        Cp = [sinh(u).*cosh(v); ...
            sinh(v); ...
            cosh(u).*cosh(v)];
    end


end


