function [geo_sol,ode_sol,ptinfo] = parallel_transport_H2(pp,qq,u,varargin)

    parm = inputParser;
    parm.KeepUnmatched = true;
    addRequired(parm,'pp');
    addRequired(parm,'qq');
    addRequired(parm,'u');
    addOptional(parm,'eps',1e-12);
    addOptional(parm,'check',true);
    addOptional(parm,'verbose',false);
    odeopt = odeset('RelTol',1e-5,'AbsTol',1e-5);
    addOptional(parm,'odeopt',odeopt);
    addOptional(parm,'check_eps',1e-2);
 
    
    if nargout >= 2
        ode_sol_on = true;
    else
        ode_sol_on = false;
    end
    
    parse(parm,pp,qq,u,varargin{:});
    eps = parm.Results.eps;
    verbose = parm.Results.verbose;
    odeopt = parm.Results.odeopt;
      

    %% utility functions
    lorentz = @(x,y) x(1,:).*y(1,:)+x(2,:).*y(2,:)-x(3,:).*y(3,:);
    lorentz_norm = @(x) sqrt(lorentz(x,x));
    l2norm = @(x) sqrt(dot(x,x));

    rtheta = @(theta,phi) [cosh(theta).*cos(phi); ...
              cosh(theta).*sin(phi); ...
              sinh(theta)];
    rphi = @(theta,phi) [-sinh(theta).*sin(phi);...
             sinh(theta).*cos(phi);...
             zeros(1,length(theta))];
    l2nml = @(theta,phi) [-sinh(theta).*cos(phi);...
           -sinh(theta).*sin(phi);...
           cosh(theta)];
    Exp = @(t,p,h) cosh(t)*p+sinh(t)*h;
    
    ptinfo.metric = lorentz;
    ptinfo.norm = lorentz_norm;
    ptinfo.r1 = rtheta;
    ptinfo.r2 = rphi;
    ptinfo.nml = @(theta,phi) [sinh(theta).*cos(phi); ...
                               sinh(theta).*sin(phi); ...
                               cosh(theta)];
    ptinfo.geodesic = Exp;
    
    
    %% prepare data
    [~,n] = size(pp);
    ode_sol = zeros(3,n);
    geo_sol = zeros(3,n); % geo solution

    [Cp,Sp] = co_convert_H2(pp);
    [Cq,Sq] = co_convert_H2(qq);
    theta_p = Sp(1,:);
    phi_p = Sp(2,:);
    theta_q = Sq(1,:);
    phi_q = Sq(2,:);

    rtheta_p = rtheta(theta_p,phi_p);
    rphi_p = rphi(theta_p,phi_p);
    l2nml_p = l2nml(theta_p,phi_p);
    rtheta_q = rtheta(theta_q,phi_q);
    rphi_q = rphi(theta_q,phi_q);
    %nml_q = nml(theta_q,phi_q);

   
    % parallel transport
    

    Ch = zeros(3,n);
    TQ = zeros(1,n);
    
    for i = 1:n
        
        if verbose; fprintf('i=%d: ',i); end;
        
        p = Cp(:,i);  q = Cq(:,i);
        if l2norm(p-q) < eps
            if verbose;  disp('CASE: p == q'); end;
            v = u(:,i);
            v0 = v;

        elseif p(3) == 1
            if verbose;  disp('p ==  (0,0,1)'); end;
            phi = Sq(2,i); % phi = phi_p = phi_q
            rp_p = [-sin(phi); cos(phi); 0];
            rt_p = [cos(phi); sin(phi);0]; %rtheta_p(:,i);
            a1 = rt_p / lorentz_norm(rt_p);
            a2 = rp_p / lorentz_norm(rp_p);

            rp_q = rphi_q(:,i); % [-sin(phi); cos(phi); 0];
            rt_q = rtheta_q(:,i);
            b1 = rt_q / lorentz_norm(rt_q);
            b2 = rp_q / lorentz_norm(rp_q);

            u1 = lorentz(a1,u(:,i));
            u2 = lorentz(a2,u(:,i));

            if abs(lorentz(b1,b2)) > eps
                warning('%d: e1(q) is not perpendicular to e2(q) in geo-method',i);
            end
            v = u1 * b1 + u2 * b2;
            v0 = v;
            
        elseif q(3) == 1
            if verbose;  disp('q ==  (0,0,1)'); end;
            phi = Sp(2,i); % phi = phi_p = phi_q
            rp_q = [-sin(phi); cos(phi); 0];
            rt_q = [cos(phi); sin(phi); 0];
            b1 = rt_q / lorentz_norm(rt_q);
            b2 = rp_q / lorentz_norm(rp_q);

            rp_p = rphi_p(:,i); % [-sin(phi); cos(phi); 0];
            rt_p = rtheta_p(:,i);
            a1 = rt_p / lorentz_norm(rt_p);
            a2 = rp_p / lorentz_norm(rp_p);

            u1 = lorentz(a1,u(:,i));
            u2 = lorentz(a2,u(:,i));

            if abs(lorentz(a1,a2)) > eps
                warning('%d: e1(p) is not perpendicular to e2(p) in geo-method',i);
            end
            v = u1 * b1 + u2 * b2;
            v0 = v;

        elseif abs(p(1)*q(2)-p(2)*q(1))<eps && p(1)*q(1)<= 0 
            if verbose;  disp('geodesic between p and q passes (0,0,1)'); end;
            %phi = Sp(2,i); % phi_p == phi_q
            rp_p = rphi_p(:,i);
            rt_p = rtheta_p(:,i);
            a1 = rt_p / lorentz_norm(rt_p);
            a2 = rp_p / lorentz_norm(rp_p);

            rp_q = rphi_q(:,i); % [-sin(phi); cos(phi); 0];
            rt_q = rtheta_q(:,i);
            b1 = rt_q / lorentz_norm(rt_q);
            b2 = rp_q / lorentz_norm(rp_q);

            u1 = lorentz(a1,u(:,i));
            u2 = lorentz(a2,u(:,i));

            v = u1 * b1 + u2 * b2;
            v0 = v;

        else % general case
            
            % find h and tq such that q = cosh(t)*p + sinh(t)*h
            
            w = cross(p,q);
            h = cross(l2nml_p(:,i),w);
            h = h / sqrt(lorentz(h,h));
            np = p;
            if abs(lorentz(h,np)) > eps || abs(dot(h,w)) > eps
                warning('h is not in tangent space');
            end

            px = p(1); hx = h(1); qx = q(1); 
            [k1,k2] = quadroot(px+hx,-2*qx,px-hx);
            if k2 < 0 && k1 < 0
                warning('the solution is negative');
            end
        
            if k2 > 0;   t2 = log(k2);  else  t2 = 0;  end
            if k1 > 0;   t1 = log(k1);  else  t1 = 0;  end
        
            if l2norm(Exp(t1,p,h)-q) > eps
                if l2norm(Exp(t2,p,h)-q) > eps
                    warning('no solution for tq');
                end
            end
            
            r1 = l2norm(Exp(t1,p,h)-q);
            r2 = l2norm(Exp(t2,p,h)-q);
            if r1 < r2
                tq = t1;
            else
                tq = t2;
            end
        
            if tq < 0
                tq = -tq;
                h = -h;
            end
            
            Ch(:,i) = h;
            TQ(1,i) = tq;
            
            % find parallel transport

            py = p(2); pz = p(3);
            %qy = q(2); qz = q(3);
            hx = h(1); hy = h(2); hz = h(3);

            uu = u(:,i);

            % geometric solution
            h1 = [hx;hy;-hz];
            p1 = [px;py;-pz];
            gp = cross(h1,p1);
            gp = gp / sqrt(lorentz(gp,gp));

            hq = sinh(tq)*p+cosh(tq)*h; % parallel transport of h to q
            gq = cross([hq(1);hq(2);-hq(3)],[q(1);q(2);-q(3)]);
            gq = gq / sqrt(lorentz(gq,gq));

            v0 = lorentz(uu,gp)*gq + lorentz(uu,h)*hq;

            if ode_sol_on
                % ODE solution
                odef = @(t,y) h2odefun(t,y,p,h);

                u1 = lorentz(rtheta_p(:,i),uu)/lorentz(rtheta_p(:,i),rtheta_p(:,i));
                u2 = lorentz(rphi_p(:,i),uu)/lorentz(rphi_p(:,i),rphi_p(:,i));

                [~,y] = ode45(odef,[0 tq],[u1;u2],odeopt);
                if any(isnan(y(:))) || any(isinf(y(:)))
                    warning('NaN or Inf for v for ode-method');
                    v = zeros(3,1);
                else
                    v1 = y(end,1);
                    v2 = y(end,2);
                    v = v1*rtheta_q(:,i) + v2*rphi_q(:,i);
                end
            end
        end
        
        geo_sol(:,i) = v0;
        
        if ode_sol_on;  ode_sol(:,i) = v; end;
       
        if ode_sol_on
            if parm.Results.check  % check solutions
                eps1 = parm.Results.check_eps;
                % are they in tangent space?
                if abs(lorentz(q,v0)) > eps1 || abs(lorentz(q,v)) > eps1
                    warning('i=%d: v0 or v not in tangent space',i);
                    disp([p q h]);
                    disp([uu v v0]);
                end

                if abs(l2norm(v0-v)) > eps1
                    warning('i=%d: v0 is not equal to v',i);
                    disp([p q h]);
                    disp([uu v v0]);
                end
            end
        end
        
    end
    
    ptinfo.geo_sol = geo_sol;
    ptinfo.ode_sol = ode_sol;
    ptinfo.TQ = TQ;
    ptinfo.h = Ch;
    
end