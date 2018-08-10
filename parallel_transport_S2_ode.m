function [ode_sol,geo_sol,ptinfo] = parallel_transport_S2_ode(pp,qq,u,varargin)
% paralell transport of u from p to q along minimizing geodesic on the unit
% sphere S2.

% pp:    the coordinate of p, could be Cartesian or math spherical
% qq:    the coordinate of q, could be Cartesian or math spherical
% u:    3-by-1 vector, the tangent vector at p to be transported

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
    %addOptional(parm,'method','ode'); %ode,geo,odeandgeo
    
    parse(parm,pp,qq,u,varargin{:});
    eps = parm.Results.eps;
    verbose = parm.Results.verbose;
    odeopt = parm.Results.odeopt;

    %% utility functions
    metric = @(x,y) x(1,:).*y(1,:)+x(2,:).*y(2,:)+x(3,:).*y(3,:);
    l2norm = @(x) sqrt(metric(x,x));

    rtheta = @(theta,phi) [-sin(theta).*sin(phi); ...
              cos(theta).*sin(phi); ...
              zeros(1,length(theta))];
    rphi = @(theta,phi) [cos(theta).*cos(phi);...
             sin(theta).*cos(phi);...
             -sin(phi)];
    nml = @(theta,phi) [cos(theta).*sin(phi);...
           sin(theta).*sin(phi);...
           cos(phi)];
    Exp = @(t,p,h) cos(t)*p+sin(t)*h;
    
    ptinfo.metric = metric;
    ptinfo.norm = l2norm;
    ptinfo.r1 = rtheta;
    ptinfo.r2 = rphi;
    ptinfo.nml = nml;
    ptinfo.geodesic = Exp;

    %%
    
    [~,n] = size(pp);
    ode_sol = zeros(3,n);
    geo_sol = zeros(3,n); % geo solution

    [Cp,Sp] = co_convert_S2(pp);
    [Cq,Sq] = co_convert_S2(qq);

    theta_p = Sp(1,:);
    phi_p = Sp(2,:);
    theta_q = Sq(1,:);
    phi_q = Sq(2,:);


    rtheta_p = rtheta(theta_p,phi_p);
    rphi_p = rphi(theta_p,phi_p);
    %nml_p = nml(theta_p,phi_p);
    rtheta_q = rtheta(theta_q,phi_q);
    rphi_q = rphi(theta_q,phi_q);

    TQ = zeros(1,n);
    Ch = zeros(3,n);

    % check parallel transport
    for i = 1:n
        
        if verbose; fprintf('i=%d: ',i); end;
        
        p = Cp(:,i);  q = Cq(:,i);
        vfound = true;
        
        if abs(l2norm(p-q)) < eps
            if verbose; disp('p == q'); end
            v = u(:,i);
            v0 = v;
        elseif abs(p(3)-1) < eps || abs(q(3)-1) < eps || abs(p(3)+1) < eps || abs(q(3)+1) < eps
            if verbose; disp('p==(0,0,1) or p==(0,0,-1) or q==(0,0,1) or q==(0,0,-1)'); end
            
            if abs(p(3)-1) < eps
                theta = Sq(1,i); % phi = phi_p = phi_q
                rp_p = [cos(theta); sin(theta); 0];
                rt_p = [-sin(theta); cos(theta); 0];
                rp_q = rphi_q(:,i);
                rt_q = rtheta_q(:,i);
            else
                theta = Sp(1,i); % phi = phi_p = phi_q
                rp_q = [cos(theta); sin(theta); 0];
                rt_q = [-sin(theta); cos(theta); 0];
                rp_p = rphi_p(:,i);
                rt_p = rtheta_p(:,i);
            end
            
            
            a1 = rt_p / l2norm(rt_p);
            a2 = rp_p / l2norm(rp_p);
            
            b1 = rt_q / l2norm(rt_q);
            b2 = rp_q / l2norm(rp_q);

            u1 = dot(a1,u(:,i));
            u2 = dot(a2,u(:,i));

            if abs(dot(b1,b2)) > eps
                warning('%d: e1 is not perpendicular to e2',i);
            end
            v = u1 * b1 + u2 * b2;
            v0 = v;
            
        elseif abs(p(1)*q(2)-p(2)*q(1))<eps && p(1)*q(1)<= 0 
            if verbose; disp('geodesic between p and q passes (0,0,1) or (0,0,-1)'); end
            %phi = Sp(2,i); % phi_p == phi_q
            rp_p = rphi_p(:,i);
            rt_p = rtheta_p(:,i);
            a1 = rt_p / l2norm(rt_p);
            a2 = rp_p / l2norm(rp_p);

            rp_q = rphi_q(:,i); % [-sin(phi); cos(phi); 0];
            rt_q = rtheta_q(:,i);
            b1 = rt_q / l2norm(rt_q);
            b2 = rp_q / l2norm(rp_q);

            u1 = dot(a1,u(:,i));
            u2 = dot(a2,u(:,i));

            v = u1 * b1 + u2 * b2;
            v0 = v;

        else % the ODE case
            if verbose; disp('compute parallel transport by ODE'); end
            w = cross(p,q);
            h = cross(p,w);
            h = h / sqrt(metric(h,h));
            np = p;
            if abs(dot(h,np)) > eps || abs(dot(h,w)) > eps
                warning('h is not in tangent space');
            end

            px = p(1); hx = h(1); qx = q(1);
            %py = p(2); hy = h(2); qy = q(2);
            %pz = p(3); hz = h(3); qz = q(3);

            [k1,k2] = quadroot(px^2+hx^2,-2*qx*px,qx^2-hx^2);
            if abs(k1) > 1 && abs(k2) >0 
                warning('solution for tq');
            end

            tcand = zeros(1,4);
            if abs(k1) <= 1
                tcand(1:2) = [acos(k1) 2*pi-acos(k1)];
            end
            if abs(k2) <= 1
                tcand(3:4) = [acos(k2) 2*pi-acos(k2)];
            end
            

            bfound = false;
            for k = 1:length(tcand)
                if abs(tcand(k)) < pi && l2norm(q-Exp(tcand(k),p,h)) < eps
                    tq = tcand(k);
                    bfound = true;
                    break;
                end
            end
            
            if ~bfound
                for k = 1:length(tcand)
                    if abs(tcand(k)) < pi && l2norm(q-Exp(-tcand(k),p,h)) < eps
                        tq = -tcand(k);
                        bfound = true;
                        break;
                    end
                end
            end
        
            if ~bfound
                for j = 1:3
                    f = @(t) cos(t)*p(j)+sin(t)*h(j)-q(j);
                    tq = fzero(f,0);
                    qq = cos(tq)*p + sin(tq)*h;
                    if abs(l2norm(qq-q)) < eps
                        break;
                    end
                end
            end

            if tq < 0
                tq = -tq;
                h = -h;
            end

            while tq > 2*pi
                tq = tq - 2*pi;
            end
            
            TQ(i) = tq;
            Ch(:,i) = h;


            uu = u(:,i); 

            % geometric solution
            v0 = parallel_transport_S2_geo(p,q,uu);

            % ODE solution

            odef = @(t,y) s2odefun(t,y,p,h);


            u1 = metric(rtheta_p(:,i),uu) / l2norm(rtheta_p(:,i))^2;
            u2 = metric(rphi_p(:,i),uu) / l2norm(rphi_p(:,i))^2;

            %odef = @(t,y) s2odefun_r0(t,y,r0);

            [~,y] = ode45(odef,[0,tq],[u1;u2],odeopt);
            if any(isnan(y(:))) || any(isinf(y(:)))
                warning('NaN or Inf at i=%d',i);
                vfound = false;
            else
                v1 = y(end,1);
                v2 = y(end,2);
                v = v1*rtheta_q(:,i) + v2*rphi_q(:,i);
            end
            
            
            
            if vfound && parm.Results.check  % check solutions
                eps1 = parm.Results.check_eps;
                % are they in tangent space?
                if abs(metric(q,v0)) > eps1 || abs(metric(q,v)) > eps1
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
        if vfound
            geo_sol(:,i) = v0;
            ode_sol(:,i) = v;
        else
            warning('v is not found');
        end
    end
    
    ptinfo.geo_sol = geo_sol;
    ptinfo.ode_sol = ode_sol;
    ptinfo.TQ = TQ;
    ptinfo.h = Ch;
end
