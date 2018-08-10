function z = s2odefun(t,y,p,h)

    px = p(1); pz = p(3);
    hx = h(1); hz = h(3);
    
    gm = cos(t)*p+sin(t)*h;
        
    % prepare ODE
    
    eps = 1e-10;

    a = 1-gm(3)^2;
    phi_prime = -(hz*cos(t)-pz*sin(t))/sqrt(a);
    g = gm(1)/sqrt(a);
    theta_gamma = acos(g);
    if gm(2) < 0
        theta_gamma = 2*pi - theta_gamma;
    end
    
    x_prime = -sin(t)*px + cos(t)*hx;
    z_prime = -sin(t)*pz + cos(t)*hz;
    g_prime = (x_prime * (a) + gm(1)*gm(3)*z_prime)/(a^(3/2));
    
    f = 1-g^2;
   
    if f < eps
        theta_prime = 0;
    else
        if theta_gamma < pi
            theta_prime = -g_prime / sqrt(f);
        else
            theta_prime = g_prime / sqrt(f);
        end
    end
    
    cot_phi = gm(3) / sqrt(a);
    sin_cos_phi = gm(3)*sqrt(a);
    
    z1 = -cot_phi * (y(1) * phi_prime + y(2) * theta_prime);
    z2 = sin_cos_phi * y(1) * theta_prime;
    
    z = [z1;z2];
end