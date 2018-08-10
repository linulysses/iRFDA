function ynew = h2odefun(t,y,p,h)

    px = p(1); py = p(2); pz = p(3);
    hx = h(1); hy = h(2); hz = h(3);
    
    gm = cosh(t)*p+sinh(t)*h;
        
    % prepare ODE
    
    x = px * cosh(t) + hx * sinh(t); 
    x_prime = px * sinh(t) + hx * cosh(t);
    z = pz * cosh(t) + hz * sinh(t);
    e = py * cosh(t) + hy * sinh(t);
    
    if z^2-1 < 0
        warning('(pz * cosh(t) + hz * sinh(t))^2-1 < 0');
    end
    
    z_prime = pz * sinh(t) + hz * cosh(t);
    c = z^2-1;
    theta_prime = z_prime / sqrt(c);
    g_prime = (x_prime*c - x*z_prime*z) / ((c)^(3/2));
    
    f = 1-x^2/c; 
    
    if f < 0
        warning('1-g^2 = %f < 0', f*1e12);
        f = 0;
    end
    
    if e >= 0
        phi_prime = - g_prime / sqrt(f);
    else
        phi_prime = g_prime / sqrt(f);
    end
    
    G212 = gm(3)/sqrt(gm(1)^2+gm(2)^2);
    G122 = -gm(3)*sqrt(gm(1)^2+gm(2)^2);
    z1 = -G122*y(2)*phi_prime;
    z2 = -G212*(y(1)*phi_prime+y(2)*theta_prime);
     
    ynew = [z1;z2];
end