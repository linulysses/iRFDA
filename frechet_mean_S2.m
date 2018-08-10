function [Cmu,Smu] = frechet_mean_S2(points)

    [Cp,Sp] = co_convert_S2(points);

    sumsquare = @(x) sum(x.^2);

    f = @(x) sumsquare(acos(cos(x(2)) * cos(Sp(2,:)) ...
                     +sin(x(2)) * sin(Sp(2,:)) .* cos(Sp(1,:)-x(1))));
    y0 = mean(Cp,2);
    y0 = y0 / norm(y0);
    [~,x0] = co_convert_S2(y0);

    options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
    [Smu,fval] = fmincon(f,x0,[],[],[],[],[0,0],[2*pi,pi],[],options);
    Cmu = co_convert_S2(Smu);  
    
    
    if abs(Cmu(3)-1) < 1e-8 || abs(Cmu(3)+1) < 1e-8
        y0 = mean(Cp,2);
        y0 = y0 / norm(y0);
        [~,x0] = co_convert_S2(y0);
        [z,fv] = fmincon(f,x0,[],[],[],[],[0,0],[2*pi,pi],[],options);
        if fv < fval
            Smu = z;
            Cmu = co_convert_S2(Smu); 
        end
    end
    
    iter = 10;
    n = size(points,2);
    while (iter > 0) && (abs(Cmu(3)-1) < 1e-8 || abs(Cmu(3)+1) < 1e-8)
        iter = iter - 1;
        x0 = Sp(:,randi([1,n],1)); % random start
        [z,fv] = fmincon(f,x0,[],[],[],[],[0,0],[2*pi,pi],[],options);
        if fv < fval
            Smu = z;
            Cmu = co_convert_S2(Smu); 
        end
    end
    
end