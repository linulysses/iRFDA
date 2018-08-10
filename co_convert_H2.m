function [cp,sp] = co_convert_H2(p)
% p: the matrix of coordinate, either 3-by-n or 2-by-m or 1-by-2 or 1-by-3

    [m,n] = size(p);
    if m == 1 || n == 1
        p = p(:);
        n = 1;
    elseif m < 2 || m > 3
        error('p must be 3-by-n or 2-by-n');
    end
        
        sp = zeros(2,n);  % Sphereical coordinate of p
        
        if m == 3  % the p is given in Cartesian coordiante
            cp = p;
            if any(cp(3)<1)
                error('p is not in H2: z<1');
            end
            
            sp(1,:) = acosh(cp(3,:)); %theta
            
            Q = (cp(3)>1); % the index of points that are not (0,0,1)
            sp(2,~Q) = 0; % special case for those (0,0,1)
            sinh_theta = sqrt(cp(3,:).^2-1);
            sgn = sign(cp(2,:));
            Q1 = (sgn == 0) & Q;
            if any(Q1)
                sp(2,Q1) = acos(cp(1,Q1)./sinh_theta(Q1));
            end
            
            Q2 = (sgn ~= 0) & Q;
            if any(Q2)
                sp(2,Q2) = acos(sgn(Q2).*cp(1,Q2)./sinh_theta(Q2)) + (1-sgn(Q2))*pi/2;
            end

        else  % p is given the math spherical coordinate
            sp = p;
            cp = [sinh(sp(1,:)).*cos(sp(2,:)); ...
                  sinh(sp(1,:)).*sin(sp(2,:)); ...
                  cosh(sp(1,:))];
        end
    end