function [x1,x2] = quadroot(a,b,c)

    x1 = 0; x2 = 0;
    
    if a == 0
        if b == 0
            if c ~= 0
                warning('the equation is c = 0 but c~=0!');
            else
                warning('the equation is 0 = 0, there are infinite solutions');
            end
        else
            x1 = -c/b;
            x2 = x1;
        end
        
    else
        d = b^2 -4*a*c;
        if  d < 0
            warning('the roots are not real: b^2-4ac=%f',d);
        
        else
            x1 = (sqrt(d)-b)/(2*a);
            x2 = (-sqrt(d)-b)/(2*a);
            
            if x1 < 0
                x0 = x1;
                x1 = x2;
                x2 = x0;
            end
        end
    end