function [cx,sx] = sample_H2(n,varargin)
    parm = inputParser;
    parm.KeepUnmatched = true;
    addRequired(parm,'n');
    addOptional(parm,'Tlow',1);
    addOptional(parm,'Thigh',2);
    addOptional(parm,'Plow',0);
    addOptional(parm,'Phigh',2*pi);
    
    parse(parm,n,varargin{:});

    Tlow = parm.Results.Tlow;
    Thigh = parm.Results.Thigh;
    Plow = parm.Results.Plow;
    Phigh = parm.Results.Phigh;

    cx = zeros(3,n);
    sx = zeros(2,n);
    
    sx(1,:) = acosh(Tlow+(Thigh-Tlow)*rand(1,n));
    sx(2,:) = Plow + (Phigh-Plow)*rand(1,n);
    
    
    cx(1,:) = sinh(sx(1,:)).*cos(sx(2,:));
    cx(2,:) = sinh(sx(1,:)).*sin(sx(2,:));
    cx(3,:) = cosh(sx(1,:));
end