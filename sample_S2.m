function [cx,sx] = sample_S2(n,varargin)
parm = inputParser;
    parm.KeepUnmatched = true;
    addRequired(parm,'n');
    addOptional(parm,'Tlow',1);
    addOptional(parm,'Thigh',2);
    addOptional(parm,'Plow',0);
    addOptional(parm,'Phigh',1.5*pi);
    
    parse(parm,n,varargin{:});

    Tlow = parm.Results.Tlow;
    Thigh = parm.Results.Thigh;
    Plow = parm.Results.Plow;
    Phigh = parm.Results.Phigh;

    cx = zeros(3,n);
    sx = zeros(2,n);
    
    sx(1,:) = rand(1,n)*(Phigh-Plow)+Plow;
    sx(2,:) = (acos(1-2*rand(1,n))/pi)*(Thigh-Tlow)+Tlow;
    
    
    cx(1,:) = cos(sx(1,:)).*sin(sx(2,:));
    cx(2,:) = sin(sx(1,:)).*sin(sx(2,:));
    cx(3,:) = cos(sx(2,:));
end