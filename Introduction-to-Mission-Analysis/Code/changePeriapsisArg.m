function [Dv, th1, thf, Dt] = changePeriapsisArg(a, e, omi, omf, mu, th0)
% [Dv, th1, thf, Dt] = changePeriapsisArg(a, e, omi, omf, mu, th0)
%CHANGE PERIAPSIS ARGUMENT Executes a change of periapsis transfer
%
%   input:   a, e: shape parameters
%             omi: initial periapsis angle(s) (1x1 or 1xn)
%             omf: final periapsis angle(s) (1x1 or 1xn)
%             th0: initial angle(s) (required only if Dt is needed) (1xn)
%
%   output:    Dv: change of periapsis arg. impulse (1xn)
%             th1: transfer points on the initial orbit (2xn)
%             thf: transfer points on the final orbit (2xn)
%              Dt: time required from th0 to th1 (2xn)
%
% The first row of the output refers to the first maneuver point from the
% pericenter, the second row refers to the opposite point
    if size(omi,1)>1||size(omf,1)>1
        error("Please input row vectors of angles")
    end
    
    if length(omi)==1||length(omf)==1
        cases=max(length(omi),length(omf));
    elseif length(omi)==length(omf)
        cases=length(omi);
    else
        warning("Vectors have different lengths, using the shortest");
        cases=min(length(omi),length(omf));
        omi=omi(1:min(cases,length(omi)));
        omf=omf(1:min(cases,length(omf)));
    end
    omi=mod(omi.*ones(1,cases),2*pi); %1xn
    omf=mod(omf.*ones(1,cases),2*pi); %1xn
    
    if nargin == 4
        mu = 398600.433; %for geocentric orbits
    end
    
    Dom = omf - omi; %1xn
    
    th1 = mod([Dom/2; pi+Dom/2],2*pi);%2xn
    thf = mod([-Dom/2; pi-Dom/2],2*pi);%2xn
    
    Dv = 2*sqrt(mu/(a*(1-e^2)))*e*sin(abs(Dom)./2);%1xn
    
    if nargin <= 5 && nargout > 3
        warning('Cannot calculate time without starting position, NaN is returned')
        Dt = NaN;
        return
    elseif nargin <= 5 && nargout <= 3
        return
    else
        Dt = timeOfFlight(a,e,th0.*[1;1],th1,mu);%2xn
    end

end