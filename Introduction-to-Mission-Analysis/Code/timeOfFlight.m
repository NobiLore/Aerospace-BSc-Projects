function Dt = timeOfFlight(a, e, thi, thf, mu)
    %TIME OF FLIGHT Calculates the time required to go from an initial angle to
    % a final angle
    %
    %   input: a, e: shape parameters
    %           thi: initial angle(s)
    %           thf: final angle(s)
    %            mu: planetary constant
    %
    %   output:  Dt: time required from thi to thf
    
    if nargin == 4
        mu = 398600.433;
    end
%{
    if size(thi,1)>1||size(thf,1)>1
        error("Please input row vectors of angles")
    end
%}
    %{
    if prod(prod(prod(size(thi)==size(thf))))==1
        rc=size(thi);
    elseif prod(size(thi)==[1 1])||prod(size(thf)==[1 1])
        rc=max(size(thi),size(thf));
    else
        warning("Angle matrices have different dimensions, using the smallest");
        rc=min(size(thi),size(thf));
        thi=thi(1:rc(1),1:rc(2));
        thf=thf(1:rc(1),1:rc(2));
    end
    %}
    Ei = 2*atan(sqrt((1-e)./(1+e)).*tan(thi/2));%.*ones(rc(1),rc(2));
    Ef = 2*atan(sqrt((1-e)./(1+e)).*tan(thf/2));%.*ones(rc(1),rc(2));
    Dt = sqrt(a.^3/mu).*(mod(Ef - Ei,2*pi) - e.*(sin(Ef) - sin(Ei)));
%{
    Ei = 2*atan(sqrt((1-e)./(1+e)).*tan(thi/2));%.*ones(rc(1),rc(2));
    Ef = 2*atan(sqrt((1-e)./(1+e)).*tan(thf/2));%.*ones(rc(1),rc(2));
    Ei-Ef
    Ef=Ef+2*pi*(Ei-Ef>1e-7);
    Dt = sqrt(a.^3/mu).*(Ef - Ei - e.*(sin(Ef) - sin(Ei)));
%}
end