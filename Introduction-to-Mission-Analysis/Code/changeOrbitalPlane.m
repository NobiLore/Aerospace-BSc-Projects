function [Dv,om2,th_m,Dt] = changeOrbitalPlane(a,e,i1,OM1,om1,i2,OM2,th0)
% [Dv,om2,th_m,Dt] = changeOrbitalPlane(a,e,i1,OM1,om1,i2,OM2,th0)
%
% Executes a change of plane transfer
%
% inputs:  a, e: shape parameters (1x1)
%       i1, OM1: initial plane parameters (1x1, 1x1 or 1xn)
%       i2, OM2: final plane parameters (1x1, 1x1 or 1xn)
%           om1: initial periapsis angle (1x1 or 1xn)
%           th0: initial angle (required only if Dt are needed) (1x1 or 1xn)
%
% outputs:   Dv: change of plane impulse (2xn)
%           om2: final periapsis angle (1xn)
%          th_m: transfer points (2xn)
%            Dt: time required from th0 to th_m (2xn)
%
% The first row of the outputs refers to the first intersection after the
% pericenter, the second row refers to the opposite point
%
% Note: Earth's planetary constant is used

mu = 398600.433; %for geocentric orbits
alpha = acos(cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(OM2-OM1)); % 1xn

u1 = atan2(sin(OM2-OM1)/sin(alpha)*sin(i2),(cos(alpha)*cos(i1)-cos(i2))/(sin(alpha)*sin(i1))); % 1xn
u2 = atan2(sin(OM2-OM1)/sin(alpha)*sin(i1),(cos(i1)-cos(alpha)*cos(i2))/(sin(alpha)*sin(i2))); % 1xn

th_m=mod(u1.*(2*(OM2>OM1)-1) - om1,2*pi); % 1xn
om2=mod(om1+(u2-u1).*(2*(OM2>OM1)-1),2*pi); % 1xn

th_m=[mod(th_m,pi); mod(th_m,pi)+pi]; % 2xn

Dv = 2*(sqrt(mu/(a*(1-e^2)))*(1+e*cos(th_m))).*sin(alpha/2); %2xn

if nargin == 7 && nargout > 3
    warning('Cannot calculate time without starting position, NaN is returned')
    Dt = NaN;
    return
else
    Dt = timeOfFlight(a,e,th0.*[1;1],th_m.*ones(size(th0)),mu); % 2xn
end

end