function varargout = changeOrbitShape(ai,ei,omi,af,ef,omf,rb)
% [Dvi, Dvf, Dt, thf, kep_tr] = changeOrbitShape(ai,ei,omi,af,ef,omf)
% executes a generalized Hohmann transfer
%
% inputs:  ai, ei: initial orbit shape parameters (1x1)
%             omi: initial orbit orientation (1xn)
%          af, ef: final orbit shape parameters (1x1, 1x1)
%             omf: final orbit orientation (1xn)
%
% outputs:    Dvi: initial impulse, signed (2xn)
%             Dvf: final impulse, signed (2xn)
%              Dt: travel time on transfer orbit (2xn)
%             thf: position on final orbit after the maneuver (2xn)
%          kep_tr: transfer orbit parameters (a,e,ω,θ1,θ2) (1xnx5)
%
% The first row of the outputs refers to the maneuver performed from the
% pericenter of the first orbit, the second row refers to the maneuver
% performed from the apocenter of the first orbit
%
% If the two apoapsi are not aligned, NaN is returned
%
% [Dv1, Dv2, Dv3, Dt, thf] = changeOrbitShape(ai,ei,omi,af,ef,omf,rb) to
% execute a bi-elliptical maneuver

if size(omi,1)>1||size(omf,1)>1
    error("Please input row vectors of angles")
end

if length(omi)==1||length(omf)==1
    cases=max(length(omi),length(omf));
    omi=omi.*ones(1,cases);
    omf=omf.*ones(1,cases);
elseif length(omi)==length(omf)
    cases=length(omi);
else
    warning("Vectors have different lengths, using the shortest");
    cases=min(length(omi),length(omf));
    omi=omi(1:min(cases,length(omi)));
    omf=omf(1:min(cases,length(omf)));
end

mu = 398600;

rai = ai*(1+ei);
rpi = ai*(1-ei);
raf = af*(1+ef);
rpf = af*(1-ef);

if nargin==6
    % primo impulso
    Dv1 = [sqrt(2*mu*(1/rpi-1/(raf+rpi)))-sqrt(2*mu*(1/rpi-1/(2*ai))); sqrt(2*mu*(1/rai-1/(rai+rpf)))-sqrt(2*mu*(1/rai-1/(2*ai)))].*(abs(omi-omf)<1e-6)+[sqrt(2*mu*(1/rpi-1/(rpf+rpi)))-sqrt(2*mu*(1/rpi-1/(2*ai))); sqrt(2*mu*(1/rai-1/(rai+raf)))-sqrt(2*mu*(1/rai-1/(2*ai)))].*(abs(abs(omi-omf)-pi)<1e-6);
    Dv1(abs(omi-omf).*[1;1]>1e-6 & abs(abs(omi-omf)-pi).*[1;1]>1e-6) = NaN;

    % secondo impulso
    Dv2 = [sqrt(2*mu*(1/raf-1/(2*af)))-sqrt(2*mu*(1/raf-1/(rpi+raf))); sqrt(2*mu*(1/rpf-1/(2*af)))-sqrt(2*mu*(1/rpf-1/(rai+rpf)))].*(abs(omi-omf)<1e-6)+[sqrt(2*mu*(1/rpf-1/(2*af)))-sqrt(2*mu*(1/rpf-1/(rpi+rpf))); sqrt(2*mu*(1/raf-1/(2*af)))-sqrt(2*mu*(1/raf-1/(rai+raf)))].*(abs(abs(omi-omf)-pi)<1e-6);
    Dv2(abs(omi-omf).*[1;1]>1e-6 & abs(abs(omi-omf)-pi).*[1;1]>1e-6) = NaN;

    % tempo di trasferimento
    Dt = [pi*sqrt((0.5*(rpi+raf))^3/mu); pi*sqrt((0.5*(rpf+rai))^3/mu)].*(abs(omi-omf)<1e-6)+[pi*sqrt((0.5*(rpi+rpf))^3/mu); pi*sqrt((0.5*(raf+rai))^3/mu)].*(abs(abs(omi-omf)-pi)<1e-6);
    Dt(abs(omi-omf).*[1;1]>1e-6 & abs(abs(omi-omf)-pi).*[1;1]>1e-6)=NaN;

    % θ a cui ci si trova sull'orbita finale dopo la manovra
    thf= [pi; 0].*(abs(omi-omf)<1e-6) + [0; pi].*(abs(abs(omi-omf)-pi)<1e-6);
    thf(abs(omi-omf).*[1;1]>1e-6 & abs(abs(omi-omf)-pi).*[1;1]>1e-6)=NaN;

    % parametri kepleriani orbita di trasferimento
    a = [(rpi+raf)/2;(rai+rpf)/2].*(abs(omi-omf)<1e-6) + [(rpi+rpf)/2;(rai+raf)/2].*(abs(abs(omi-omf)-pi)<1e-6);
    a(abs(omi-omf).*[1;1]>1e-6 & abs(abs(omi-omf)-pi).*[1;1]>1e-6)=NaN;
    e = abs([(raf-rpi)/(rpi+raf);(rpf-rai)/(rpf+rai)]).*(abs(omi-omf)<1e-6)+abs([(rpf-rpi)/(rpi+rpf);(raf-rai)/(raf+rai)]).*(abs(abs(omi-omf)-pi)<1e-6);
    e(abs(omi-omf).*[1;1]>1e-6 & abs(abs(omi-omf)-pi).*[1;1]>1e-6)=NaN;
    om = mod([omi+pi*(rpi>raf);omi+pi*(rai<rpf)],2*pi).*(abs(omi-omf)<1e-6) + mod([omi+pi*(rpi>rpf); omi+pi*(rai<raf)],2*pi).*(abs(abs(omi-omf)-pi)<1e-6);
    om(abs(omi-omf).*[1;1]>1e-6 & abs(abs(omi-omf)-pi).*[1;1]>1e-6)=NaN;

    % θ iniziale sull'orbita di trasferimento
    th_m_i = [pi*(rpi>raf); pi*(rai>rpf)].*(abs(omi-omf)<1e-6) + [pi*(rpi>rpf); pi*(rai>raf)].*(abs(abs(omi-omf)-pi)<1e-6);
    th_m_i(abs(omi-omf).*[1;1]>1e-6 & abs(abs(omi-omf)-pi).*[1;1]>1e-6)=NaN;
    
    % θ finale sull'orbita di trasferimento
    th_m_f = [pi*(rpi<raf); pi*(rai<rpf)].*(abs(omi-omf)<1e-6) + [pi*(rpi<rpf); pi*(rai<raf)].*(abs(abs(omi-omf)-pi)<1e-6);
    th_m_f(abs(omi-omf).*[1;1]>1e-6 & abs(abs(omi-omf)-pi).*[1;1]>1e-6)=NaN;

    kepout(:,:,1)=a;
    kepout(:,:,2)=e;
    kepout(:,:,3)=om;
    kepout(:,:,4)=th_m_i;
    kepout(:,:,5)=th_m_f;
    varargout{1}=Dv1;
    varargout{2}=Dv2;
    varargout{3}=Dt;
    varargout{4}=thf;
    varargout{5}=kepout;
elseif nargin==7
    if abs(omi-omf) < 1e-6
        varargout{1} = sqrt(2*mu*( (1/rpi) - (1/(rb+rpi)) )) - sqrt(2*mu*( (1/rpi) - (1/(2*ai)) ));
        varargout{2} = sqrt(2*mu*( (1/rb) - (1/(rb+rpf)) )) - sqrt(2*mu*( (1/rb) - (1/(rb+rpi)) ));
        varargout{3} = sqrt(2*mu*( (1/rpf) - (1/(2*af)) )) - sqrt(2*mu*( (1/rpf) - (1/(rb+rpf)) ));
        varargout{4} = pi*( sqrt( ((0.5*(rpi+rb))^3)/mu ) + sqrt( ((0.5*(rpf+rb))^3)/mu ) );
        varargout{5} = 0;
    else
        error('le anomalie di pericentro iniziale e finale non sono compatibili con una manovra biellittica')
    end
end