%%% STD STRATEGY %%%
%  -attesa fino al punto di manovra piu` conveniente
%  -manovra di cambio piano
% -attesa fino al punto di manovra
% -cambio di anomalia del pericentro
% -attesa fino al punto di manovra
% -trasferimento pericentro-apocentro
% -attesa fino al punto di arrivo
%   assegnato

clear
clc
close all

mu=398600.433;

% orbita di partenza
r_A=[-4944.4602;4989.7097;6462.6224]; v_A=[-4.7010;-4.6070;0.3576];
[a_A,e_A,i_A,OM_A,om_A,th_A]=car2kep(r_A,v_A,mu);
kep_0=[a_A,e_A,i_A,OM_A,om_A,th_A,th_A];
kep_A=[a_A,e_A,i_A,OM_A,om_A,th_A,NaN];

% orbita di arrivo
a_B=14870; e_B=0.3004; i_B=0.9583; OM_B=2.493; om_B=2.9; th_B=0.2117;
[r_B,v_B]=kep2car(a_B,e_B,i_B,OM_B,om_B,th_B,mu);
kep_B=[a_B,e_B,i_B,OM_B,om_B,pi,th_B];

% cambio piano
[dv_pl,om_2,th_pl,dt_pl]=changeOrbitalPlane(a_A,e_A,i_A,OM_A,om_A,i_B,OM_B,th_A);
[dv_pl,i]=min(dv_pl);
th_pl=th_pl(i);
dt_pl=dt_pl(i);
kep_A(7)=th_pl;
kep_2=[a_A,e_A,i_B,OM_B,om_2,th_pl,NaN];

% cambio periasse
[dv_pa,th_pa_i,th_pa_f,dt_pa]=changePeriapsisArg(a_A,e_A,om_2,om_B,mu,th_pl);
dt_pa=dt_pa(2);
th_pa_i=th_pa_i(2);
th_pa_f=th_pa_f(2);
kep_2(7)=th_pa_i;
kep_3=[a_A,e_A,i_B,OM_B,om_B,th_pa_f,0];

% attesa per cambio forma
dt_wait_1=timeOfFlight(a_A,e_A,th_pa_f,0,mu);

% cambio forma
[dv_sh_1,dv_sh_2,dt_sh,th_f_sh,kep_4]=changeOrbitShape(a_A,e_A,om_B,a_B,e_B,om_B);
dv_sh_1=dv_sh_1(1);
dv_sh_2=dv_sh_2(1);
dt_sh=dt_sh(1);
th_f_sh=th_f_sh(1);
kep_B(6)=th_f_sh;
kep_4=[kep_4(1,1,1),kep_4(1,1,2),i_B,OM_B,kep_4(1,1,3),kep_4(1,1,4),kep_4(1,1,5)];

% attesa fino al punto d'arrivo
dt_wait_2=timeOfFlight(a_B,e_B,th_f_sh,th_B,mu);

%% risultati
close all
fprintf("Δv totale richiesto dalla manovra: %f km/s\n",abs(dv_pl)+abs(dv_pa)+abs(dv_sh_1)+abs(dv_sh_2));
fprintf("Tempo richiesto dalla manovra: %f h\n",(dt_pl+dt_pa+dt_sh(1)+dt_wait_2+dt_wait_1));
kep=[kep_0;kep_A;kep_2;kep_3;kep_4;kep_B];
orbitDraw(kep);
xlabel("km");
ylabel("km");
zlabel("km"); 