%% dati iniziali
clear
clc
close all
r_A=[-4944.4602;4989.7097;6462.6224];
v_A=[-4.7010;-4.6070;0.3576];
mu=398600;
Rmin=6378+150;
[a_A,e_A,i_A,OM_A,om_A,th_A]=car2kep(r_A,v_A,mu);

a_B=14870;
e_B=0.3004;
i_B=0.9583;
OM_B=2.493;
om_B=2.9;
th_B=0.2117;
[r_B,v_B]=kep2car(a_B,e_B,i_B,OM_B,om_B,th_B,mu);

Orbits=zeros(1,7);
%{
%figure
%plot3(r_A(1),r_A(2),r_A(3),"mo");
hold on
grid on
axis equal
plot3(r_B(1),r_B(2),r_B(3),"kx");
%plotOrbit([a_A,e_A,i_A,OM_A,om_A],linspace(0,2*pi,300),mu,"m",1);
%plotOrbit([a_B,e_B,i_B,OM_B,om_B],linspace(0,2*pi,300),mu,"k",0);
%}
%% cambio forma
clc
% p→a p→p
% a→p a→a

%     all  opp
% p→   1    2
% a→   3    4

om_2=[om_A, mod(om_A+pi,2*pi)];
[dv1,dv2,dt_f,th_sh,kep_sh]=changeOrbitShape(a_A,e_A,om_A,a_B,e_B,om_2);
dt_f=dt_f+timeOfFlight(a_A,e_A,th_A,[0;pi],mu).*ones(2);
fprintf("Δv cambio forma: %f km/s, %f km/s\n",dv1(2,2),dv2(2,2));
fprintf(" t cambio forma: %f s, %f s\n",timeOfFlight(a_A,e_A,th_A,pi,mu),dt_f(2,2));

% orbite dopo cambio forma

for i=1:2
    %plotOrbit([a_B,e_B,i_A,OM_A,om_2(i)],linspace(0,2*pi,200),mu,"y",0);
end
% orbite di trasferimento per cambio forma

for i=1:2
    for j=1:2
        %plotOrbit([kep_sh(i,j,1),kep_sh(i,j,2),i_A,OM_A,kep_sh(i,j,3)],linspace(kep_sh(i,j,4),kep_sh(i,j,5)+2*pi*(kep_sh(i,j,4)>kep_sh(i,j,5)),200),mu,"r",1);
    end
end

%% cambio orientazione
% calcolo ω target
h_A=cross(r_A,v_A);
h_B=cross(r_B,v_B);
e_B_v=cross(v_B,h_B)/mu-r_B/norm(r_B);
N_rel=cross(h_B,h_A);
N_rel=N_rel/norm(N_rel);
phi=acos(h_A'*h_B/(norm(h_A)*norm(h_B)));
R=[N_rel(1)^2+(1-N_rel(1)^2)*cos(phi), prod(N_rel([1,2]))*(1-cos(phi))-N_rel(3)*sin(phi), prod(N_rel([1,3]))*(1-cos(phi))+N_rel(2)*sin(phi);
   prod(N_rel([1,2]))*(1-cos(phi))+N_rel(3)*sin(phi), N_rel(2)^2+(1-N_rel(2)^2)*cos(phi), prod(N_rel([2,3]))*(1-cos(phi))-N_rel(1)*sin(phi);
   prod(N_rel([1,3]))*(1-cos(phi))-N_rel(2)*sin(phi), prod(N_rel([2,3]))*(1-cos(phi))+N_rel(1)*sin(phi), N_rel(3)^2+(1-N_rel(3)^2)*cos(phi)];
e_tr=R*e_B_v;
N_A=cross([0;0;1],h_A)/norm(cross([0;0;1],h_A));
om_per=mod(atan2(h_A'/norm(h_A)*cross(N_A,e_tr)/norm(e_tr),N_A'*e_tr/norm(e_tr)),2*pi);

%dv_per=zeros(2);
%       (p→ all) (p→ opp) (a→ all) (a→ opp)
% θ min     1        2        3        4
% θ max     5        6        7        8

% om = 1 2 1 2
%  θ = 1 2 3 4
[dv_per, th_per_1, th_per_2, dt_per]=changePeriapsisArg(a_B,e_B,[om_2 om_2],om_per,398600,[th_sh(1,:) th_sh(2,:)]);
%   dv_per = 1 2 1 2
% th_per_1 = 1 2 1 2
%            5 6 5 6
% th_per_2 = 1 2 1 2
%            5 6 5 6
%   dt_per = 1 2 3 4
%            5 6 7 8
dv_per = dv_per(1:2).*ones(2,2,2);
th_per_1(1,:,2) = th_per_1(2,:,1);
th_per_1 = th_per_1(1,1:2,:).*ones(2,2,2);
th_per_2(1,:,2) = th_per_2(2,:,1);
th_per_2 = th_per_2(1,1:2,:).*ones(2,2,2);
dt_per(1,:,2) = dt_per(2,:,1);
dt_per = [dt_per(1,1:2,:); dt_per(1,3:4,:)];

fprintf("Δv cambio periasse: %f km/s\n",dv_per(2,2,2));
fprintf(" t cambio periasse: %f s\n",dt_f(2,2)+dt_per(2,2,2));
%     θ min             θ max
%     all  opp          all  opp
% p→   1    2       p→   5    6
% a→   3    4       a→   7    8
for i=2
    for j=1:2
        for k=1:2
            %plotOrbit([a_B,e_B,i_A,OM_A,om_per],th_per_2(i,j,k),mu,"g",1);
        end
    end
end
%plotOrbit([a_B,e_B,i_A,OM_A,om_per],linspace(0,2*pi,200),mu,"c",1);

%% cambio piano

%       (θm all) (θM all) (θm opp) (θM opp)
% θ min     1        2        3        4
% θ max     5        6        7        8
[dv_pl,om_f,th3,dt_pl] = changeOrbitalPlane(a_B,e_B,i_A,OM_A,om_per,i_B,OM_B,[th_per_2(1,1,1) th_per_2(1,1,2) th_per_2(1,2,1) th_per_2(1,2,2)]);
%     θ min             θ max
%      θm  θM            θm  θM
% all   1   2       all   5   6
% opp   3   4       opp   7   8

%  dv = 1
%       5
dv_pl(1,1,2)=dv_pl(2);
dv_pl=dv_pl(1,1,:).*ones(2,2,2);
dt_pl(1,:,2) = dt_pl(2,:,1);
dt_pl = [dt_pl(1,1:2,:); dt_pl(1,3:4,:)];
fprintf("Δv cambio piano: %f km/s\n",dv_pl(2,2,2));
fprintf(" t cambio piano: %f s\n",dt_f(2,2)+dt_per(2,2,2)+dt_pl(2,2,2));
% th3 = 1
%       5
th3(1,1,2)=th3(2);
th3=th3(1,1,:).*ones(2,2,2);
for i=1:2
    %plotOrbit([a_B,e_B,i_B,OM_B,om_f],th3(1,1,i),mu,"b",1);
end

%% attesa finale
dt_wait=timeOfFlight(a_B,e_B,th3,th_B,mu);
%     θ min             θ max
%      θm  θM            θm  θM
% all   1   2       all   5   6
% opp   3   4       opp   7   8

%% organizzazione dati e combinazione

delta_v=zeros(2,2,2,2); % p/a all/opp θm/θM θmin/θMAX
delta_t=zeros(2,2,2,2); % p/a all/opp θm/θM θmin/θMAX

delta_v=(abs(dv1)+abs(dv2)).*ones(2,2,2,2); % cambio forma
delta_v=delta_v+abs(dv_per).*ones(2,2,2,2); % cambio orientazione
dv_pl_aux(1,1:2,1:2,1:2)=dv_pl;
delta_v=delta_v+dv_pl_aux.*ones(2,2,2,2); % cambio piano

delta_t=dt_f.*ones(2,2,2,2); % cambio forma
delta_t=delta_t+abs(dt_per).*ones(2,2,2,2); % cambio orientazione
dt_pl_aux(1,1:2,1:2,1:2)=dt_pl;
delta_t=delta_t+dt_pl_aux.*ones(2,2,2,2); % cambio piano

dt_pl_aux(1,1:2,1:2,1:2)=dt_wait;
delta_t=delta_t+dt_pl_aux.*ones(2,2,2,2); % attesa finale
fprintf(" t finale: %f\n",delta_t(2,2,2,2));

fprintf("Δv totale: %f km/s\n",delta_v(2,2,2,2));
fprintf("Δt totale: %f s\n",delta_t(2,2,2,2));

figure
for i=1:2
    for j=1:2
        for k=1:2
            plot(delta_v(1,i,j,k),delta_t(1,i,j,k),"x","Color",[(3*i-2)/5,j-1,k-1]);
            hold on
            plot(delta_v(2,i,j,k),delta_t(2,i,j,k),"*","Color",[(3*i-2)/5,j-1,k-1]);
        end
    end
end
grid on
legend("p,a,m,m","a,a,m,m","p,a,m,M","a,a,m,M", "p,a,M,m","a,a,M,m","p,a,M,M","a,a,M,M" , "p,o,m,m","a,o,m,m","p,o,m,M","a,o,m,M", "p,o,M,m","a,o,M,m","p,o,M,M","a,o,M,M")
xlabel("Δv")
ylabel("Δt")

%% tempo minimo
%{
figure
plot3(r_A(1),r_A(2),r_A(3),"ko");
hold on
grid on
axis equal
%plotOrbit([a_A,e_A,i_A,OM_A,om_A],linspace(th_A,2*pi,300),mu,"k",0);

%plotOrbit([kep_sh(1,2,1),kep_sh(1,2,2),i_A,OM_A,kep_sh(1,2,3)],linspace(kep_sh(1,2,4),kep_sh(1,2,5)+2*pi*(kep_sh(1,2,4)>kep_sh(1,2,5)),200),mu,"b",1);

%plotOrbit([a_B,e_B,i_A,OM_A,om_2(2)],linspace(th_sh(1,2),th_per_1(1,2,1),200),mu,"c",0);

%plotOrbit([a_B,e_B,i_A,OM_A,om_per],th_per_2(1,2,1),mu,"g",1);

%plotOrbit([a_B,e_B,i_A,OM_A,om_per],linspace(th_per_2(1,2,1),th3(2,1,1)+2*pi,200),mu,"y",1);

%plotOrbit([a_B,e_B,i_B,OM_B,om_f],th3(2,1,1),mu,"r",1);

%plotOrbit([a_B,e_B,i_B,OM_B,om_B],linspace(th3(1,1,1),th_B+2*pi,300),mu,"m",0);
plot3(r_B(1),r_B(2),r_B(3),"mx");
%}
%% velocità minima
%figure
%plot3(r_A(1),r_A(2),r_A(3),"ko");
hold on
grid on
axis equal
kep_A=[a_A,e_A,i_A,OM_A,om_A,th_A,pi];
kep_2=[kep_sh(2,2,1),kep_sh(2,2,2),i_A,OM_A,kep_sh(2,2,3),kep_sh(2,2,4),kep_sh(2,2,5)];
kep_3=[a_B,e_B,i_A,OM_A,om_2(2),th_sh(2,2),th_per_1(2,2,2)];
kep_4=[a_B,e_B,i_A,OM_A,om_per,th_per_2(2,2,2),th3(2,2,2)];
kep_B=[a_B,e_B,i_B,OM_B,om_B,th3(2,2,2),th_B];
orbitDraw([kep_A;kep_2;kep_3;kep_4;kep_B])