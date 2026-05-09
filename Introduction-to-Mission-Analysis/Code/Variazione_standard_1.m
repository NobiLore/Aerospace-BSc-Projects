%% data
clear
clc
close all

mu=398600;
Rmin=6378+150;

r_A=[-4944.4602;4989.7097;6462.6224]; v_A=[-4.7010;-4.6070;0.3576];
[a_A,e_A,i_A,OM_A,om_A,th_A]=car2kep(r_A,v_A,mu);
kep_A=[a_A,e_A,i_A,OM_A,om_A,th_A,NaN]; %orbita iniziale

a_B=14870; e_B=0.3004; i_B=0.9583; OM_B=2.493; om_B=2.9; th_B=0.2117;
[r_B,v_B]=kep2car(a_B,e_B,i_B,OM_B,om_B,th_B,mu);
kep_B=[a_B,e_B,i_B,OM_B,om_B,NaN,th_B];

dt_wait=timeOfFlight(a_A,e_A,th_A,pi,mu);

kep_A(7)=pi;
[dv_sh_1,dv_sh_2,dt_sh,th_sh,kep_sh]=changeOrbitShape(a_A,e_A,om_A,a_B,0,om_A);
dv_sh_1=dv_sh_1(2);
dv_sh_2=dv_sh_2(2);
dt_sh=dt_sh(2);
th_sh=th_sh(2);
kep_sh=reshape(kep_sh(2,:,:),1,5);
kep_2=[kep_sh(1:2),i_A,OM_A,kep_sh(3:5)]; %orbita per cambio forma

[dv_pl,om_pl,th_pl,dt_pl]=changeOrbitalPlane(a_B,0,i_A,OM_A,om_A,i_B,OM_B,th_sh);
[dt_pl,i]=min(dt_pl);
dv_pl=dv_pl(i);
th_pl=th_pl(i);
kep_3=[a_B,0,i_A,OM_A,om_A,0,th_pl]; %orbita dopo cambio forma prima di cambio piano
kep_4=[a_B,0,i_B,OM_B,om_pl,th_pl,-om_pl+om_B]; %orbita dopo cambio piano prima di cambio forma

th_pl_2=th_pl+om_pl-om_B;

dt_wait_2=timeOfFlight(a_B,0,th_pl_2,0,mu);

[dv_sh_1_,dv_sh_2_,dt_sh_,th_sh_,kep_sh_]=changeOrbitShape(a_B,0,om_B,a_B,e_B,om_B);
dv_sh_1_=dv_sh_1_(1);
dv_sh_2_=dv_sh_2_(1);
dt_sh_=dt_sh_(1);
th_sh_=th_sh_(1);
kep_sh_=reshape(kep_sh_(1,:,:),1,5);
kep_5=[kep_sh_(1:2),i_B,OM_B,kep_sh_(3:5)]; %orbita per cambio forma
kep_B(6)=pi;

dt_wait_3=timeOfFlight(a_B,e_B,th_sh_,th_B,mu);

dv=abs(dv_sh_1)+abs(dv_sh_2)+dv_pl+abs(dv_sh_1_)+abs(dv_sh_2_);
dt=dt_wait+dt_sh+dt_pl+dt_wait_2+dt_sh_+dt_wait_3;
%%
close all
orbitDraw([kep_A;kep_2;kep_3;kep_4;kep_5;kep_B],["Ocompletions","yes"])

xlabel("km")
ylabel("km")
zlabel("km")