%% dati
clear
clc
close all

mu=398600;
Rmin=6378+150;

r_A=[-4944.4602;4989.7097;6462.6224]; v_A=[-4.7010;-4.6070;0.3576];
[a_A,e_A,i_A,OM_A,om_A,th_A]=car2kep(r_A,v_A,mu);
kepA=[a_A,e_A,i_A,OM_A,om_A,th_A];

a_B=14870; e_B=0.3004; i_B=0.9583; OM_B=2.493; om_B=2.9; th_B=0.2117;
[r_B,v_B]=kep2car(a_B,e_B,i_B,OM_B,om_B,th_B,mu);
kepB=[a_B,e_B,i_B,OM_B,om_B,th_B];

n=1e4;
ver=0;
for i=1:10
    om_min(i)=2*pi*rand(1);
    range=om_min(i)+pi*(-1:2:1);
    om_min(i)=NaN;
    
    while range(2)-range(1)>1e-8
        om=linspace(range(1),range(2),n);
        [dv,dt,kep_tr]=due_impulsi(kepA,kepB,om,Rmin,mu);
        [dt_min(i),om_min(i)]=min(dt,[],'all');
        ver=floor(om_min(i)/n)+1;
        kep_tr_min(1:7,i)=kep_tr(:,mod(om_min(i),n)+om_min(i)*(om_min(i)==n),ver);
        om_min(i)=om(mod(om_min(i),n)+om_min(i)*(om_min(i)==n));
        range=om_min(i)+5*(range(2)-range(1))/(n-1)*(-1:2:1);
    end
end

clc
[dt_min,i]=min(dt_min);
om_min=mod(om_min(i),2*pi);
disp("Orbita più veloce trovata. Parametri:")
fprintf("a: %f km\ne: %f\ni: %f rad\nΩ: %f rad\nω: %f rad\nθi: %f rad\nθf: %f rad\n",kep_tr_min(1,i),kep_tr_min(2,i),kep_tr_min(3,i),kep_tr_min(4,i),kep_tr_min(5,i),mod(kep_tr(6,i),2*pi),mod(kep_tr_min(7,i),2*pi))
fprintf("Tempo di percorrenza: %f min\n", dt_min/60);
fprintf("Δv richiesti: %f, %f km/s\n",dv(1,i,ver),dv(2,i,ver));
fprintf("Δv totale richiesto: %f km/s\n",dv(3,i,ver));

%%
%close all
orbitDraw([kepA(1:5),th_A+1e-5,th_A;kep_tr_min(:,i)';kepB(1:5),th_B+1e-5,th_B;ones(2,1),zeros(2,6)],["Ocompletion","No"]);