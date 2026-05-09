 %% dati
clear
clc
% dati esterni
mu=398600;

% orbita di partenza
r_A=[-4944.4602;4989.7097;6462.6224]; v_A=[-4.7010;-4.6070;0.3576];
[a_A,e_A,i_A,OM_A,om_A,th_A]=car2kep(r_A,v_A,mu);

% orbita di arrivo
a_B=14870; e_B=0.3004; i_B=0.9583; OM_B=2.493; om_B=2.9; th_B=0.2117;
[r_B,v_B]=kep2car(a_B,e_B,i_B,OM_B,om_B,th_B,mu);

%% parametri orbitali

% vettori normali ai piani
h_A=cross(r_A,v_A);
h_B=cross(r_B,v_B);

% angolo tra i piani
alpha=acos(h_A'*h_B/(norm(h_A)*norm(h_B)));

% vettori eccentricità
e_A_v=cross(v_A,h_A)/mu-r_A/norm(r_A);
e_B_v=cross(v_B,h_B)/mu-r_B/norm(r_B);

% asse nodale relativo
N_rel=cross(h_A,h_B)/norm(cross(h_A,h_B));

% semilati retti
p_A=a_A*(1-e_A^2);
p_B=a_B*(1-e_B^2);

%     Δω
%      θm  θM
% Δω1
% Δω2
delta_w=[atan2(h_A'/norm(h_A)*cross(e_A_v,N_rel)/e_A, e_A_v'*N_rel/e_A); 
         atan2(h_B'/norm(h_B)*cross(e_B_v,N_rel)/e_B, e_B_v'*N_rel/e_B)];
delta_w=[mod(delta_w,2*pi) mod(delta_w+pi,2*pi)];

% parametri orbite di trasferimento
theta1 = @(r,p1,e1,w) mod(2*pi*(w>pi)+acos((r.^2.*e1.^2.*cos(w)+2*r.*e1.*(p1-r)+(p1-r).^2.*cos(w))./(r.^2.*e1.^2+2*r.*e1.*(p1-r).*cos(w)+(p1-r).^2))*(2*(w<pi)-1), 2*pi); % θ di manovra sulla prima orbita
theta2 = @(r,p1,e1,w) mod(theta1(r,p1,e1,w)-w,2*pi); % θ di arrivo post manovra sul'orbita di trasferimento
e2 = @(r,p1,e1,w) r.*e1.*sin(theta1(r,p1,e1,w))./(r.*e1.*sin(theta1(r,p1,e1,w))+p1.*sin(theta2(r,p1,e1,w)));
a2 = @(r,p1,e1,w) r./(1+e2(r,p1,e1,w));

%% calcoli di velocità e tempo
n=1e3;                          % numero di raggi di apogeo
r=linspace(max(a_A*(1+e_A),a_B*(1+e_B))+0.01,100e3,n);  % raggi di apogeo
delta_v=zeros(4,n,2);
delta_t=zeros(5,n,2);

dv = @(v2,v1) sqrt(sum((v2-v1).^2));

for i=1:2

    % velocità sulla prima orbita nel punto di tangenza
    [r1,v1]=kep2car(a_A,e_A,i_A,OM_A,om_A,theta1(r,p_A,e_A,delta_w(1,i)),mu);

    % velocità sulla prima orbita di trasferimento nel punto di tangenza
    [r2,v2]=kep2car(a2(r,p_A,e_A,delta_w(1,i)),e2(r,p_A,e_A,delta_w(1,i)),i_A,OM_A,mod(om_A+delta_w(1,i),2*pi),theta2(r,p_A,e_A,delta_w(1,i)),mu);

    % velocità sulla prima orbita di trasferimento all'apocentro
    [r3,v3]=kep2car(a2(r,p_A,e_A,delta_w(1,i)),e2(r,p_A,e_A,delta_w(1,i)),i_A,OM_A,mod(om_A+delta_w(1,i),2*pi),pi*ones(1,n),mu);

    % velocità sulla seconda orbita di trasferimento all'apocentro
    [r4,v4]=kep2car(a2(r,p_B,e_B,delta_w(2,i)),e2(r,p_B,e_B,delta_w(2,i)),i_B,OM_B,mod(om_B+delta_w(2,i),2*pi),pi*ones(1,n),mu);

    % velocità sulla seconda orbita di trasferimento nel punto di tangenza
    [r5,v5]=kep2car(a2(r,p_B,e_B,delta_w(2,i)),e2(r,p_B,e_B,delta_w(2,i)),i_B,OM_B,mod(om_B+delta_w(2,i),2*pi),theta2(r,p_B,e_B,delta_w(2,i)),mu);

    % velocità sull'orbita di arrivo nel punto di tangenza
    [r6,v6]=kep2car(a_B,e_B,i_B,OM_B,om_B,theta1(r,p_B,e_B,delta_w(2,i)),mu);
    
    % calcolo Δv e Δt
    delta_v(1,:,i)=dv(v2,v1);
    delta_v(2,:,i)=dv(v4,v3);
    delta_v(3,:,i)=dv(v6,v5);
    delta_v(4,:,i)=sum(delta_v(1:3,:,i));

    delta_t(1,:,i)=timeOfFlight(a_A,e_A,th_A,theta1(r,p_A,e_A,delta_w(1,i)),mu);
    delta_t(2,:,i)=timeOfFlight(a2(r,p_A,e_A,delta_w(1,i)),e2(r,p_A,e_A,delta_w(1,i)),theta2(r,p_A,e_A,delta_w(1,i)),pi,mu);
    delta_t(3,:,i)=timeOfFlight(a2(r,p_B,e_B,delta_w(2,i)),e2(r,p_B,e_B,delta_w(2,i)),pi,theta2(r,p_B,e_B,delta_w(2,i)),mu);
    delta_t(4,:,i)=timeOfFlight(a_B,e_B,theta1(r,p_B,e_B,delta_w(2,i)),th_B,mu);
    delta_t(5,:,i)=sum(delta_t(1:4,:,i));
end

%% risultati

% parametri orbite di trasferimento
kep_A=[a_A,e_A,i_A,OM_A,om_A,th_A,theta1(r(end),p_A,e_A,delta_w(1,2))];
kep_1=[a2(r(end),p_A,e_A,delta_w(1,2)),e2(r(end),p_A,e_A,delta_w(1,2)),i_A,OM_A,mod(om_A+delta_w(1,2),2*pi),theta2(r(end),p_A,e_A,delta_w(1,2)),pi];
kep_2=[a2(r(end),p_B,e_B,delta_w(2,2)),e2(r(end),p_B,e_B,delta_w(2,2)),i_B,OM_B,mod(om_B+delta_w(2,2),2*pi),pi,theta2(r(end),p_B,e_B,delta_w(2,2))];
kep_B=[a_B,e_B,i_B,OM_B,om_B,theta1(r(end),p_B,e_B,delta_w(2,2)),th_B];

close all
orbitDraw([kep_A;kep_1;kep_2;kep_B],["Ocompletion","Yes"])

figure
plot(r,delta_v(1,:,1),"r--",r,delta_v(2,:,1),"g--",r,delta_v(3,:,1),"b--",r,delta_v(4,:,1),"k--",r,delta_v(1,:,2),"r",r,delta_v(2,:,2),"g",r,delta_v(3,:,2),"b",r,delta_v(4,:,2),"k")
grid on
legend("Δv_1^A","Δv_2^A","Δv_3^A","Δv_t_o_t^A","Δv_1^B","Δv_2^B","Δv_3^B","Δv_t_o_t^B",'Location','best','fontsize',11,'NumColumns',4)
ylabel("Δv (km/s)",'fontsize',14)
xlabel("r_a (km)",'fontsize',14)
ylim([0,8])
set(gcf,'color','w');

%delta_t=delta_t./3600;
figure
plot(r,delta_t(1,:,1),"r--",r,delta_t(2,:,1),"g--",r,delta_t(3,:,1),"c--",r,delta_t(4,:,1),"b--",r,delta_t(5,:,1),"k--",r,delta_t(1,:,2),"r",r,delta_t(2,:,2),"g",r,delta_t(3,:,2),"c",r,delta_t(4,:,2),"b",r,delta_t(5,:,2),"k")
grid on
legend({'Δt_1^A','Δt_2^A','Δt_3^A','Δt_4^A','Δt_t_o_t^A','Δt_1^B','Δt_2^B','Δt_3^B','Δt_4^B','Δt_t_o_t^B'},'Location','best','fontsize',11,'NumColumns',2)
xlabel("r_a (km)",'fontsize',14)
ylabel("Δt (s)",'fontsize',14)
set(gcf,'color','w');

figure
plot(delta_v(end,:,1),delta_t(end,:,1),"k--",delta_v(end,:,2),delta_t(end,:,2),"k")
grid on
legend("A","B",'fontsize',11,'Location','best','NumColumns',2)
xlabel("Δv (km/s)",'fontsize',14)
ylabel("Δt (s)",'fontsize',14)
set(gcf,'color','w');
%delta_t=delta_t.*3600;