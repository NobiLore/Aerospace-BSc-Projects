%% dati
clear
clc
close all
r_A=[-4944.4602;4989.7097;6462.6224];
v_A=[-4.7010;-4.6070;0.3576];
mu=398600;
Rmin=6378+150;
[a_A,e_A,i_A,OM_A,om_A,th_A]=car2kep(r_A,v_A,mu);
kep_A=[a_A,e_A,i_A,OM_A,om_A,th_A,NaN];

a_B=14870;
e_B=0.3004;
i_B=0.9583;
OM_B=2.493;
om_B=2.9;
th_B=0.2117;
[r_B,v_B]=kep2car(a_B,e_B,i_B,OM_B,om_B,th_B,mu);
kep_B=[a_B,e_B,i_B,OM_B,om_B,NaN,th_B];

%% calcoli
n_1=1;
n_2=1;
n_3=1;
%th1=[th_A linspace(0,2*pi*(1-1/(n_1-1)),n_1-1)];
th1=linspace(4.16,4.17,n_1);
th1=4.1606;
%th1=th_A;
%th2=[th_B linspace(0,2*pi*(1-1/(n_2-1)),n_2-1)];
%th2=linspace(0,2*pi*(1-1/(n_2)),n_2);
th2=linspace(4.2,4.25,n_2);
th2=4.2212;
%th2=th_B;
om=linspace(0,2*pi*(1-1/(n_3)),n_3);
om=linspace(4.9,5,n_3);
om=4.9749;
V_T=NaN*ones(n_1,n_2,2);

[r1,v1]=kep2car(a_A,e_A,i_A,OM_A,om_A,th1,mu);
[r2,v2]=kep2car(a_B,e_B,i_B,OM_B,om_B,th2,mu);

h=zeros(3,1);
N=zeros(3,1);
K=[0;0;1];
e_v=zeros(3,1);
theta=zeros(1,2);
e_sc=0;
I=0;
a=0;
OM=0;
v_tr=zeros(3,2);
deltaV=0;

deltaV_min=1e5;
Vmin=zeros(11,1); % θ1; θ2; a; e; i; Ω; ω; θ'; θ"; Δv; Δt
delta_t_min=1e5;
tmin=zeros(11,1); % θ1; θ2; a; e; i; Ω; ω; θ'; θ"; Δv; Δt

wait1=waitbar(0,'0');
wait2=waitbar(0,'0');
wait3=waitbar(0,'0');
for i=1:n_1
    waitbar(i/n_1,wait1,i);
    for j=1:n_2
        waitbar(j/n_2,wait2,j);
        v_min_ij=NaN;
        t_min_ij=NaN;
        % se piano orbitale è valido
        if(abs(acos(r1(:,i)'*r2(:,j)/(norm(r1(:,i)*norm(r2(:,j)))))-pi/2)<pi/2-1e-10)
            % versore h
            h=cross(r1(:,i),r2(:,j))/(norm(cross(r1(:,i),r2(:,j))));
            for ver=1:2
                h=-h;
                % inclinazione (i)
                I=acos(h(3));
                % asse nodale (N)
                N=cross(K,h)/norm(cross(K,h));
                % RAAN (Ω)
                OM=atan2(N(2),N(1));
                % dipendenti da ω
                for k=1:n_3
                    waitbar(k/n_3,wait3,k);
                    % versore e
                    e_v=cos(om(k))*N+sin(om(k))*cross(h,N);
                    % angoli θ
                    theta(1)=atan2(h'*cross(e_v,r1(:,i))/norm(r1(:,i)),r1(:,i)'*e_v/norm(r1(:,i)));
                    theta(2)=atan2(h'*cross(e_v,r2(:,j))/norm(r2(:,j)),r2(:,j)'*e_v/norm(r2(:,j)));
                    % eccentricità (e)
                    e_sc=(norm(r2(:,j))-norm(r1(:,i)))/(norm(r1(:,i))*cos(theta(1))-norm(r2(:,j))*cos(theta(2)));
                    
                    % se l'orbita di trasferimento è valida
                    if(e_sc>=0 && e_sc<1)
                        % semiasse maggiore (a)
                        a=norm(r1(:,i))*(1+e_sc*cos(theta(1)))/(1-e_sc^2);
                        % orbita non compenetra la Terra
                        if (a*(1-e_sc)>=Rmin)
                            % velocità sull'orbita di trasferimento
                            [r,v_tr]=kep2car(a,e_sc,I,OM,om(k),theta,mu);
                            % Δv complessivo
                            deltaV=norm(v_tr(:,1)-v1(:,i))+norm(v_tr(:,2)-v2(:,j));
                            % Δt complessivo
                            delta_t=timeOfFlight(a_A,e_A,th_A,th1(i),mu);
                            delta_t=delta_t+timeOfFlight(a,e_sc,theta(1),theta(2),mu);
                            delta_t=delta_t+timeOfFlight(a_B,e_B,th2(j),th_B,mu);

                            %plot(delta_t,deltaV,"k.");
                            
                            if(deltaV<deltaV_min)
                                deltaV_min=deltaV;
                                Vmin=[th1(i);th2(j);a;e_sc;I;OM;om(k);theta(1);theta(2);deltaV;delta_t];
                            end
                            %{
                            if(delta_t<delta_t_min)
                                delta_t_min=delta_t;
                                tmin=[th1(i);th2(j);a;e_sc;I;OM;om(k);theta(1);theta(2);deltaV;delta_t];
                            end
                            %}
                            v_min_ij=min(deltaV,v_min_ij);
                            %t_min_ij=min(delta_t,t_min_ij);
                        end
                    end
                end
            end
        end
        V_T(i,j,:)=[v_min_ij, t_min_ij];
    end
end
clearvars -except a_A e_A i_A OM_A om_A th_A r_A v_A a_B e_B i_B OM_B om_B th_B r_B v_B delta_t_min deltaV_min mu tmin Vmin th1 th2 V_T n_1 n_2 n_3 Rmin kep_A kep_B
%% risultati
%clc
format shorteng
fprintf("Δv minimo:\n");
fprintf("Δv: %f km/s\n",Vmin(10));
fprintf("Δt: %fs \n",Vmin(11));
fprintf("θ1: %f rad \nθ2: %f rad \nω: %f rad \n",Vmin(1),Vmin(2),Vmin(7));
fprintf("Eccentricità trasferimento: %f \n",Vmin(4));
fprintf("Semiasse maggiore trasferimento: %f km\n",Vmin(3));

fprintf("\n");
fprintf("Δt minimo:\n");
fprintf("Δt: %fs \n",tmin(11));
fprintf("Δv: %f km/s\n",tmin(10));
fprintf("θ1: %f° \nθ2: %f° \nω: %f° \n",rad2deg(tmin(1)),rad2deg(tmin(2)),rad2deg(tmin(7)));
fprintf("Eccentricità trasferimento: %f \n",tmin(4));
fprintf("Semiasse maggiore trasferimento: %f km\n",tmin(3));

%% plot
close all
kep_2=Vmin(3:9)';
kep_A(7)=Vmin(1);
kep_B(6)=Vmin(2);
orbitDraw([kep_A;kep_2;kep_B]);
xlabel("km");
ylabel("km");
zlabel("km");
%{
figure

plotOrbit([a_A,e_A,i_A,OM_A,om_A],linspace(0,2*pi,300),mu,"r",0);
hold on
grid on
axis equal
plot3(r_A(1),r_A(2),r_A(3),"ro");
plotOrbit([a_B,e_B,i_B,OM_B,om_B],linspace(0,2*pi,300),mu,"b",0);
plot3(r_B(1),r_B(2),r_B(3),"bx");
%}
%% dfg
% plotOrbit(Vmin(3:7,1)',linspace(Vmin(8),2*pi*(Vmin(8)>Vmin(9))+Vmin(9),100),mu,"c",1);
% plotOrbit(tmin(3:7,1)',linspace(tmin(8),2*pi*(tmin(8)>tmin(9))+tmin(9),100),mu,"g",1);
%{
for i=28:29
    [m,j]=min(V_T(i,end,:,:,end),[],4);
    [m,k]=min(m,[],3);
    plotOrbit(V_T(i,end,k,j(k),3:7),linspace(V_T(i,end,k,j(k),8),V_T(i,end,k,j(k),9)+2*pi*(V_T(i,end,k,j(k),8)>V_T(i,end,k,j(k),9)),200),mu,[i/n_2,0,1-i/n_2],0);
end
%}
xlim([-3e4,3e4])
ylim([-3e4, 3e4])
zlim([-3e4, 3e4])
%surf(th1,th2,-min(V_T(:,:,:,1),[],3),"EdgeColor","none");
%%
%{
figure
surf(th1(2:end),th2(2:end),-V_T(2:end,2:end,1),"EdgeColor","none");
title("delta V");
xlabel("θ1");
ylabel("θ2");

figure
surf(th1(2:end),th2(2:end),-V_T(2:end,2:end,2),"EdgeColor","none");
title("delta t");
xlabel("θ1");
ylabel("θ2");
%}

%zlim([0,1e5])
%{
figure
wait=waitbar(0,'0');
for i=1:n_1
    waitbar(i/n_1,wait,i);
    for j=1:n_2
        for k=1:n_3
            plot(V_T(i,j,k,2),V_T(i,j,k,1),".","MarkerEdgeColor",[i/n_1,j/n_2,k/n_3]);
        end
    end
end
%}
%%
%{
figure
plot(th1,V_T(:,338,1),"r")
hold on
grid on
ylabel("Δv")
plot(th2,V_T(332,:,1),"g")
legend("θ1","θ2")
%%
figure
plot(1,1)
hold off
pause(0.5)
wait1=waitbar(0,'0');
for i=330:340
    waitbar(i/n_2,wait1,i)
    plot(th1(2:end),V_T(2:end,i,1));
    hold on
    ylim([5,20])
    pause(2)
end
legend("330","331","332","333","334","335","336","337","338","339","340")
%}