function [dv,dt,kep_transfer]=due_impulsi(kep1,kep2,om,Rmin,mu)
    %[dv,dt,kep_transfer]=due_impulsi(kep1,kep2,om,Rmin,mu)
    %calculates the transfer time and keplerian parameters of the transfer
    %orbit between two orbits
    %
    % kep1=[a1,e1,i1,OM1,om1,th1i,th1m]
    % kep2=[a2,e2,i2,OM2,om2,th2m,th2f]
    % om: row vector of transfer orbit ω

    if(length(kep1)==6) 
        kep1(7)=kep1(6);
    end
    if(length(kep2)==6) 
        kep2(7)=kep2(6);
    end


    [r1,v1]=kep2car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(7),mu);
    [r2,v2]=kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),mu);
    
    % versore h
    h=cross(r1,r2)/(norm(cross(r1,r2))); %3x1
    h(:,:,2)=-h; %3x1x2
    % inclinazione (i)
    I=acos(h(3,:,:)); %1x1x2
    % asse nodale (N)
    N=cross([0;0;1].*ones(3,1,2),h)./vecnorm(cross([0;0;1].*ones(3,1,2),h)); %3x1x2
    % RAAN (Ω)
    OM=atan2(N(2,:,:),N(1,:,:)); %1x1x2
    % dipendenti da ω
    % versore e
    e_v=cos(om.*ones(3,1,2)).*N+sin(om.*ones(3,1,2)).*cross(h,N); %3xnx2
    % angoli θ
    theta(1,1:length(om),1:2)=atan2(pagemtimes(h,'transpose',cross(e_v,r1.*ones(3,length(om),2)),'none')/norm(r1),sum(r1.*ones(3,length(om),2).*e_v)/norm(r1)); %1xnx2
    theta(2,:,:)=atan2(pagemtimes(h,'transpose',cross(e_v,r2.*ones(3,length(om),2)),'none')/norm(r2),sum(r2.*ones(3,length(om),2).*e_v)/norm(r2)); %1xnx2
    % eccentricità (e)
    e_sc=(norm(r2)-norm(r1))./(norm(r1)*cos(theta(1,:,:))-norm(r2)*cos(theta(2,:,:))); %1xnx2
    e_sc(e_sc<0)=NaN;
    e_sc(e_sc>=1-1e-5)=NaN;

    % semiasse maggiore (a)
    a=norm(r1)*(1+e_sc.*cos(theta(1)))./(1-e_sc.^2); %1xnx2
    a(a.*(1-e_sc)<Rmin)=NaN;
    dv(1,length(om),1:2)=0;
    dt(1,length(om),1:2)=0;
    for p=1:2
        for t=1:2
            [~,v_tr]=kep2car(reshape(a(:,:,p),1,1,length(om)),reshape(e_sc(:,:,p),1,1,length(om)),I(1,1,p),OM(1,1,p),reshape(om,1,1,length(om)),reshape(theta(t,:,p),1,1,length(om)),mu);
            dv(t,1:length(om),p)=reshape(vecnorm(v_tr-v1*(t==1)-v2*(t==2)),1,length(om),1);
        end
        kep_transfer(:,:,p)=[a(:,:,p);e_sc(:,:,p);I(1,1,p).*ones(1,length(om),1);OM(1,1,p).*ones(1,length(om),1);om;theta(1,1:length(om),p);theta(2,1:length(om),p)];
    end
    dv(3,:,:)=sum(dv,1);
    dt=timeOfFlight(kep1(1),kep1(2),kep1(6),kep1(7),mu)+timeOfFlight(a,e_sc,theta(1,:,:),theta(2,:,:),mu)+timeOfFlight(kep2(1),kep2(2),kep2(6),kep2(7),mu);
%{
    % orbita non compenetra la Terra
    for c=1:length(om)
        for p=1:2
            if (a(1,c,p)*(1-e_sc(1,c,p))>=Rmin)
                % velocità sull'orbita di trasferimento
                [~,v_tr]=kep2car(a(1,c,p),e_sc(1,c,p),I(1,1,p),OM(1,1,p),om(c),theta(:,c,p)',mu);
                % Δv complessivo
                dv(p,c)=norm(v_tr(:,1)-v1)+norm(v_tr(:,2)-v2);
                % Δt complessivo
                dt(p,c)=timeOfFlight(kep1(1),kep1(2),kep1(6),kep1(7),mu)+timeOfFlight(a(1,c,p),e_sc(1,c,p),theta(1,c,p),theta(2,c,p),mu)+timeOfFlight(kep2(1),kep2(2),kep2(6),kep2(7),mu);
                kep_transfer(:,c,p)=[a(1,c,p),e_sc(1,c,p),I(1,1,p),OM(1,1,p),om(c),theta(1,c,p),theta(2,c,p)];
            else
                dv(p,c)=NaN;
                dt(p,c)=NaN;
                kep_transfer(1:7,c,p)=NaN*ones(1,7);
            end
        end
    end
%}
end