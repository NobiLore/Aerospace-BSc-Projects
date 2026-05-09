function [r,v] = kep2car(a, e, i, OM, om, th, mu)
    p=a.*(1-e.^2);%1x1xn
    r=(ones(3,1).*p./(1+e.*cos(th))).*[cos(th); sin(th); th*0];%3xnxn
  %{  
    d=zeros(1,length(om));
    for j=1:3
        d(j,:)=diag(reshape(r(j,:,:),length(om),length(om),1));
    end
    r=reshape(d,3,1,length(om));%3x1xn
%}
    v=sqrt(mu./p).*[e*0-sin(th);e+cos(th);(th+e)*0];%3xnxn
%{
    d=zeros(1,length(om));
    for j=1:3
        d(j,:)=diag(reshape(v(j,:,:),length(om),length(om),1));
    end
    v=reshape(d,3,1,length(om));%3x1xn
    %}
    R_W=[cos(OM) sin(OM) 0;
        -sin(OM) cos(OM) 0;
            0      0     1];%3x3
    R_i=[1   0       0;
         0 cos(i)  sin(i);
         0 -sin(i) cos(i)];%3x3
    R_w=[cos(om) sin(om) om*0;
        -sin(om) cos(om) om*0;
          om*0    om*0  om*0+1];%3x3xm
    r=pagemtimes(pagemtimes(pagemtimes(R_W,'transpose',R_i,'transpose'),'none',R_w,'transpose'),r);%3xnxn
    v=pagemtimes(pagemtimes(pagemtimes(R_W,'transpose',R_i,'transpose'),'none',R_w,'transpose'),v);
end