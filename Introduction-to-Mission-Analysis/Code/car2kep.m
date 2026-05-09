function [a, e, i, OM, om, th] = car2kep(r, v, mu)
    a=mu/(2*mu/norm(r)-norm(v)^2);
    h=cross(r,v);
    e_vect=cross(v,h)/mu-r/norm(r);%((norm(v)^2-mu/norm(r))*r-(r'*v)*v)/mu;
    e=norm(e_vect);
    i=acos(h(3)/norm(h));
    N=cross([0;0;1],h)+[1;0;0]*prod((cross([0;0;1],h)==[0;0;0]));
    e_vect=e_vect*(e~=0)+N*(e==0);
    OM=acos(N(1)/norm(N))*(2*(N(2)>=0)-1)+2*pi*(N(2)<0);
    if(N'*e_vect/(norm(N)*norm(e_vect))>1+1e-8)
        error("Errore 1")
    end
    om=real(acos(N'*e_vect/(norm(N)*norm(e_vect)))*(2*(e_vect(3)>=0)-1)+2*pi*(e_vect(3)<0));
    v_r=r'*v/norm(r);
    if(e_vect'*r/(norm(e_vect)*norm(r))>1+1e-8)
        error("Errore 2");
    end
    th=real(acos(e_vect'*r/(norm(e_vect)*norm(r)))*(2*(v_r>=0)-1)+2*pi*(v_r<0));
end