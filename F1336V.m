function G = F1336V(POS1,POS2,ORR2,G0)

D=[POS1(1)-POS2(1),POS1(2)-POS2(2)];


%angle between BS max gain vector and displacement vector
phi=180.*acos(dot(D,ORR2(1:2))./norm(D))./pi;


tilt=180.*asin(ORR2(3))./pi;
theta=90+tilt-180.*atan(norm(D)./POS2(3))./pi;

Ka=0.7;
Kp=0.7;
Kh=0.7;
Kv=0.3;

phi3=65;

%theta3=107.6.*10.^(-0.1.*G0); %BAD BAD BAD 
theta3=31000.*10.^(-0.1.*G0)./phi3;

xh=abs(phi)./phi3;
xv=abs(theta)./theta3;
G180=-15+10.*log10(1+8.*Kp)-15.*log(180./theta3);

if xh<=0.5
    Ghr=-12.*xh.^2;
elseif xh>0.5
    Ghr=-12.*xh.^(2-Kh)-3.*(1-0.5^(-Kh));
end

if Ghr<G180
        Ghr=G180;
end

xk=sqrt(1-0.36*Kv);
C=10.*log10((180./theta3).^1.5*(4.^-1.5+Kv)./(1+8.*Kp))./log10(22.5./theta3);

if xv<xk
    Gvr=-12.*xv.^2;
elseif xk<=xv<4
    Gvr=-12+log10(xv.^(-1.5)+Kv);
elseif 4<=xv<(90./theta3)
    Gvr=12-C.*log10(4)-10.*log10(4.^-1.5+Kv)-C.*log10(xv);
elseif xv>=(90./theta3)
    Gvr=G180;
end

R=(Ghr+12.*(180./65).^(2-Kh)-3.*(1-0.5^(-Kh)))./(12.*(180./65).^(2-Kh)-3.*(1-0.5^(-Kh)));

G=G0+Ghr+R.*Gvr;

end
