
function G = F1336V(POS1,POS2,ORR2,G0,tilt)
%function G=F1336(phi,theta,G0)
%
%function returns the gain of a sectoral antenna in dBi
%in the direction of a specified point
%
%Based on the ITU model F.1336 recomendation 3.1
%INPUTS
% POS1 : xyz position of the area of interest (could be rx or tx antenna
% for example) 1x3 vector 
%
% POS2 : xyz position of the F1336 sectoral antenna 1x3 vector (must be
% same units as POS1
%
% ORR2 : xyz orrientation unit vector representing direction of max gain of
% the sectoral antenna 1x3 vector
%
%All other parameters must be manually set bellow
%
%
%peak=0; %average sidelobes
peak=1; %peak sidelobes

Ka=0.7;
Kp=0.7;
Kh=0.7;
Kv=0.3;

phi3=65;

%Displacement vector
d=[POS1(1)-POS2(1),POS1(2)-POS2(2),POS1(3)-POS2(3)];

%angle between BS max gain vector and displacement vector in xy plane
phih=acosd(dot(d(1:2),ORR2(1:2))./norm(d(1:2))./norm(ORR2(1:2)));

%elevation angle
thetah=acosd(dot([norm(d(1:2)),d(3)],[norm(ORR2(1:2)),ORR2(3)])./norm(d)./norm(ORR2)).*d(3)./abs(d(3));

theta=asind(sind(thetah).*cosd(-tilt)+cosd(thetah).*cosd(phih).*sind(-tilt));
phi=acosd((-sind(thetah).*sind(-tilt)+cosd(thetah).*cosd(phih).*cosd(-tilt))./cosd(theta));

%theta3=107.6.*10.^(-0.1.*G0); %BAD BAD BAD 
theta3=31000.*10.^(-0.1.*G0)./phi3;




xh=abs(phi)./phi3;
xv=abs(theta)./theta3;
    if peak==0

    G180=-15+10.*log10(1+8.*Ka)-15.*log10(180./theta3);
    xk=sqrt(1.33-0.33*Kv);
    C=10.*log10((180./theta3).^1.5*(4.^-1.5+Kv)./(1+8.*Ka))./log10(22.5./theta3);

    elseif peak==1
  
    G180=-12+10.*log10(1+8.*Kp)-15.*log10(180./theta3);
    xk=sqrt(1-0.36*Kv);
    C=10.*log10((180./theta3).^1.5*(4.^-1.5+Kv)./(1+8.*Ka))./log10(22.5./theta3);

    end

R=(Ghr(xh,Kh,G180)-Ghr(180/phi3,Kh,G180))./(-Ghr(180/phi3,Kh,G180));

G=G0+Ghr(xh,Kh,G180)+R.*Gvr(xv,xk,Kv,theta3,C,G180,peak);
    

end

function Ghr = Ghr(xh,Kh,G180)
if xh<=0.5
    Ghr=-12.*xh.^2;
elseif xh>0.5
    Ghr=-12.*xh.^(2-Kh)-3.*(1-0.5^(-Kh));
end

if Ghr<G180
        Ghr=G180;
end

end

function Gvr = Gvr(xv,xk,Kv,theta3,C,G180,peak)
if xv<xk
    Gvr=-12.*xv.^2;
elseif xk<=xv && xv<4
    Gvr=-12+10*log10(xv.^(-1.5)+Kv);
    if peak==0
        Gvr=Gvr-3;
    end
elseif 4<=xv && xv<(90/theta3)
    Gvr=-12+C.*log10(4)+10.*log10(4.^-1.5+Kv)-C.*log10(xv);
    if peak==1
        Gvr=Gvr-3;
    end
elseif xv>=(90/theta3)
    Gvr=G180;
end
end
