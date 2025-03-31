function G = F1336(phih,thetah,G0,tilt)

%function G=F1336(phi,theta,G0,tilt)
%
%function returns the gain of a sectoral antenna in dBi
%Based on the ITU model F.1336 recomendation 3.1
%INPUTS
%phi (degrees) horizontal plane azimuth angle
%theta (degrees) vertical plane elevation angle
%
%
%All other parameters must be manually set bellow
%
%
peak=0; %average sidelobes
%peak=1; %peak sidelobes

Ka=0.7;
Kp=0.7;
Kh=0.7;
Kv=0.3;

phi3=65;

%theta3=107.6.*10.^(-0.1.*G0); %BAD BAD BAD 
theta3=31000.*10.^(-0.1.*G0)./phi3;



%theta=asind(sind(theta)*cosd(-tilt)+cosd(theta)*cosd(phi)*sind(-tilt));
%phi=acosd(-(sind(theta)*sind(-tilt)+cosd(theta)*cosd(phi)*cosd(-tilt))/cosd(theta));
theta=asind(sind(thetah).*cosd(-tilt)+cosd(thetah).*cosd(phih).*sind(-tilt));
phi=acosd((-sind(thetah).*sind(-tilt)+cosd(thetah).*cosd(phih).*cosd(-tilt))./cosd(theta));


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