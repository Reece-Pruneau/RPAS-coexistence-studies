
clear all
lat1=40;
lat2=40;
lon1=-77;
lon2=-77;

fc=1900e6; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)

n=150;

hr=30;
ht=linspace(30,300,27);
d=logspace(0,5,n);


L528=zeros(30,n);


hold on
for i=1:27
    i
for t=1:n
temp=tl_p528(d(t)./1000,hr,ht(i),708,1,10);
L528(i,t)=temp.A__db;
%L452(t)=tl_p452(0.708, 50, d(1:t), h(1:t), g(1:t), z(1:t), 1.5, 30, lon1, lat1, lon2, lat2, -3, 15, 1, inf, inf, 1013, 25);

%L1546(t)=P1546FieldStrMixed(708,50,ht(i),hr,10,'Rural',d(t)./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', ht(i));
Lfs(t)=-20.*log10(lambda./(4.*pi.*sqrt((hr-ht(i)).^2+d(t).^2)));
end
plot(d,L528,'b-*')
%plot(d,L1546,'r-o')
plot(d,Lfs,'g-o')
%legend('p528','p1546','FSPL')
end
hold off
