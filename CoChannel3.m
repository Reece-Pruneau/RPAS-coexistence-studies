
tic

fc=700000000; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)
bw=6000000; %bandwidth (Hz)
tilt=-3;  %mechanical downtilt (degrees)
SINR_target=10; %(dB)
nf=5; %noise factor
L_feeder=-3; %(dB)
L_body=-1; %(dB)
UE_gain=-3; %(dBi)
L_entry=-10; %(dB)
h_BSa=40; %antenna height (m)
h_BSb=40;  %antenna height (m)
d_sep=6000; % (m) seperation distance between both networks
% Note d_sep=0 coresponds to cells touching, not overlaping
r_a=4000; %cell radius (m)
n_a=100; %number of monte carlo runs
max_dBmUE=23; %max transmiting power of UEs (dBm)
min_dBmUE=-40; %min transmiting power of UEs (dBm)
nf=5; %noise figure (dB)



%count number of UE over or under tx power range
n_UEmaxA=0;
n_UEmaxB=0;
n_UEminA=0;
n_UEminB=0;

global HEXphi
global HEXtheta
global HEXd

%Noise
noisefloor=10.*log10(1.380649e-23*290*bw*1000)+nf;  %dBm
N=10.^(noisefloor./10); %mW

outline=NaN(7,2,60);

%Interfering network A

%position and antenna orrientation
BS_a=[0,0,h_BSa;0,1,sin(pi*tilt/180)];
BS_a1=[0,0,h_BSa;0,1,sin(pi*tilt/180)];
BS_a2=[0,0,h_BSa;-0.866,-0.5,sin(pi*tilt/180)];
BS_a3=[0,0,h_BSa;0.866,-0.5,sin(pi*tilt/180)];

%Center positions of each hex cell
x1=0+BS_a(1,1);
y1=r_a+BS_a(1,2);
x2=r_a.*sqrt(3)/2+BS_a(1,1);
y2=-1/2.*r_a+BS_a(1,2);
x3=-sqrt(3)/2.*r_a+BS_a(1,1);
y3=-1/2.*r_a+BS_a(1,2);


%matrix for storing drone position in all three sectors
UE_a=zeros(3,3,n_a);

%generate n_a random xy drone positions for each cell
[UE_a(1:2,1,:),outline(:,:,1)]=randHEX(r_a,y1,x1,n_a);
[UE_a(1:2,2,:),outline(:,:,2)]=randHEX(r_a,y2,x2,n_a);
[UE_a(1:2,3,:),outline(:,:,3)]=randHEX(r_a,y3,x3,n_a);

%generate n_a random drone altitudes
for i=1:n_a
UE_a(3,1,i)=(rand*270)+30;
UE_a(3,2,i)=(rand*270)+30;
UE_a(3,3,i)=(rand*270)+30;
end

%Victim Network B
BS_b=[0,4*r_a+d_sep,h_BSb;0,-1,sin(pi*tilt/180)]; %position and antenna orrientation base station
UE_b=zeros(3,57,n_a);

UE_b(3,:,:)=1.5;
for q=1
[UE_b(1:2,1,:),outline(:,:,4)]=randHEX(r_a,BS_b(1,2)+0,BS_b(1,1)+0,n_a);

[UE_b(1:2,2,:),outline(:,:,5)]=randHEX(r_a,BS_b(1,2)+1.5*r_a,BS_b(1,1)-3*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,3,:),outline(:,:,6)]=randHEX(r_a,BS_b(1,2)+1.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,n_a);
[UE_b(1:2,4,:),outline(:,:,7)]=randHEX(r_a,BS_b(1,2)+1.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,n_a);
[UE_b(1:2,5,:),outline(:,:,8)]=randHEX(r_a,BS_b(1,2)+1.5*r_a,BS_b(1,1)+3*sqrt(3)*r_a/2,n_a);

[UE_b(1:2,6,:),outline(:,:,9)]=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)-3*sqrt(3)*r_a,n_a);
[UE_b(1:2,7,:),outline(:,:,10)]=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)-2*sqrt(3)*r_a,n_a);
[UE_b(1:2,8,:),outline(:,:,11)]=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)-sqrt(3)*r_a,n_a);
[UE_b(1:2,9,:),outline(:,:,12)]=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)+0*sqrt(3)*r_a,n_a);
[UE_b(1:2,10,:),outline(:,:,13)]=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)+sqrt(3)*r_a,n_a);
[UE_b(1:2,11,:),outline(:,:,14)]=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)+2*sqrt(3)*r_a,n_a);
[UE_b(1:2,12,:),outline(:,:,15)]=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)+3*sqrt(3)*r_a,n_a);

[UE_b(1:2,13,:),outline(:,:,16)]=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)-7*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,14,:),outline(:,:,17)]=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)-5*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,15,:),outline(:,:,18)]=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)-3*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,16,:),outline(:,:,19)]=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,n_a);
[UE_b(1:2,17,:),outline(:,:,20)]=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,n_a);
[UE_b(1:2,18,:),outline(:,:,21)]=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)+3*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,19,:),outline(:,:,22)]=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)+5*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,20,:),outline(:,:,23)]=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)+7*sqrt(3)*r_a/2,n_a);

[UE_b(1:2,21,:),outline(:,:,24)]=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)-3*sqrt(3)*r_a,n_a);
[UE_b(1:2,22,:),outline(:,:,25)]=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)-2*sqrt(3)*r_a,n_a);
[UE_b(1:2,23,:),outline(:,:,26)]=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)-sqrt(3)*r_a,n_a);
[UE_b(1:2,24,:),outline(:,:,27)]=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)+0*sqrt(3)*r_a,n_a);
[UE_b(1:2,25,:),outline(:,:,28)]=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)+sqrt(3)*r_a,n_a);
[UE_b(1:2,26,:),outline(:,:,29)]=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)+2*sqrt(3)*r_a,n_a);
[UE_b(1:2,27,:),outline(:,:,30)]=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)+3*sqrt(3)*r_a,n_a);

[UE_b(1:2,28,:),outline(:,:,31)]=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)-7*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,29,:),outline(:,:,32)]=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)-5*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,30,:),outline(:,:,33)]=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)-3*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,31,:),outline(:,:,34)]=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,n_a);
[UE_b(1:2,32,:),outline(:,:,35)]=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,n_a);
[UE_b(1:2,33,:),outline(:,:,36)]=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)+3*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,34,:),outline(:,:,37)]=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)+5*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,35,:),outline(:,:,38)]=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)+7*sqrt(3)*r_a/2,n_a);

[UE_b(1:2,36,:),outline(:,:,39)]=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)-3*sqrt(3)*r_a,n_a);
[UE_b(1:2,37,:),outline(:,:,40)]=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)-2*sqrt(3)*r_a,n_a);
[UE_b(1:2,38,:),outline(:,:,41)]=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)-sqrt(3)*r_a,n_a);
[UE_b(1:2,39,:),outline(:,:,42)]=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)+0*sqrt(3)*r_a,n_a);
[UE_b(1:2,40,:),outline(:,:,43)]=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)+sqrt(3)*r_a,n_a);
[UE_b(1:2,41,:),outline(:,:,44)]=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)+2*sqrt(3)*r_a,n_a);
[UE_b(1:2,42,:),outline(:,:,45)]=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)+3*sqrt(3)*r_a,n_a);

[UE_b(1:2,43,:),outline(:,:,46)]=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)-7*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,44,:),outline(:,:,47)]=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)-5*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,45,:),outline(:,:,48)]=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)-3*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,46,:),outline(:,:,49)]=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,n_a);
[UE_b(1:2,47,:),outline(:,:,50)]=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,n_a);
[UE_b(1:2,48,:),outline(:,:,51)]=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)+3*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,49,:),outline(:,:,52)]=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)+5*sqrt(3)*r_a/2,n_a);
[UE_b(1:2,50,:),outline(:,:,53)]=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)+7*sqrt(3)*r_a/2,n_a);

[UE_b(1:2,51,:),outline(:,:,54)]=randHEX(r_a,BS_b(1,2)+12.0*r_a,BS_b(1,1)-2*sqrt(3)*r_a,n_a);
[UE_b(1:2,52,:),outline(:,:,55)]=randHEX(r_a,BS_b(1,2)+12.0*r_a,BS_b(1,1)-sqrt(3)*r_a,n_a);
[UE_b(1:2,53,:),outline(:,:,56)]=randHEX(r_a,BS_b(1,2)+12.0*r_a,BS_b(1,1)+0*sqrt(3)*r_a,n_a);
[UE_b(1:2,54,:),outline(:,:,57)]=randHEX(r_a,BS_b(1,2)+12.0*r_a,BS_b(1,1)+sqrt(3)*r_a,n_a);
[UE_b(1:2,55,:),outline(:,:,58)]=randHEX(r_a,BS_b(1,2)+12.0*r_a,BS_b(1,1)+2*sqrt(3)*r_a,n_a);

[UE_b(1:2,56,:),outline(:,:,59)]=randHEX(r_a,BS_b(1,2)+13.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,n_a);
[UE_b(1:2,57,:),outline(:,:,60)]=randHEX(r_a,BS_b(1,2)+13.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,n_a);
end


%WANTED SIGNAL
S_const=10.^((SINR_target+noisefloor)./10);


I=zeros(1,n_a); %aggregate Interference power (mW)
I_self=zeros(1,n_a);

%Calculate amount of self interference
for i=1:n_a


%Calculate amount of self interference

for w=2:length(UE_b(1,:))
    %calculate LOS displacement vector in xy plane
D=[UE_b(1,w)-BS_b(1,1),UE_b(2,w)-BS_b(1,2)];

%angle between BS max gain vector and displacement vector in xy plane
phi=180.*acos(dot(D,BS_b(2,1:2))./norm(D))./pi;
%elevation angle
theta=90+tilt-180.*atan(norm(D)./BS_b(1,3))./pi;

if rand<=0.5
           L_entry=-10;
       else
           L_entry=0;
       end

I_dBmUE=SINR_target+noisefloor-F1336(HEXphi(round(rand*999999)+1),HEXtheta(round(rand*999999)+1),18)-UE_gain-L_body-L_entry+P1546FieldStrMixed(700,50,1.5,40,0,'Suburban',HEXd(round(rand*999999)+1)./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', 1.5);

if I_dBmUE>max_dBmUE
    I_dBmUE=max_dBmUE;
    n_UEmaxB=n_UEmaxB+1;
end
if I_dBmUE<min_dBmUE
    I_dBmUE=min_dBmUE;
    n_UEminB=n_UEminB+1;
end

I_self(i)=I_self(i)+10.^((I_dBmUE+F1336(phi,theta,18)+UE_gain+L_body+L_entry-P1546FieldStrMixed(700,50,1.5,40,0,'Suburban',norm(D)./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', 1.5))./10);

end





%Drone Interference


%calculate LOS displacement vector for each drone in the 3 sectors
%repeat for n_a runs
D1=[UE_a(1,1,i)-BS_b(1,1),UE_a(2,1,i)-BS_b(1,2),UE_a(3,1,i)-BS_b(1,3)];
D2=[UE_a(1,2,i)-BS_b(1,1),UE_a(2,2,i)-BS_b(1,2),UE_a(3,2,i)-BS_b(1,3)];
D3=[UE_a(1,3,i)-BS_b(1,1),UE_a(2,3,i)-BS_b(1,2),UE_a(3,3,i)-BS_b(1,3)];

d1=[UE_a(1,1,i)-BS_a1(1,1),UE_a(2,1,i)-BS_a1(1,2),UE_a(3,1,i)-BS_a1(1,3)];
d2=[UE_a(1,2,i)-BS_a2(1,1),UE_a(2,2,i)-BS_a2(1,2),UE_a(3,2,i)-BS_a2(1,3)];
d3=[UE_a(1,3,i)-BS_a3(1,1),UE_a(2,3,i)-BS_a3(1,2),UE_a(3,3,i)-BS_a3(1,3)];


I_dBmUE1=SINR_target+noisefloor-F1336V(UE_a(:,1,i),BS_a1(1,:),BS_a1(2,:),18)-20.*log10(lambda./(4.*pi.*norm(d1)+1e-3))-UE_gain;
I_dBmUE2=SINR_target+noisefloor-F1336V(UE_a(:,2,i),BS_a2(1,:),BS_a2(2,:),18)-20.*log10(lambda./(4.*pi.*norm(d2)+1e-3))-UE_gain;
I_dBmUE3=SINR_target+noisefloor-F1336V(UE_a(:,3,i),BS_a3(1,:),BS_a3(2,:),18)-20.*log10(lambda./(4.*pi.*norm(d3)+1e-3))-UE_gain;

if I_dBmUE1>max_dBmUE
    I_dBmUE1=max_dBmUE;
    n_UEmaxA=n_UEmaxA+1;
end
if I_dBmUE1<min_dBmUE
    I_dBmUE1=min_dBmUE;
    n_UEminA=n_UEminA+1;
end
if I_dBmUE2>max_dBmUE
    I_dBmUE2=max_dBmUE;
    n_UEmaxA=n_UEmaxA+1;
end
if I_dBmUE2<min_dBmUE
    I_dBmUE2=min_dBmUE;
    n_UEminA=n_UEminA+1;
end
if I_dBmUE3>max_dBmUE
    I_dBmUE3=max_dBmUE;
    n_UEmaxA=n_UEmaxA+1;
end
if I_dBmUE3<min_dBmUE
    I_dBmUE3=min_dBmUE;
    n_UEminA=n_UEminA+1;
end

I(i)=10.^((I_dBmUE1+F1336V(UE_a(:,1,i),BS_b(1,:),BS_b(2,:),18)+20.*log10(lambda./(4.*pi.*norm(D1)))+UE_gain)./10)+10.^((I_dBmUE2+F1336V(UE_a(:,2,i),BS_b(1,:),BS_b(2,:),18)+20.*log10(lambda./(4.*pi.*norm(D2)))+UE_gain)./10)+10.^((I_dBmUE3+F1336V(UE_a(:,2,i),BS_b(1,:),BS_b(2,:),18)+20.*log10(lambda./(4.*pi.*norm(D3)))+UE_gain)./10);


end

I_0=mean(I_self);
SINR_0=S_const./(N+I_self);
SINR_const=S_const./(N+I_self+I);
bitloss=100.*(log2(1+SINR_0)-log2(1+SINR_const))./(log2(1+SINR_0)+log2(1+SINR_const));

average_throughput_loss=mean(bitloss)

figure(1)
subplot(1,2,1)
hold on
for i=1:60
    plot(outline(:,2,i),outline(:,1,i),'b')
    %plot(UE(1,i),UE(2,i),'k.')
    xlabel('x position (m)')
    ylabel('y position (m)')
end
plot(outline(:,2,4),outline(:,1,1),'r')
hold off

subplot(1,2,2)
histogram(bitloss)
title('Distribution of Bitrate loss')
xlabel('Bitrate loss (%)')
ylabel(['instances out of ', num2str(n_a)])



toc
