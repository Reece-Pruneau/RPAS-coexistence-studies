
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
h_BSa=20; %antenna height (m)
h_BSb=20;  %antenna height (m)
d_sep=6000; % (m) seperation distance between both networks
% Note d_sep=0 coresponds to cells touching, not overlaping
r_a=4000; %cell radius (m)
n_a=100; %number of monte carlo runs
max_dBmUE=23; %max transmiting power of UEs (dBm)
min_dBmUE=-40; %min transmiting power of UEs (dBm)



%count number of UE over or under tx power range
n_UEmaxA=0;
n_UEmaxB=0;
n_UEminA=0;
n_UEminB=0;

global HEXphi
global HEXtheta
global HEXd

%Noise
noisefloor=-98.977; %dBm
N=0.6.*10.^(noisefloor./10); %mW

outline=NaN(7,2,57+3);
%Interfering network A

%position and antenna orrientation
BS_a=[0,0,h_BSa;0,1,sin(pi*tilt/180)];
BS_a1=[0,0,h_BSa;0,1,sin(pi*tilt/180)];
BS_a2=[0,0,h_BSa;-0.866,-0.5,sin(pi*tilt/180)];
BS_a3=[0,0,h_BSa;0.866,-0.5,sin(pi*tilt/180)];

%Center positions of each hex cell
x1=0+BS_a(1,1);
y1=r_a+BS_a(2,1);
x2=r_a.*sqrt(3)/2+BS_a(1,1);
y2=-1/2.*r_a+BS_a(2,1);
x3=-sqrt(3)/2.*r_a+BS_a(1,1);
y3=-1/2.*r_a+BS_a(2,1);

%Victim Network B
BS_b=[0,4*r_a+d_sep,h_BSb;0,-1,sin(pi*tilt/180)]; %position and antenna orrientation base station


%PLOT NEWORKS
%{
dd=7*r_a;
UE=NaN(2,60);
for q=1
[UE(:,1),outline(:,:,1)]=randHEX(r_a,dd+0,BS_b(1,1)+0,1);

[UE(:,2),outline(:,:,2)]=randHEX(r_a,dd+1.5*r_a,BS_b(1,1)-3*sqrt(3)*r_a/2,1);
[UE(:,3),outline(:,:,3)]=randHEX(r_a,dd+1.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,1);
[UE(:,4),outline(:,:,4)]=randHEX(r_a,dd+1.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,1);
[UE(:,5),outline(:,:,5)]=randHEX(r_a,dd+1.5*r_a,BS_b(1,1)+3*sqrt(3)*r_a/2,1);

[UE(:,6),outline(:,:,6)]=randHEX(r_a,dd+3.0*r_a,BS_b(1,1)-3*sqrt(3)*r_a,1);
[UE(:,7),outline(:,:,7)]=randHEX(r_a,dd+3.0*r_a,BS_b(1,1)-2*sqrt(3)*r_a,1);
[UE(:,8),outline(:,:,8)]=randHEX(r_a,dd+3.0*r_a,BS_b(1,1)-sqrt(3)*r_a,1);
[UE(:,9),outline(:,:,9)]=randHEX(r_a,dd+3.0*r_a,BS_b(1,1)+0*sqrt(3)*r_a,1);
[UE(:,10),outline(:,:,10)]=randHEX(r_a,dd+3.0*r_a,BS_b(1,1)+sqrt(3)*r_a,1);
[UE(:,11),outline(:,:,11)]=randHEX(r_a,dd+3.0*r_a,BS_b(1,1)+2*sqrt(3)*r_a,1);
[UE(:,12),outline(:,:,12)]=randHEX(r_a,dd+3.0*r_a,BS_b(1,1)+3*sqrt(3)*r_a,1);

[UE(:,13),outline(:,:,13)]=randHEX(r_a,dd+4.5*r_a,BS_b(1,1)-7*sqrt(3)*r_a/2,1);
[UE(:,14),outline(:,:,14)]=randHEX(r_a,dd+4.5*r_a,BS_b(1,1)-5*sqrt(3)*r_a/2,1);
[UE(:,15),outline(:,:,15)]=randHEX(r_a,dd+4.5*r_a,BS_b(1,1)-3*sqrt(3)*r_a/2,1);
[UE(:,16),outline(:,:,16)]=randHEX(r_a,dd+4.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,1);
[UE(:,17),outline(:,:,17)]=randHEX(r_a,dd+4.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,1);
[UE(:,18),outline(:,:,18)]=randHEX(r_a,dd+4.5*r_a,BS_b(1,1)+3*sqrt(3)*r_a/2,1);
[UE(:,19),outline(:,:,19)]=randHEX(r_a,dd+4.5*r_a,BS_b(1,1)+5*sqrt(3)*r_a/2,1);
[UE(:,20),outline(:,:,20)]=randHEX(r_a,dd+4.5*r_a,BS_b(1,1)+7*sqrt(3)*r_a/2,1);

[UE(:,21),outline(:,:,21)]=randHEX(r_a,dd+6.0*r_a,BS_b(1,1)-3*sqrt(3)*r_a,1);
[UE(:,22),outline(:,:,22)]=randHEX(r_a,dd+6.0*r_a,BS_b(1,1)-2*sqrt(3)*r_a,1);
[UE(:,23),outline(:,:,23)]=randHEX(r_a,dd+6.0*r_a,BS_b(1,1)-sqrt(3)*r_a,1);
[UE(:,24),outline(:,:,24)]=randHEX(r_a,dd+6.0*r_a,BS_b(1,1)+0*sqrt(3)*r_a,1);
[UE(:,25),outline(:,:,25)]=randHEX(r_a,dd+6.0*r_a,BS_b(1,1)+sqrt(3)*r_a,1);
[UE(:,26),outline(:,:,26)]=randHEX(r_a,dd+6.0*r_a,BS_b(1,1)+2*sqrt(3)*r_a,1);
[UE(:,27),outline(:,:,27)]=randHEX(r_a,dd+6.0*r_a,BS_b(1,1)+3*sqrt(3)*r_a,1);

[UE(:,28),outline(:,:,28)]=randHEX(r_a,dd+7.5*r_a,BS_b(1,1)-7*sqrt(3)*r_a/2,1);
[UE(:,29),outline(:,:,29)]=randHEX(r_a,dd+7.5*r_a,BS_b(1,1)-5*sqrt(3)*r_a/2,1);
[UE(:,30),outline(:,:,30)]=randHEX(r_a,dd+7.5*r_a,BS_b(1,1)-3*sqrt(3)*r_a/2,1);
[UE(:,31),outline(:,:,31)]=randHEX(r_a,dd+7.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,1);
[UE(:,32),outline(:,:,32)]=randHEX(r_a,dd+7.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,1);
[UE(:,33),outline(:,:,33)]=randHEX(r_a,dd+7.5*r_a,BS_b(1,1)+3*sqrt(3)*r_a/2,1);
[UE(:,34),outline(:,:,34)]=randHEX(r_a,dd+7.5*r_a,BS_b(1,1)+5*sqrt(3)*r_a/2,1);
[UE(:,35),outline(:,:,35)]=randHEX(r_a,dd+7.5*r_a,BS_b(1,1)+7*sqrt(3)*r_a/2,1);

[UE(:,36),outline(:,:,36)]=randHEX(r_a,dd+9.0*r_a,BS_b(1,1)-3*sqrt(3)*r_a,1);
[UE(:,37),outline(:,:,37)]=randHEX(r_a,dd+9.0*r_a,BS_b(1,1)-2*sqrt(3)*r_a,1);
[UE(:,38),outline(:,:,38)]=randHEX(r_a,dd+9.0*r_a,BS_b(1,1)-sqrt(3)*r_a,1);
[UE(:,39),outline(:,:,39)]=randHEX(r_a,dd+9.0*r_a,BS_b(1,1)+0*sqrt(3)*r_a,1);
[UE(:,40),outline(:,:,40)]=randHEX(r_a,dd+9.0*r_a,BS_b(1,1)+sqrt(3)*r_a,1);
[UE(:,41),outline(:,:,41)]=randHEX(r_a,dd+9.0*r_a,BS_b(1,1)+2*sqrt(3)*r_a,1);
[UE(:,42),outline(:,:,42)]=randHEX(r_a,dd+9.0*r_a,BS_b(1,1)+3*sqrt(3)*r_a,1);

[UE(:,43),outline(:,:,43)]=randHEX(r_a,dd+10.5*r_a,BS_b(1,1)-7*sqrt(3)*r_a/2,1);
[UE(:,44),outline(:,:,44)]=randHEX(r_a,dd+10.5*r_a,BS_b(1,1)-5*sqrt(3)*r_a/2,1);
[UE(:,45),outline(:,:,45)]=randHEX(r_a,dd+10.5*r_a,BS_b(1,1)-3*sqrt(3)*r_a/2,1);
[UE(:,46),outline(:,:,46)]=randHEX(r_a,dd+10.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,1);
[UE(:,47),outline(:,:,47)]=randHEX(r_a,dd+10.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,1);
[UE(:,48),outline(:,:,48)]=randHEX(r_a,dd+10.5*r_a,BS_b(1,1)+3*sqrt(3)*r_a/2,1);
[UE(:,49),outline(:,:,49)]=randHEX(r_a,dd+10.5*r_a,BS_b(1,1)+5*sqrt(3)*r_a/2,1);
[UE(:,50),outline(:,:,50)]=randHEX(r_a,dd+10.5*r_a,BS_b(1,1)+7*sqrt(3)*r_a/2,1);

[UE(:,51),outline(:,:,51)]=randHEX(r_a,dd+12.0*r_a,BS_b(1,1)-2*sqrt(3)*r_a,1);
[UE(:,52),outline(:,:,52)]=randHEX(r_a,dd+12.0*r_a,BS_b(1,1)-sqrt(3)*r_a,1);
[UE(:,53),outline(:,:,53)]=randHEX(r_a,dd+12.0*r_a,BS_b(1,1)+0*sqrt(3)*r_a,1);
[UE(:,54),outline(:,:,54)]=randHEX(r_a,dd+12.0*r_a,BS_b(1,1)+sqrt(3)*r_a,1);
[UE(:,55),outline(:,:,55)]=randHEX(r_a,dd+12.0*r_a,BS_b(1,1)+2*sqrt(3)*r_a,1);

[UE(:,56),outline(:,:,56)]=randHEX(r_a,dd+13.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,1);
[UE(:,57),outline(:,:,57)]=randHEX(r_a,dd+13.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,1);

[UE(:,58),outline(:,:,58)]=randHEX(r_a,y1,x1,1);
[UE(:,59),outline(:,:,59)]=randHEX(r_a,y2,x2,1);
[UE(:,60),outline(:,:,60)]=randHEX(r_a,y3,x3,1);
end
hold on
for i=2:60
    plot(outline(:,2,i),outline(:,1,i),'b')
    plot(UE(1,i),UE(2,i),'k.')
    xlabel('x position (m)')
    ylabel('y position (m)')
end
plot(outline(:,2,1),outline(:,1,1),'r')
hold off
%}

S_dBmUE=20; %transmiting power of UEs (dBm)

I_self=zeros(1,n_a);

%Self interfering cells
for w=1:n_a
for q=1
UE_b(1:2,1)=randHEX(r_a,BS_b(1,2)+0,BS_b(1,1)+0,1);

UE_b(1:2,2)=randHEX(r_a,BS_b(1,2)+1.5*r_a,BS_b(1,1)-3*sqrt(3)*r_a/2,1);
UE_b(1:2,3)=randHEX(r_a,BS_b(1,2)+1.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,1);
UE_b(1:2,4)=randHEX(r_a,BS_b(1,2)+1.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,1);
UE_b(1:2,5)=randHEX(r_a,BS_b(1,2)+1.5*r_a,BS_b(1,1)+3*sqrt(3)*r_a/2,1);

UE_b(1:2,6)=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)-3*sqrt(3)*r_a,1);
UE_b(1:2,7)=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)-2*sqrt(3)*r_a,1);
UE_b(1:2,8)=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)-sqrt(3)*r_a,1);
UE_b(1:2,9)=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)+0*sqrt(3)*r_a,1);
UE_b(1:2,10)=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)+sqrt(3)*r_a,1);
UE_b(1:2,11)=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)+2*sqrt(3)*r_a,1);
UE_b(1:2,12)=randHEX(r_a,BS_b(1,2)+3.0*r_a,BS_b(1,1)+3*sqrt(3)*r_a,1);

UE_b(1:2,13)=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)-7*sqrt(3)*r_a/2,1);
UE_b(1:2,14)=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)-5*sqrt(3)*r_a/2,1);
UE_b(1:2,15)=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)-3*sqrt(3)*r_a/2,1);
UE_b(1:2,16)=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,1);
UE_b(1:2,17)=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,1);
UE_b(1:2,18)=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)+3*sqrt(3)*r_a/2,1);
UE_b(1:2,19)=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)+5*sqrt(3)*r_a/2,1);
UE_b(1:2,20)=randHEX(r_a,BS_b(1,2)+4.5*r_a,BS_b(1,1)+7*sqrt(3)*r_a/2,1);

UE_b(1:2,21)=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)-3*sqrt(3)*r_a,1);
UE_b(1:2,22)=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)-2*sqrt(3)*r_a,1);
UE_b(1:2,23)=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)-sqrt(3)*r_a,1);
UE_b(1:2,24)=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)+0*sqrt(3)*r_a,1);
UE_b(1:2,25)=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)+sqrt(3)*r_a,1);
UE_b(1:2,26)=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)+2*sqrt(3)*r_a,1);
UE_b(1:2,27)=randHEX(r_a,BS_b(1,2)+6.0*r_a,BS_b(1,1)+3*sqrt(3)*r_a,1);

UE_b(1:2,28)=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)-7*sqrt(3)*r_a/2,1);
UE_b(1:2,29)=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)-5*sqrt(3)*r_a/2,1);
UE_b(1:2,30)=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)-3*sqrt(3)*r_a/2,1);
UE_b(1:2,31)=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,1);
UE_b(1:2,32)=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,1);
UE_b(1:2,33)=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)+3*sqrt(3)*r_a/2,1);
UE_b(1:2,34)=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)+5*sqrt(3)*r_a/2,1);
UE_b(1:2,35)=randHEX(r_a,BS_b(1,2)+7.5*r_a,BS_b(1,1)+7*sqrt(3)*r_a/2,1);

UE_b(1:2,36)=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)-3*sqrt(3)*r_a,1);
UE_b(1:2,37)=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)-2*sqrt(3)*r_a,1);
UE_b(1:2,38)=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)-sqrt(3)*r_a,1);
UE_b(1:2,39)=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)+0*sqrt(3)*r_a,1);
UE_b(1:2,40)=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)+sqrt(3)*r_a,1);
UE_b(1:2,41)=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)+2*sqrt(3)*r_a,1);
UE_b(1:2,42)=randHEX(r_a,BS_b(1,2)+9.0*r_a,BS_b(1,1)+3*sqrt(3)*r_a,1);

UE_b(1:2,43)=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)-7*sqrt(3)*r_a/2,1);
UE_b(1:2,44)=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)-5*sqrt(3)*r_a/2,1);
UE_b(1:2,45)=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)-3*sqrt(3)*r_a/2,1);
UE_b(1:2,46)=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,1);
UE_b(1:2,47)=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,1);
UE_b(1:2,48)=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)+3*sqrt(3)*r_a/2,1);
UE_b(1:2,49)=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)+5*sqrt(3)*r_a/2,1);
UE_b(1:2,50)=randHEX(r_a,BS_b(1,2)+10.5*r_a,BS_b(1,1)+7*sqrt(3)*r_a/2,1);

UE_b(1:2,51)=randHEX(r_a,BS_b(1,2)+12.0*r_a,BS_b(1,1)-2*sqrt(3)*r_a,1);
UE_b(1:2,52)=randHEX(r_a,BS_b(1,2)+12.0*r_a,BS_b(1,1)-sqrt(3)*r_a,1);
UE_b(1:2,53)=randHEX(r_a,BS_b(1,2)+12.0*r_a,BS_b(1,1)+0*sqrt(3)*r_a,1);
UE_b(1:2,54)=randHEX(r_a,BS_b(1,2)+12.0*r_a,BS_b(1,1)+sqrt(3)*r_a,1);
UE_b(1:2,55)=randHEX(r_a,BS_b(1,2)+12.0*r_a,BS_b(1,1)+2*sqrt(3)*r_a,1);

UE_b(1:2,56)=randHEX(r_a,BS_b(1,2)+13.5*r_a,BS_b(1,1)-sqrt(3)*r_a/2,1);
UE_b(1:2,57)=randHEX(r_a,BS_b(1,2)+13.5*r_a,BS_b(1,1)+sqrt(3)*r_a/2,1);
end

%Calculate amount of self interference

for i=2:length(UE_b(1,:))
    %calculate LOS displacement vector in xy plane
D=[UE_b(1,i)-BS_b(1,1),UE_b(2,i)-BS_b(1,2)];

%angle between BS max gain vector and displacement vector in xy plane
phi=180.*acos(dot(D,BS_b(2,1:2))./norm(D))./pi;
%elevation angle
theta=90+tilt-180.*atan(norm(D)./BS_b(1,3))./pi;

if rand<=0.5
           L_entry=-10;
       else
           L_entry=0;
       end

I_dBmUE=SINR_target+noisefloor-F1336(HEXphi(round(rand*1000000)),HEXtheta(round(rand*1000000)),18)-UE_gain-L_body-L_entry+P1546FieldStrMixed(700,50,1.5,40,0,'Suburban',HEXd(round(rand*1000000))./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', 1.5);

if I_dBmUE>max_dBmUE
    I_dBmUE=max_dBmUE;
    n_UEmaxB=n_UEmaxB+1;
end
if I_dBmUE<min_dBmUE
    I_dBmUE=min_dBmUE;
    n_UEminB=n_UEminB+1;
end

I_self(w)=I_self(w)+10.^((I_dBmUE+F1336(phi,theta,18)+UE_gain+L_body+L_entry-P1546FieldStrMixed(700,50,1.5,40,0,'Suburban',norm(D)./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', 1.5))./10);

end
end

I_0=mean(I_self);


%WANTED SIGNAL
S_const=10.^((SINR_target+noisefloor)./10);

%Generate Signal Strength Map (no power control)
n_t=100;
S=NaN(n_t); %(mW)
y_t=linspace(0,-2*r_a,n_t);
x_t=linspace(-sqrt(3)*r_a/2,sqrt(3)*r_a/2,n_t);
for j=1:n_t

    for i=1:n_t
       %x_t=BS_b(1,1)+(i.*4.*r_a./(sqrt(3)*n_t));
       %y_t=BS_b(2,1)-(j.*2.*r_a./n_t);

        %for Gain of BS antenna in direction of UE
       D=[x_t(i),-y_t(j),-h_BSb]; %displacment vector
        phi=acos(dot(D,BS_b(2,:))./norm(D)); %Gmax->BS_b->UE angle


       %d=sqrt((x_t-BS_b(1,1)).^2+(y_t-BS_b(1,2)).^2+(0-BS_b(1,3)).^2); %seperation distance from UE to Base station
       d=sqrt(x_t(i).^2+y_t(j).^2+h_BSb.^2);
       %Ls=20.*log10(lambda./(4.*pi.*d));%P1546FieldStrMixed(700,50,1.5,40,10,'Suburban',d./1000,'Land',0); %p.1546 path loss
       Ls=-P1546FieldStrMixed(700,50,1.5,40,0,'Suburban',d./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', 1.5); %p.1546 path loss

       if rand<=0.5
           L_entry=-10;
       else
           L_entry=0;
       end

       %normalize position to unit hexagon
      % x_tt=x_t-BS_b(1,1);
      % x=(x_tt*(1-sqrt(3)*r_a/2/abs(x_tt)))./r_a;
       %y=(y_t-BS_b(1,2)+r_a)/r_a;
       x=x_t(i)/r_a;
       y=(y_t(j)+r_a)/r_a;

        if x<=0 && y>=0
        if y<=(x/sqrt(3)+1)

            S(i,j)=10.^((S_dBmUE+AntGain(phi,18,"F.1336")+Ls+L_body+L_entry)./10);

        end
    elseif x>0 && y>0
        if y<=(-x/sqrt(3)+1)

            S(i,j)=10.^((S_dBmUE+AntGain(phi,18,"F.1336")+Ls+L_body+L_entry)./10);

        end
    elseif x<0 && y<0
        if y>=-(x/sqrt(3)+1)

           S(i,j)=10.^((S_dBmUE+AntGain(phi,18,"F.1336")+Ls+L_body+L_entry)./10);

        end
    elseif x>0 && y<0
        if y>=-(-x/sqrt(3)+1)

           S(i,j)=10.^((S_dBmUE+AntGain(phi,18,"F.1336")+Ls+L_body+L_entry)./10);

        end
    end



    end
end

L_entry=-10;

UE_a=zeros(3,3,n_a);



%generate n_a random xy drone positions for each cell
UE_a(1:2,1,:)=randHEX(r_a,y1,x1,n_a);
UE_a(1:2,2,:)=randHEX(r_a,y2,x2,n_a);
UE_a(1:2,3,:)=randHEX(r_a,y3,x3,n_a);

%generate n_a random drone altitudes
for i=1:n_a
UE_a(3,1,i)=(rand*270)+30;
UE_a(3,2,i)=(rand*270)+30;
UE_a(3,3,i)=(rand*270)+30;
end



%while (abs(throughput_loss./target)-1)>0.1
 SINR=NaN(n_t,n_t,n_a);
throughput_loss=NaN(n_t,n_t,n_a);

I=zeros(1,n_a); %aggregate Interference power (mW)


%Drone Interference
for i=1:n_a

%calculate LOS displacement vector for each drone in the 3 sectors
%repeat for n_a runs
D1=[UE_a(1,1,i)-BS_b(1,1),UE_a(2,1,i)-BS_b(1,2),UE_a(3,1,i)-BS_b(1,3)];
D2=[UE_a(1,2,i)-BS_b(1,1),UE_a(2,2,i)-BS_b(1,2),UE_a(3,2,i)-BS_b(1,3)];
D3=[UE_a(1,3,i)-BS_b(1,1),UE_a(2,3,i)-BS_b(1,2),UE_a(3,3,i)-BS_b(1,3)];

d1=[UE_a(1,1,i)-BS_a1(1,1),UE_a(2,1,i)-BS_a1(1,2),UE_a(3,1,i)-BS_a1(1,3)];
d2=[UE_a(1,2,i)-BS_a2(1,1),UE_a(2,2,i)-BS_a2(1,2),UE_a(3,2,i)-BS_a2(1,3)];
d3=[UE_a(1,3,i)-BS_a3(1,1),UE_a(2,3,i)-BS_a3(1,2),UE_a(3,3,i)-BS_a3(1,3)];


I_dBmUE1=SINR_target+noisefloor-F1336V(UE_a(:,1,i),BS_a1(1,:),BS_a1(2,:),18)-20.*log10(lambda./(4.*pi.*norm(d1)))-UE_gain;
I_dBmUE2=SINR_target+noisefloor-F1336V(UE_a(:,2,i),BS_a2(1,:),BS_a2(2,:),18)-20.*log10(lambda./(4.*pi.*norm(d2)))-UE_gain;
I_dBmUE3=SINR_target+noisefloor-F1336V(UE_a(:,3,i),BS_a3(1,:),BS_a3(2,:),18)-20.*log10(lambda./(4.*pi.*norm(d3)))-UE_gain;

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

SINR(:,:,i)=S./(I(i)+N+I_0);
throughput_loss(:,:,i)=(log2(1+S./(N+I_0))-log2(1+SINR(:,:,i)))./(log2(1+S./(N+I_0))+log2(1+SINR(:,:,i)));
end
mean_throughput_loss=sum(throughput_loss,3)./n_a;
peak=max(max(mean_throughput_loss));

SINR_0=S_const./(N+I_self);
SINR_const=S_const./(N+I_self+I);
bitloss=100.*(log2(1+SINR_0)-log2(1+SINR_const))./(log2(1+SINR_0)+log2(1+SINR_const));
%end

drones_x=NaN(1,n_a*3);
drones_x(1:n_a)=UE_a(1,1,:);
drones_x(n_a+1:2*n_a)=UE_a(1,2,:);
drones_x(2*n_a+1:3*n_a)=UE_a(1,3,:);

drones_y=NaN(1,n_a*3);
drones_y(1:n_a)=UE_a(2,1,:);
drones_y(n_a+1:2*n_a)=UE_a(2,2,:);
drones_y(2*n_a+1:3*n_a)=UE_a(2,3,:);



figure(1)
subplot(2,2,1)
contourf(mean_throughput_loss)
title('Average Throughput Loss')
xlabel('y step')
ylabel('x step')
colorbar

subplot(2,2,2)
title('Victim UE positions for all time steps')
xlabel('x position (m)')
ylabel('y position (m)')
plot(UE_b(1,:),UE_b(2,:),'o')


subplot(2,2,3)
histogram(bitloss)
title('Distribution of Bitrate loss')
xlabel('Bitrate loss (%)')
ylabel(['instances out of', num2str(n_a)])

subplot(2,2,4)
title('Drone positions for all time steps')
xlabel('x position (m)')
ylabel('y position (m)')
plot(drones_x,drones_y,'.')

toc
