
tic

fc=700000000; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)
bw=6000000; %bandwidth (Hz)
tilt=-6;  %mechanical downtilt (degrees)
SINR_target=10; %(dB)
nf=5; %noise factor
L_feeder=-3; %(dB)
L_body=-1; %(dB)
UE_gain=-3; %(dBi)
L_entry=-20; %(dB)
h_BSa=30; %antenna height (m)
h_BSb=30;  %antenna height (m)
G_0=15; %BS antenna gain (dBi)
d_sep=3900; % (m) seperation distance between both networks
% Note d_sep=0 coresponds to cells touching, not overlaping
r_a=300; %cell radius (m)
n_a=1; %number of monte carlo runs
max_dBmUE=23; %max transmiting power of UEs (dBm)
min_dBmUE=-40; %min transmiting power of UEs (dBm)
nf=5; %noise figure (dB)





%Noise
noisefloor=10.*log10(1.380649e-23*290*bw*1000)+5;  %dBm
N=10.^(noisefloor./10); %mW

outline=NaN(7,2,4);

%Interfering network A

%position and antenna orrientation

BS_a1=[0,0,h_BSa;0,1,sin(pi*tilt/180)];


%Center positions of each hex cell
x1=0+BS_a1(1,1);
y1=r_a+BS_a1(1,2);


%matrix for storing drone position in all three sectors
UE_a=zeros(3,1,n_a);

%generate n_a random xy drone positions for each cell


UE_a(:,1,1)=[0,2*r_a,300];
outline(:,:,1)=[BS_a1(1,2)+r_a,BS_a1(1,1)]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];




%Victim Network B
BS_b1=[0,d_sep,h_BSb;0,-1,sin(pi*tilt/180)]; %position and antenna orrientation base station
BS_b2=[0,d_sep,h_BSb;-0.866,0.5,sin(pi*tilt/180)];
BS_b3=[0,d_sep,h_BSb;0.866,0.5,sin(pi*tilt/180)];

UE_b=zeros(3,2,n_a);

UE_b(:,1,:)=[sqrt(3)*r_a+BS_b1(1,1),r_a+BS_b1(1,2),1.5];
UE_b(:,2,:)=[-sqrt(3)*r_a+BS_b1(1,1),r_a+BS_b1(1,2),1.5];

outline(:,:,2)=[BS_b1(1,2)+0.5*r_a,BS_b1(1,1)-sqrt(3)*r_a/2]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
outline(:,:,3)=[BS_b1(1,2)+0.5*r_a,BS_b1(1,1)+sqrt(3)*r_a/2]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
outline(:,:,4)=[BS_b1(1,2)-r_a,BS_b1(1,1)]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];



hold on
for i=2:3
    plot(outline(:,2,i),outline(:,1,i),'b')
    %plot(UE(1,i),UE(2,i),'k.')
    xlabel('x position (m)')
    ylabel('y position (m)')
end
plot(outline(:,2,1),outline(:,1,1),'r')
plot(outline(:,2,4),outline(:,1,4),'g')
plot(UE_a(1,1,1),UE_a(2,1,1),'ro','MarkerSize',15)
plot(UE_b(1,1,1),UE_b(2,1,1),'bo','MarkerSize',15)
plot(UE_b(1,2,1),UE_b(2,2,1),'bo','MarkerSize',15)
plot([BS_a1(1,1),BS_b1(1,1)],[BS_a1(1,2),BS_b1(1,2)],'k--')
legend('victim network','victim network','aerial network','reference sector','Aerial UE','Terestrial UE','','seperation distance')
hold off

%count number of UE over or under tx power range
n_UEmaxA=0;
n_UEmaxB=0;
n_UEminA=0;
n_UEminB=0;

%WANTED SIGNAL
S_const=10.^((SINR_target+noisefloor)./10);


I=zeros(1,n_a); %aggregate Interference power (mW)
I_self=zeros(1,n_a);

%Calculate amount of self interference
for i=1:n_a


%Calculate amount of self interference

for w=2:length(UE_b(1,:))
    %calculate LOS displacement vector in xy plane
D=[UE_b(1,w)-BS_b1(1,1),UE_b(2,w)-BS_b1(1,2)];

%angle between BS max gain vector and displacement vector in xy plane
phi=180.*acos(dot(D,BS_b1(2,1:2))./norm(D))./pi;
%elevation angle
theta=90+tilt-180.*atan(norm(D)./BS_b1(1,3))./pi;

if rand<=0.5
           L_entry=-10;
       else
           L_entry=0;
       end

I_dBmUE=SINR_target+noisefloor-F1336(HEXphi(round(rand*999999)+1),HEXtheta(round(rand*999999)+1),G_0)-UE_gain-L_body-L_entry+P1546FieldStrMixed(700,50,1.5,30,15,'Urban',HEXd(round(rand*999999)+1)./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', 1.5);

if I_dBmUE>max_dBmUE
    I_dBmUE=max_dBmUE;
    n_UEmaxB=n_UEmaxB+1;
end
if I_dBmUE<min_dBmUE
    I_dBmUE=min_dBmUE;
    n_UEminB=n_UEminB+1;
end

I_self(i)=I_self(i)+10.^((I_dBmUE+F1336(phi,theta,G_0)+UE_gain+L_body+L_entry-P1546FieldStrMixed(700,50,1.5,30,15,'Urban',norm(D)./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', 1.5))./10);

end





%Drone Interference


%calculate LOS displacement vector for each drone in the 3 sectors
%repeat for n_a runs
D1=[UE_a(1,1,i)-BS_b1(1,1),UE_a(2,1,i)-BS_b1(1,2),UE_a(3,1,i)-BS_b1(1,3)];


d1=[UE_a(1,1,i)-BS_a1(1,1),UE_a(2,1,i)-BS_a1(1,2),UE_a(3,1,i)-BS_a1(1,3)];



I_dBmUE1=SINR_target+noisefloor-F1336V(UE_a(:,1,i),BS_a1(1,:),BS_a1(2,:),G_0)-20.*log10(lambda./(4.*pi.*norm(d1)))-UE_gain;

if I_dBmUE1>max_dBmUE
    I_dBmUE1=max_dBmUE;
    n_UEmaxA=n_UEmaxA+1;
end
if I_dBmUE1<min_dBmUE
    I_dBmUE1=min_dBmUE;
    n_UEminA=n_UEminA+1;
end


I(i)=10.^((I_dBmUE1+F1336V(UE_a(:,1,i),BS_b1(1,:),BS_b1(2,:),G_0)+20.*log10(lambda./(4.*pi.*norm(D1)))+UE_gain)./10);

end

I_0=mean(I_self);
SINR_0=S_const./(N+I_self);
SINR_const=S_const./(N+I_self+I);
bitloss=100.*(log2(1+SINR_0)-log2(1+SINR_const))./(log2(1+SINR_0)+log2(1+SINR_const));

average_throughput_loss=mean(bitloss)




toc
