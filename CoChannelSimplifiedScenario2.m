%Main function for simplified scenario 2
%Urban and co-channel
%PhiDistribution.m does not need to be run prior to executing this function

%***********
%to switch between free space and P1546 propagation 
%losses for terestrial UEs, 
%manualy comment/uncomment lines 164-169 & 203-208
%***********



%clear all stored values that may interphere
clear all

%start run timer
tic


fc=708000000; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)
bw=6000000; %bandwidth (Hz)
tilt=-3;  %mechanical downtilt (degrees)
SINR_target=10; %(dB)
L_feeder=-3; %(dB)
L_body=-4; %(dB)
UE_gain=-3; %(dBi)
L_entry=0; %(dB)
h_BSa=30; %BS antenna height (m)
h_BSb=30;  %BS antenna height (m)
h_UEb=1.5; %UE antenna height (m)
G_0=16; %BS antenna gain (dBi)
d_sep=106000; % (m) seperation distance between both networks
r_a=500; %cell radius (m)
n_a=1; %number of monte carlo runs
max_dBmUE=23; %max transmiting power of UEs (dBm)
min_dBmUE=-40; %min transmiting power of UEs (dBm)
nf=5; %noise figure (dB)





%Noise
noisefloor=10.*log10(1.380649e-23*290*bw*1000)+5;  %dBm
N=10.^(noisefloor./10); %mW


%Matrix for storing positions of hexagon corners
%so we can plot the outline later
%7 because there are 7 corner (first + last overlap)
%2 because x&y
%6 because we have two tri-sectors
outline=NaN(7,2,6);

%Interfering network A
%position and antenna orrientation of BS
BS_a1=[0,0,h_BSa;0,1,0];
BS_a2=[0,0,h_BSa;0.866,-0.5,0];
BS_a3=[0,0,h_BSa;-0.866,-0.5,0];


%matrix for storing drone position in all sectors

UE_a=zeros(3,1,n_a);

UE_a(:,1,1)=[0,2*r_a,300];
UE_a(:,2,1)=[sqrt(3)*r_a,-1*r_a,300];
UE_a(:,3,1)=[-sqrt(3)*r_a,-1*r_a,300];

%outline(:,:,1)=[BS_a1(1,2)+r_a,BS_a1(1,1)]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
outline(:,:,1)=[BS_a1(1,2)-0.5*r_a,BS_a1(1,1)-sqrt(3)*r_a/2]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
outline(:,:,2)=[BS_a1(1,2)-0.5*r_a,BS_a1(1,1)+sqrt(3)*r_a/2]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
outline(:,:,3)=[BS_a1(1,2)+r_a,BS_a1(1,1)]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];




%Victim Network B
%three BS antennas at same position and 120 deg spread
BS_b1=[0,d_sep,h_BSb;0,-1,0]; %position and antenna orrientation base station
BS_b2=[0,d_sep,h_BSb;-0.866,0.5,0];
BS_b3=[0,d_sep,h_BSb;0.866,0.5,0];

UE_b=zeros(3,3,n_a);

%2 self intefering UEs (no UE in reference cell)
UE_b(:,2,1)=[-sqrt(3)*r_a+BS_b1(1,1),r_a+BS_b1(1,2),h_UEb];
UE_b(:,3,1)=[sqrt(3)*r_a+BS_b1(1,1),r_a+BS_b1(1,2),h_UEb];

outline(:,:,4)=[BS_b1(1,2)+0.5*r_a,BS_b1(1,1)-sqrt(3)*r_a/2]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
outline(:,:,5)=[BS_b1(1,2)+0.5*r_a,BS_b1(1,1)+sqrt(3)*r_a/2]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
outline(:,:,6)=[BS_b1(1,2)-r_a,BS_b1(1,1)]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];


%Plotting hex outlines and UE positions
hold on
for i=1:6
    plot(outline(:,2,i),outline(:,1,i),'b')
    %plot(UE(1,i),UE(2,i),'k.')
    xlabel('x position (m)')
    ylabel('y position (m)')
end
plot(UE_a(1,:,1),UE_a(2,:,1),'ro','MarkerSize',15)
plot(UE_b(1,2,1),UE_b(2,2,1),'bo','MarkerSize',15)
plot(UE_b(1,3,1),UE_b(2,3,1),'bo','MarkerSize',15)
plot([BS_a1(1,1),BS_b1(1,1)],[BS_a1(1,2),BS_b1(1,2)],'k--')
%legend('victim network','victim network','aerial network','reference sector','Aerial UE','Terestrial UE','','seperation distance')
hold off


%count number of UE over or under tx power range
n_UEmaxA=0;
n_UEmaxB=0;
n_UEminA=0;
n_UEminB=0;

%WANTED SIGNAL
S_const=10.^((SINR_target+noisefloor)./10); %(mW)


I=zeros(1,n_a); %aggregate Interference power (mW)
I_self=zeros(1,n_a); %aggregate self interference power (mW)


for i=1:n_a


%Calculate amount of self interference

%note there are two sets of displacement vectors and angles
%UPPER CASE is wrt to reference cell BS
%lower case is wrt to own BS

%Displacement vector from UE to reference cell

%D1=[UE_b(1,1,i)-BS_b1(1,1),UE_b(2,1,i)-BS_b1(1,2),UE_b(3,1,i)-BS_b1(1,3)];
D2=[UE_b(1,2,i)-BS_b1(1,1),UE_b(2,2,i)-BS_b1(1,2),UE_b(3,2,i)-BS_b1(1,3)];
D3=[UE_b(1,3,i)-BS_b1(1,1),UE_b(2,3,i)-BS_b1(1,2),UE_b(3,3,i)-BS_b1(1,3)];

%angle between BS max gain vector and displacement vector in xy plane
PHI2b=acosd(dot(D2(1:2),BS_b1(2,1:2))./norm(D2(1:2)));
PHI3b=acosd(dot(D3(1:2),BS_b1(2,1:2))./norm(D3(1:2)));

%elevation angle
THETA2b=acosd(dot([norm(D2(1:2)),D2(3)],[norm(BS_b1(2,1:2)),BS_b1(2,3)])./norm(BS_b1(2,:))./norm(D2)).*D2(3)./abs(D2(3));
THETA3b=acosd(dot([norm(D3(1:2)),D3(3)],[norm(BS_b1(2,1:2)),BS_b1(2,3)])./norm(BS_b1(2,:))./norm(D3)).*D3(3)./abs(D3(3));

%displacement vector to own BS
%d1=[UE_b(1,1,i)-BS_b1(1,1),UE_a(2,1,i)-BS_a1(1,2),UE_a(3,1,i)-BS_a1(1,3)];
d2=[UE_b(1,2,i)-BS_b2(1,1),UE_b(2,2,i)-BS_b2(1,2),UE_b(3,2,i)-BS_b2(1,3)];
d3=[UE_b(1,3,i)-BS_b3(1,1),UE_b(2,3,i)-BS_b3(1,2),UE_b(3,3,i)-BS_b3(1,3)];

%angle between own BS max gain vector and displacement vector in xy plane
phi2b=acosd(dot(d2(1:2),BS_b2(2,1:2))./norm(d2(1:2)));
phi3b=acosd(dot(d3(1:2),BS_b3(2,1:2))./norm(d3(1:2)));

%elevation angle
theta2b=acosd(dot([norm(d2(1:2)),d2(3)],[norm(BS_b2(2,1:2)),BS_b2(2,3)])./norm(BS_b2(2,:))./norm(d2)).*d2(3)./abs(d2(3));
theta3b=acosd(dot([norm(d3(1:2)),d3(3)],[norm(BS_b3(2,1:2)),BS_b3(2,3)])./norm(BS_b3(2,:))./norm(d3)).*d3(3)./abs(d3(3));

%Antenna gain to own BS
AntGain_b2_b2=F1336(phi2b,theta2b,G_0,tilt);
AntGain_b3_b3=F1336(phi3b,theta3b,G_0,tilt);

%propagation loss to Own BS
%free space
%L_b2_b2=20.*log10(lambda./(4.*pi.*norm(d2)));
%L_b3_b3=20.*log10(lambda./(4.*pi.*norm(d3)));

%P1546
L_b2_b2=-P1546FieldStrMixed(708,50,h_BSb,h_UEb,0,'Urban',norm(d2(1:2))./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', h_BSb);
L_b3_b3=-P1546FieldStrMixed(708,50,h_BSb,h_UEb,0,'Urban',norm(d2(1:2))./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', h_BSb);


%I_dBmUE1=SINR_target+noisefloor-F1336V(UE_b(:,1,i),BS_b1(1,:),BS_b1(2,:),18)-20.*log10(lambda./(4.*pi.*norm(d1)+1e-3))-UE_gain;
I_dBmUE4=SINR_target+noisefloor-AntGain_b2_b2-L_b2_b2-UE_gain-L_feeder-L_body;
I_dBmUE5=SINR_target+noisefloor-AntGain_b3_b3-L_b3_b3-UE_gain-L_feeder-L_body;



if I_dBmUE4>max_dBmUE
    I_dBmUE4=max_dBmUE;
    n_UEmaxB=n_UEmaxB+1;
end
if I_dBmUE4<min_dBmUE
    I_dBmUE4=min_dBmUE;
    n_UEminB=n_UEminB+1;
end

if I_dBmUE5>max_dBmUE
    I_dBmUE5=max_dBmUE;
    n_UEmaxB=n_UEmaxB+1;
end
if I_dBmUE5<min_dBmUE
    I_dBmUE5=min_dBmUE;
    n_UEminB=n_UEminB+1;
end


%Antenna gain to reference BS
AntGain_b2_b1=F1336(PHI2b,THETA2b,G_0,tilt);
AntGain_b3_b1=F1336(PHI3b,THETA3b,G_0,tilt);

%propagation loss to reference BS
%free space
%L_b2_b1=20.*log10(lambda./(4.*pi.*norm(D2)));
%L_b3_b1=20.*log10(lambda./(4.*pi.*norm(D3)));

%P1546
L_b2_b1=-P1546FieldStrMixed(708,50,h_BSb,h_UEb,0,'Urban',norm(D2(1:2))./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', h_BSb);
L_b3_b1=-P1546FieldStrMixed(708,50,h_BSb,h_UEb,0,'Urban',norm(D3(1:2))./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', h_BSb);


I_self(i)=10.^((I_dBmUE4+AntGain_b2_b1+L_b2_b1+UE_gain+L_body+L_feeder)./10);
I_self(i)=I_self(i)+10.^((I_dBmUE5+AntGain_b3_b1+L_b3_b1+UE_gain+L_body+L_feeder)./10);






%Drone Interference


%calculate LOS displacement vector for each drone in the 3 sectors
%repeat for n_a runs
D1a=[UE_a(1,1,i)-BS_b1(1,1),UE_a(2,1,i)-BS_b1(1,2),UE_a(3,1,i)-BS_b1(1,3)];
D2a=[UE_a(1,2,i)-BS_b1(1,1),UE_a(2,2,i)-BS_b1(1,2),UE_a(3,2,i)-BS_b1(1,3)];
D3a=[UE_a(1,3,i)-BS_b1(1,1),UE_a(2,3,i)-BS_b1(1,2),UE_a(3,3,i)-BS_b1(1,3)];

%angle between BS max gain vector and displacement vector in xy plane
PHI1a=180.*acos(dot(D1a(1:2),BS_b1(2,1:2))./norm(D1a(1:2))./norm(BS_b1(2,1:2)))./pi;
PHI2a=180.*acos(dot(D2a(1:2),BS_b1(2,1:2))./norm(D2a(1:2))./norm(BS_b1(2,1:2)))./pi;
PHI3a=180.*acos(dot(D3a(1:2),BS_b1(2,1:2))./norm(D3a(1:2))./norm(BS_b1(2,1:2)))./pi;

%elevation angle
THETA1a=180.*acos(dot([norm(D1a(1:2)),D1a(3)],[norm(BS_b1(2,1:2)),BS_b1(2,3)])./norm(BS_b1(2,:))./norm(D1a))./pi;
THETA2a=180.*acos(dot([norm(D2a(1:2)),D2a(3)],[norm(BS_b1(2,1:2)),BS_b1(2,3)])./norm(BS_b1(2,:))./norm(D2a))./pi;
THETA3a=180.*acos(dot([norm(D3a(1:2)),D3a(3)],[norm(BS_b1(2,1:2)),BS_b1(2,3)])./norm(BS_b1(2,:))./norm(D3a))./pi;

d1a=[UE_a(1,1,i)-BS_a1(1,1),UE_a(2,1,i)-BS_a1(1,2),UE_a(3,1,i)-BS_a1(1,3)];
d2a=[UE_a(1,2,i)-BS_a1(1,1),UE_a(2,2,i)-BS_a1(1,2),UE_a(3,2,i)-BS_a1(1,3)];
d3a=[UE_a(1,3,i)-BS_a1(1,1),UE_a(2,3,i)-BS_a1(1,2),UE_a(3,3,i)-BS_a1(1,3)];

%angle between BS max gain vector and displacement vector in xy plane
phi1a=180.*acos(dot(d1a(1:2),BS_a1(2,1:2))./norm(d1a(1:2))./norm(BS_a1(2,1:2)))./pi;
phi2a=180.*acos(dot(d2a(1:2),BS_a2(2,1:2))./norm(d2a(1:2))./norm(BS_a2(2,1:2)))./pi;
phi3a=180.*acos(dot(d3a(1:2),BS_a3(2,1:2))./norm(d3a(1:2))./norm(BS_a3(2,1:2)))./pi;

%elevation angle
theta1a=180.*acos(dot([norm(d1a(1:2)),d1a(3)],[norm(BS_a1(2,1:2)),BS_a1(2,3)])./norm(BS_a1(2,:))./norm(d1a))./pi;
theta2a=180.*acos(dot([norm(d2a(1:2)),d2a(3)],[norm(BS_a2(2,1:2)),BS_a2(2,3)])./norm(BS_a2(2,:))./norm(d2a))./pi;
theta3a=180.*acos(dot([norm(d3a(1:2)),d3a(3)],[norm(BS_a3(2,1:2)),BS_a3(2,3)])./norm(BS_a3(2,:))./norm(d3a))./pi;

AntGain_a1_a1=F1336(phi1a,theta1a,G_0,tilt);
AntGain_a2_a2=F1336(phi2a,theta2a,G_0,tilt);
AntGain_a3_a3=F1336(phi3a,theta3a,G_0,tilt);

L_a1_a1=20.*log10(lambda./(4.*pi.*norm(d1a)));
L_a2_a2=20.*log10(lambda./(4.*pi.*norm(d2a)));
L_a3_a3=20.*log10(lambda./(4.*pi.*norm(d3a)));

I_dBmUE1=SINR_target+noisefloor-AntGain_a1_a1-L_a1_a1-UE_gain-L_feeder;
I_dBmUE2=SINR_target+noisefloor-AntGain_a2_a2-L_a2_a2-UE_gain-L_feeder;
I_dBmUE3=SINR_target+noisefloor-AntGain_a3_a3-L_a3_a3-UE_gain-L_feeder;

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

%Antenna gain to reference BS
AntGain_a1_b1=F1336(PHI1a,THETA1a,G_0,tilt);
AntGain_a2_b1=F1336(PHI2a,THETA2a,G_0,tilt);
AntGain_a3_b1=F1336(PHI3a,THETA3a,G_0,tilt);

%propagation loss to reference BS
L_a1_b1=20.*log10(lambda./(4.*pi.*norm(D1a)));
L_a2_b1=20.*log10(lambda./(4.*pi.*norm(D2a)));
L_a3_b1=20.*log10(lambda./(4.*pi.*norm(D3a)));

%Interference power
I(i)=10.^((I_dBmUE1+AntGain_a1_b1+L_a1_b1+UE_gain+L_feeder)./10); %(mW)
I(i)=I(i)+10.^((I_dBmUE2+AntGain_a2_b1+L_a2_b1+UE_gain+L_feeder)./10);
I(i)=I(i)+10.^((I_dBmUE3+AntGain_a3_b1+L_a3_b1+UE_gain+L_feeder)./10);
end



%Averaging over all n_a but in this case n_a=1
I_0=mean(I_self);

%Baseline SINR
SINR_0=S_const./(N+I_self);

%SINR with drone interference
SINR_const=S_const./(N+I_self+I);

%throughput (b/s)
throughput_0 = 0.4.*bw.*log2(1+SINR_0);
throughput = 0.4.*bw.*log2(1+SINR_const);

%throughput loss (%)
bitloss=100.*(log2(1+SINR_0)-log2(1+SINR_const))./(log2(1+SINR_0));

%average over all n_a
average_throughput_loss=mean(bitloss)


%end timer
toc
