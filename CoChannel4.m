
%Main function for version 4 of co channel study
%Urban and co-channel
%PhiDistribution.m does not need to be run prior to executing this function

%***********
%to switch between free space and P1546 propagation 
%losses for terestrial UEs, 
%manualy comment/uncomment lines 
%***********



%clear all stored values that may interphere
clear all

%start run timer
tic

for d=1:80
    d_sep=d*1000
   
%for a=1:60 
 

%Low Band RURAL
fc=708e6; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)
bw=5e6; %bandwidth (Hz)
tilt=-3;  %mechanical downtilt (degrees)
L_feeder=-3; %(dB)
L_body=-4; %(dB)
SINR_target=3.2; %(dB)
UE_gain=-3; %(dBi)
L_entry=-14.1; %(dB)
indoor_ratio=0.5; %ratio of terestrial UEs indoors, between 0 and 1 (0 being all outdoor) 
h_BSa=30; %BS antenna height (m)
h_BSb=30;  %BS antenna height (m)
h_BSc=30;
h_UEb=1.5; %UE antenna height (m)
G_0=15; %BS antenna gain (dBi)
d_sep=11000; % (m) seperation distance between both networks
n_a=2000; %number of monte carlo runs
offset=0; %degrees
r_a=2500; %cell radius (m)


%{
%Low Band URBAN
fc=708e6; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)
bw=5e6; %bandwidth (Hz)
tilt=-3;  %mechanical downtilt (degrees)
L_feeder=-3; %(dB)
L_body=-4; %(dB)
SINR_target=10;%18.5; %(dB)
UE_gain=-3; %(dBi)
L_entry=-14.1; %(dB)
indoor_ratio=0.7; %ratio of terestrial UEs indoors, between 0 and 1 (0 being all outdoor) 
h_BSa=30; %BS antenna height (m)
h_BSb=30;  %BS antenna height (m)
h_BSc=30;
h_UEb=1.5; %UE antenna height (m)
G_0=15; %BS antenna gain (dBi)
d_sep=70000; % (m) seperation distance between both networks
n_a=2000; %number of monte carlo runs
offset=0; %degrees
r_a=1000; %cell radius (m)
%}
%{
%Mid Band RURAL
fc=1900e6; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)
bw=5e6; %bandwidth (Hz)
tilt=-3;  %mechanical downtilt (degrees)
L_feeder=-3; %(dB)
L_body=-4; %(dB)
SINR_target=13.2; %(dB)
UE_gain=0; %(dBi)
L_entry=-14.9; %(dB)
indoor_ratio=0.5; %ratio of terestrial UEs indoors, between 0 and 1 (0 being all outdoor) 
h_BSa=30; %BS antenna height (m)
h_BSb=30;  %BS antenna height (m)
h_BSc=30;
h_UEb=1.5; %UE antenna height (m)
G_0=18; %BS antenna gain (dBi)
d_sep=11000; % (m) seperation distance between both networks
n_a=2000; %number of monte carlo runs
offset=60; %degrees
r_a=2500; %cell radius (m)
%}
%{
%Mid Band URBAN
fc=1900e6; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)
bw=5e6; %bandwidth (Hz)
tilt=-10;  %mechanical downtilt (degrees)
L_feeder=-3; %(dB)
L_body=-4; %(dB)
SINR_target=10;%25.3; %(dB)
UE_gain=0; %(dBi)
L_entry=-14.9; %(dB)
indoor_ratio=0.7; %ratio of terestrial UEs indoors, between 0 and 1 (0 being all outdoor) 
h_BSa=25; %BS antenna height (m)
h_BSb=25;  %BS antenna height (m)
h_BSc=25;
h_UEb=1.5; %UE antenna height (m)
G_0=16; %BS antenna gain (dBi)
d_sep=11000; % (m) seperation distance between both networks
n_a=2000; %number of monte carlo runs
offset=0; %degrees
%r_a=250; %cell radius (m)
r_a=1000;

%}
max_dBmUE=23+10*log10(bw/5e6);
%max_dBmUE=36.95 - UE_gain; %max transmiting power of UEs (dBm)
%according to RSS 130

min_dBmUE=-40+10*log10(bw/5e6); %min transmiting power of UEs (dBm)

nf=5; %noise figure (dB)

%min received power per resource block
P_0=10*log10(bw)-173.8+10*log10(180e3./bw)+nf+SINR_target;

%number of allocated resource blocks
M_pusch=10*log10(round((0.95*bw)./180e3));

%coupling loss precentile
CLp=max_dBmUE-P_0-M_pusch;

%Building entry loss ITU rec P.2109-2
%p=0.5; %probability to exceed 
%F=0; %inverse cumulative normal distribution at p=0.5
%f=fc/1e9; %center freq in GHz
%L_entry=10*log10(10.^(0.1*F*(9.6+2*log10(f))+0.212*abs(tilt)+12.64+3.72*log10(f)+0.96*log10(f).^2)+10.^(0.1*F*(4.5-2*log10(f))+9.1-3*log10(f))+10.^(-0.3));

%Noise
noisefloor=10.*log10(1.380649e-23*290*bw*1000)+5;  %dBm
N_thermal=10.^(noisefloor./10); %mW


%Matrix for storing positions of hexagon corners
%so we can plot the outline later
%7 because there are 7 corner (first + last overlap)
%2 because x&y
%3 because we have tri-sectors
outline_a=NaN(7,2,19,3);
outline_b=NaN(7,2,19,3);


%Interfering network A
%antenna orrientations of BS
BS_aV=[0,1,0;
       0.866,-0.5,0;
       -0.866,-0.5,0];


%Matrix for storing BS positions
BS_a=NaN(3,19);

%set BS Z positions 
BS_a(3,:)=h_BSa;

%set BS X positions
BS_a(1,1:3)=-5.1962*r_a;
BS_a(1,4:7)=-2.5981*r_a;
BS_a(1,8:12)=0*r_a;
BS_a(1,13:16)=2.5981*r_a;
BS_a(1,17:19)=5.1962*r_a;

%set BS Y positions
BS_a(2,8)=0;
BS_a(2,[4,13])=-1.5*r_a;
BS_a(2,[1,9,17])=-3*r_a;
BS_a(2,[5,14])=-4.5*r_a;
BS_a(2,[2,10,18])=-6*r_a;
BS_a(2,[6,15])=-7.5*r_a;
BS_a(2,[3,11,19])=-9*r_a;
BS_a(2,[7,16])=-10.5*r_a;
BS_a(2,12)=-12*r_a;


%matrix for storing drone position in all sectors
%(xyz,sector,BS,run)
UE_a=zeros(3,3,19,n_a);
%note, sector numbering starts from top sector and goes clockwise

%Generate n_a aerial UE positions per sector (also generates outline points
%for sector hexagonal outline
for n=1:19
[UE_a(1:2,1,n,:),outline_a(:,:,n,1)]=randHEX(r_a,BS_a(2,n)+r_a,BS_a(1,n),n_a);
[UE_a(1:2,2,n,:),outline_a(:,:,n,2)]=randHEX(r_a,BS_a(2,n)-0.5*r_a,BS_a(1,n)+0.866*r_a,n_a);
[UE_a(1:2,3,n,:),outline_a(:,:,n,3)]=randHEX(r_a,BS_a(2,n)-0.5*r_a,BS_a(1,n)-0.866*r_a,n_a);
end

%outline(:,:,1)=[BS_a1(1,2)+r_a,BS_a1(1,1)]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
%outline(:,:,1)=[BS_a1(1,2)-0.5*r_a,BS_a1(1,1)-sqrt(3)*r_a/2]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
%outline(:,:,2)=[BS_a1(1,2)-0.5*r_a,BS_a1(1,1)+sqrt(3)*r_a/2]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
%outline(:,:,3)=[BS_a1(1,2)+r_a,BS_a1(1,1)]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];


%Seperation distance sweep


%Victim Network B
%BS have same orientations for whole network 120 deg spread
BS_bV=[0,-1,0;
       0.866,0.5,0;
       -0.866,0.5,0];


%Matrix for storing BS positions
BS_b=NaN(3,19); %(xyz , BS number)

%set BS Z positions
BS_b(3,:)=h_BSb;

%set BS X positions
BS_b(1,1:3)=-5.1962*r_a;
BS_b(1,4:7)=-2.5981*r_a;
BS_b(1,8:12)=0*r_a;
BS_b(1,13:16)=2.5981*r_a;
BS_b(1,17:19)=5.1962*r_a;

%set BS Y positions
BS_b(2,8)=d_sep;
BS_b(2,[4,13])=d_sep+1.5*r_a;
BS_b(2,[1,9,17])=d_sep+3*r_a;
BS_b(2,[5,14])=d_sep+4.5*r_a;
BS_b(2,[2,10,18])=d_sep+6*r_a;
BS_b(2,[6,15])=d_sep+7.5*r_a;
BS_b(2,[3,11,19])=d_sep+9*r_a;
BS_b(2,[7,16])=d_sep+10.5*r_a;
BS_b(2,12)=d_sep+12*r_a;

%Matrix for storing terestrial UE positions
UE_b=zeros(3,3,19,n_a);

%set height of terrestrial UEs
UE_b(3,:,:,:)=h_UEb;

%Generate n_a terrestrial UE positions per sector
for n=1:19
[UE_b(1:2,1,n,:),outline_b(:,:,n,1)]=randHEX(r_a,BS_b(2,n)-r_a,BS_b(1,n),n_a);
[UE_b(1:2,2,n,:),outline_b(:,:,n,2)]=randHEX(r_a,BS_b(2,n)+0.5*r_a,BS_b(1,n)+0.866*r_a,n_a);
[UE_b(1:2,3,n,:),outline_b(:,:,n,3)]=randHEX(r_a,BS_b(2,n)+0.5*r_a,BS_b(1,n)-0.866*r_a,n_a);
end


%Rotation of network A in XY with BS_b(:,10) as axis
%offset=0; %degrees
R=[cosd(offset),-sind(offset);
    sind(offset),cosd(offset)];

%move center of B network to origin
BS_b(2,:)=BS_b(2,:)-(d_sep);
UE_b(2,:,:,:)=UE_b(2,:,:,:)-(d_sep);
outline_b(:,2,:,:)=outline_b(:,2,:,:)-(d_sep);

%Apply rotation opperator
for s=1:3
    temp=R*squeeze(BS_bV(s,1:2))';
    BS_bV(s,1:2)=temp';
    for n=1:19
        toomp=R*squeeze(BS_b(1:2,n));
        BS_b(1:2,n)=toomp;
        for p=1:7
            broomp=R*squeeze(outline_b(p,:,n,s)');
            outline_b(p,:,n,s)=broomp';
        end
        for i=1:n_a
          keith=R*squeeze(UE_b(1:2,s,n,i));
          UE_b(1:2,s,n,i)=keith;
        end
    end
end

%move network to original position
BS_b(2,:)=BS_b(2,:)+(d_sep);
UE_b(2,:,:,:)=UE_b(2,:,:,:)+(d_sep);
outline_b(:,2,:,:)=outline_b(:,2,:,:)+(d_sep);

%{
%Axis of rotation
offset=90; %degrees
x0=0;
y0=6*r_a+d_sep;

UE_b(1,:,:,:) = (cosd(offset).*(UE_b(1,:,:,:)-x0)) - (sind(offset).*(UE_b(2,:,:,:)-y0)) + x0;
UE_b(2,:,:,:) = (sind(offset).*(UE_b(1,:,:,:)-x0)) + (cosd(offset).*(UE_b(2,:,:,:)-y0)) + y0;

BS_b(1,:) = (cosd(offset).*(BS_b(1,:)-x0)) - (sind(offset).*(BS_b(2,:)-y0)) + x0;
BS_b(2,:) = (sind(offset).*(BS_b(1,:)-x0)) + (cosd(offset).*(BS_b(2,:)-y0)) + y0;

outline_b(:,1,:,:) = (cosd(offset).*(outline_b(:,1,:,:)-x0)) - (sind(offset).*(outline_b(:,2,:,:)-y0)) + x0;
outline_b(:,2,:,:) = (sind(offset).*(outline_b(:,1,:,:)-x0)) + (cosd(offset).*(outline_b(:,2,:,:)-y0)) + y0;
%}
%count number of UEs over or under tx power range
n_UEmaxB=zeros(3,19);
n_UEminB=zeros(3,19);
n_UEmaxA=zeros(3,19);
n_UEminA=zeros(3,19);

%WANTED SIGNAL
S_const=10.^((SINR_target+noisefloor)./10); %(mW)

Signal=zeros(3,19,n_a); %incase power control cannot meet target SINR

%Precalculate Loss
%for d=1:300000
%P1546(d)=-P1546FieldStrMixed(fc./1000000,50,h_UEb,h_BSb,10,'Rural',d./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', h_UEb);
%end
%Dat=load("P1546_708MHz_50%_1.5m_30m_10m_rural_1to3e5m.mat");
%P1546=Dat.P1546;

Dat=load("P452_10percent_Below22m.npy.mat");
P452=Dat.P452_10percent_Below22m;

temp=load("P528_0to300km_30m_to_300m_708MHz_0_10%.mat");
L528=temp.L528;
d528=temp.d;
h528=temp.ht;

%Initialize variables
d_b=NaN(3,19,n_a,3);
phi_b=NaN(3,19,n_a);
theta_b=NaN(3,19,n_a);
AntGain_b_b=NaN(3,19,n_a);
L_b_b=NaN(3,19,n_a);

d_a=NaN(3,19,n_a,3);
phi_a=NaN(3,19,n_a);
theta_a=NaN(3,19,n_a);
AntGain_a_a=NaN(3,19,n_a);
L_a_a=NaN(3,19,n_a);
indoor=NaN(3,19,n_a);

I_dBmUEa=NaN(3,19,n_a);
I_dBmUEb=NaN(3,19,n_a);

%UE_a(3,:,:,:)=h_UEb;
%Loop to calculate tx power for all UEs using power control
for i=1:n_a %iterate over each monte carlo run
for n=1:19  %iterate over each BS
for s=1:3   %iterate over each sector

   
if rand>0.33
   UE_a(3,s,n,i)=1.5;
else
UE_a(3,s,n,i)=rand*300; %generate random drone altitude
end
    

%displacement vector to own BS
d_b(s,n,i,:)=[UE_b(1,s,n,i)-BS_b(1,n),UE_b(2,s,n,i)-BS_b(2,n),UE_b(3,s,n,i)-BS_b(3,n)];

%angle between own BS max gain vector and displacement vector in xy plane
phi_b(s,n,i)=acosd(dot(squeeze(d_b(s,n,i,1:2)),BS_bV(s,1:2))./norm(BS_bV(s,1:2))./vecnorm(d_b(s,n,i,1:2)));

%elevation angle
theta_b(s,n,i)=acosd(dot([vecnorm(d_b(s,n,i,1:2)),d_b(s,n,i,3)],[norm(BS_bV(s,1:2)),BS_bV(s,3)])./norm(BS_bV(s,:))./vecnorm(d_b(s,n,i,:))).*d_b(s,n,i,3)./abs(d_b(s,n,i,3));

%Antenna gain to own BS
AntGain_b_b(s,n,i)=F1336(phi_b(s,n,i),theta_b(s,n,i),G_0,tilt);

%propagation loss to Own BS
%free space
%L_b_b(s,n,i)=20.*log10(lambda./(4.*pi.*vecnorm(d_b(s,n,i,:))));

%P1546
L_b_b(s,n,i)=-P452(round(vecnorm(d_b(s,n,i,1:2))),1);


if rand<indoor_ratio
indoor(s,n,i)=0;
I_dBmUEb(s,n,i)=min([max_dBmUE,M_pusch+P_0-AntGain_b_b(s,n,i)-L_b_b(s,n,i)-UE_gain-L_feeder-L_body]);

else
indoor(s,n,i)=1;
I_dBmUEb(s,n,i)=min([max_dBmUE,M_pusch+P_0-AntGain_b_b(s,n,i)-L_b_b(s,n,i)-UE_gain-L_feeder-L_body-L_entry]);

end

if I_dBmUEb(s,n,i)>=max_dBmUE
    I_dBmUEb(s,n,i)=max_dBmUE;
    n_UEmaxB(s,n)=n_UEmaxB(s,n)+1;
end
if I_dBmUEb(s,n,i)<=min_dBmUE
    I_dBmUEb(s,n,i)=min_dBmUE;
    n_UEminB(s,n)=n_UEminB(s,n)+1;
end



%displacement vector to own BS
d_a(s,n,i,:)=[UE_a(1,s,n,i)-BS_a(1,n),UE_a(2,s,n,i)-BS_a(2,n),UE_a(3,s,n,i)-BS_a(3,n)];

%angle between own BS max gain vector and displacement vector in xy plane
phi_a(s,n,i)=acosd(dot(squeeze(d_a(s,n,i,1:2)),BS_aV(s,1:2))./norm(BS_aV(s,1:2))./vecnorm(d_a(s,n,i,1:2)));

%elevation angle
theta_a(s,n,i)=acosd(dot([vecnorm(d_a(s,n,i,1:2)),d_a(s,n,i,3)],[norm(BS_aV(s,1:2)),BS_aV(s,3)])./norm(BS_aV(s,:))./vecnorm(d_a(s,n,i,:))).*d_a(s,n,i,3)./abs(d_a(s,n,i,3));

%Antenna gain to own BS
AntGain_a_a(s,n,i)=F1336(phi_a(s,n,i),theta_a(s,n,i),G_0,tilt);

%propagation loss to Own BS
%free space
%L_a_a(s,n,i)=20.*log10(lambda./(4.*pi.*vecnorm(d_a(s,n,i,:))));

if UE_a(3,s,n,i)>22
%ITU P528
%[~, Hindex] = min(abs(h528-UE_a(3,s,n,i)));
%[~, Dindex] = min(abs(d528-vecnorm(d_a(s,n,i,1:2))));
%temp=tl_p528(vecnorm(d_a(s,n,i,1:2))./1000,h_BSa,UE_a(3,s,n,i),fc./1e6,1,10);
%L_a_a(s,n,i)=temp.A__db;
    L_a_a(s,n,i)=-L528(Hindex,Dindex);
elseif UE_a(3,s,n,i)<22 && UE_a(3,s,n,i)>=19
    L_a_a(s,n,i)=-P452(round(vecnorm(d_a(s,n,i,1:2))+1),5);
elseif UE_a(3,s,n,i)<19 && UE_a(3,s,n,i)>=14
    L_a_a(s,n,i)=-P452(round(vecnorm(d_a(s,n,i,1:2))+1),4);
elseif UE_a(3,s,n,i)<14 && UE_a(3,s,n,i)>=9
    L_a_a(s,n,i)=-P452(round(vecnorm(d_a(s,n,i,1:2))+1),3);
elseif UE_a(3,s,n,i)<9 && UE_a(3,s,n,i)>=4
    L_a_a(s,n,i)=-P452(round(vecnorm(d_a(s,n,i,1:2))+1),2);
elseif UE_a(3,s,n,i)<4
    L_a_a(s,n,i)=-P452(round(vecnorm(d_a(s,n,i,1:2))+1),1);
end

I_dBmUEa(s,n,i)=min([max_dBmUE,M_pusch+P_0-AntGain_a_a(s,n,i)-L_a_a(s,n,i)-UE_gain-L_feeder]);

if I_dBmUEa(s,n,i)>=max_dBmUE
    I_dBmUEa(s,n,i)=max_dBmUE;
    n_UEmaxA(s,n)=n_UEmaxA(s,n)+1;
end
if I_dBmUEa(s,n,i)<=min_dBmUE
    I_dBmUEa(s,n,i)=min_dBmUE;
    n_UEminA(s,n)=n_UEminA(s,n)+1;
end


end

end


end



  
%note there are two sets of displacement vectors and angles
%UPPER CASE is wrt to reference cell BS 
%lower case is wrt to own BS (calculated above already)
D_b=NaN(3,19,3,19,n_a,3);
PHI_b=NaN(3,19,3,19,n_a);
THETA_b=NaN(3,19,3,19,n_a);
AntGain_b_r=NaN(3,19,3,19,n_a);
L_b_r=NaN(3,19,3,19,n_a);

D_a=NaN(3,19,3,19,n_a,3);
PHI_a=NaN(3,19,3,19,n_a);
THETA_a=NaN(3,19,3,19,n_a);
AntGain_a_r=NaN(3,19,3,19,n_a);
L_a_r=NaN(3,19,3,19,n_a);

I_b=zeros(3,19,3,19,n_a); %self interference from tertrial network
I_a=zeros(3,19,3,19,n_a); %Interference from aerial network


%for N=1:19   %Change reference BS
%for S=1:3    %Change reference sector
 N=8;
 S=1;

for n=1:19  %iterate over each BS in both networks
for s=1:3   %iterate over each sector
for i=1:n_a
    

    

%Calculate self interference

%displacement vector to reference BS
D_b(S,N,s,n,i,:)=[UE_b(1,s,n,i)-BS_b(1,N),UE_b(2,s,n,i)-BS_b(2,N),UE_b(3,s,n,i)-BS_b(3,N)];

%angle between ref BS max gain vector and displacement vector in xy plane
PHI_b(S,N,s,n,i)=acosd(dot(squeeze(D_b(S,N,s,n,i,1:2)),BS_bV(S,1:2))./norm(BS_bV(S,1:2))./vecnorm(D_b(S,N,s,n,i,1:2)));

%elevation angle
THETA_b(S,N,s,n,i)=acosd(dot([vecnorm(D_b(S,N,s,n,i,1:2)),D_b(S,N,s,n,i,3)],[norm(BS_bV(S,1:2)),BS_bV(S,3)])./norm(BS_bV(S,:))./vecnorm(D_b(S,N,s,n,i,:))).*D_b(S,N,s,n,i,3)./abs(D_b(S,N,s,n,i,3));

%Antenna gain to ref BS
AntGain_b_r(S,N,s,n,i)=F1336(PHI_b(S,N,s,n,i),THETA_b(S,N,s,n,i),G_0,tilt);

%propagation loss to ref BS
%free space
%L_b_r(S,N,s,n,i)=20.*log10(lambda./(4.*pi.*vecnorm(D_b(S,N,s,n,i,:))));

%P1546
L_b_r(S,N,s,n,i)=-P452(round(vecnorm(D_b(S,N,s,n,i,1:2))+1),1);



if indoor(s,n,i)==0
I_b(S,N,s,n,i)=10.^((I_dBmUEb(s,n,i)+AntGain_b_r(S,N,s,n,i)+L_b_r(S,N,s,n,i)+UE_gain+L_body+L_feeder)./10); %mW
elseif indoor(s,n,i)==1
I_b(S,N,s,n,i)=10.^((I_dBmUEb(s,n,i)+AntGain_b_r(S,N,s,n,i)+L_b_r(S,N,s,n,i)+UE_gain+L_body+L_feeder+L_entry)./10); %mW
end

if s==S && n==N %don't count self interference from the reference sectors own UE
Signal(S,N,i)=I_b(S,N,s,n,i);
I_b(S,N,s,n,i)=0;

end %end if statment since all sectors in A interfere to each reference sector
  
% 50% LOADING factor
if rand<0.5 
   I_b(S,N,s,n,i)=0;
end

%Calculate Interference

    
%displacement vector to ref BS
D_a(S,N,s,n,i,:)=[UE_a(1,s,n,i)-BS_b(1,N),UE_a(2,s,n,i)-BS_b(2,N),UE_a(3,s,n,i)-BS_b(3,N)];

%angle between ref BS max gain vector and displacement vector in xy plane
PHI_a(S,N,s,n,i)=acosd(dot(squeeze(D_a(S,N,s,n,i,1:2)),BS_bV(S,1:2))./norm(BS_bV(S,1:2))./vecnorm(D_a(S,N,s,n,i,1:2)));

%elevation angle
THETA_a(S,N,s,n,i)=acosd(dot([vecnorm(D_a(S,N,s,n,i,1:2)),D_a(S,N,s,n,i,3)],[norm(BS_bV(S,1:2)),BS_bV(S,3)])./norm(BS_bV(S,:))./vecnorm(D_a(S,N,s,n,i,:))).*D_a(S,N,s,n,i,3)./abs(D_a(S,N,s,n,i,3));

%Antenna gain to ref BS
AntGain_a_r(S,N,s,n,i)=F1336(PHI_a(S,N,s,n,i),THETA_a(S,N,s,n,i),G_0,tilt);

%propagation loss to ref BS
%free space
%L_a_r(S,N,s,n,i)=20.*log10(lambda./(4.*pi.*vecnorm(D_a(S,N,s,n,i,:))));

if UE_a(3,s,n,i)>22
%ITU P528
[~, Hindex] = min(abs(h528-UE_a(3,s,n,i)));
[~, Dindex] = min(abs(d528-vecnorm(D_a(S,N,s,n,i,1:2))));

%temp=tl_p528(vecnorm(D_a(s,n,i,1:2))./1000,h_BSa,UE_a(3,s,n,i),fc./1e6,1,10);
%L_a_r(S,N,s,n,i)=temp.A__db;
    L_a_r(S,N,s,n,i)=-L528(Hindex,Dindex);
elseif UE_a(3,s,n,i)<22 && UE_a(3,s,n,i)>=19
    L_a_r(S,N,s,n,i)=-P452(round(vecnorm(D_a(S,N,s,n,i,1:2))+1),5);
elseif UE_a(3,s,n,i)<19 && UE_a(3,s,n,i)>=14
    L_a_r(S,N,s,n,i)=-P452(round(vecnorm(D_a(S,N,s,n,i,1:2))+1),4);
elseif UE_a(3,s,n,i)<14 && UE_a(3,s,n,i)>=9
    L_a_r(S,N,s,n,i)=-P452(round(vecnorm(D_a(S,N,s,n,i,1:2))+1),3);
elseif UE_a(3,s,n,i)<9 && UE_a(3,s,n,i)>=4
    L_a_r(S,N,s,n,i)=-P452(round(vecnorm(D_a(S,N,s,n,i,1:2))+1),2);
elseif UE_a(3,s,n,i)<4
    L_a_r(S,N,s,n,i)=-P452(round(vecnorm(d_a(s,n,i,1:2))+1),1);
end

I_a(S,N,s,n,i)=10.^((I_dBmUEa(s,n,i)+AntGain_a_r(S,N,s,n,i)+L_a_r(S,N,s,n,i)+UE_gain+L_feeder)./10); %mW

% 50% LOADING factor
if rand<0.5 
   I_a(S,N,s,n,i)=0;
end

end

end
end

%end
%end







%Aggregate all interference (sum along n,s)
I_0=squeeze(sum(sum(I_b,3),4));
I=squeeze(sum(sum(I_a,3),4));

%I over Thermal
IoT_map=squeeze(mean(10*log10((I+N_thermal)./N_thermal),3));

%I over N
IoN=10*log10(I./N_thermal);

IoN_map=squeeze(mean(IoN,3));

%baseline I over N 
IoN_0=10*log10(I_0./N_thermal);

IoN_0_map=squeeze(mean(IoN_0,3));

%Baseline SINR
SINR_0=Signal./(N_thermal+I_0);

%SINR with drone interference
SINR=Signal./(N_thermal+I_0+I);


%throughput (b/s)
throughput_0 = 0.4.*bw.*log2(1+SINR_0);
throughput = 0.4.*bw.*log2(1+SINR);

%throghput threshold
throughput_max=0.4.*bw.*log2(1+10^2.2);
throughput_min=0.4.*bw.*log2(1+10^-1);

idx=throughput>throughput_max;
throughput(idx)=throughput_max;

idx=throughput<throughput_min;
throughput(idx)=0;

idx=throughput_0>throughput_max;
throughput_0(idx)=throughput_max;

idx=throughput_0<throughput_min;
throughput_0(idx)=0;

%throughput loss (%)
bitloss=100.*(log2(1+SINR_0)-log2(1+SINR))./(log2(1+SINR_0));

%loss_scatter(d,:)=bitloss(1,8,:);
%average over all n_a (gives value for each sector in B network)
throughput_loss_map=squeeze(mean(bitloss,3));

%average over whole network (sectors & BS)
mean_throughput_loss(d)=squeeze(mean(mean(throughput_loss_map,1),2));
%end


%[CDF_IoN_P,CDF_IoN_V]=homemade_ecdf(reshape(IoN,[1,57*n_a]));
%[CDF_IoN_0_P,CDF_IoN_0_V]=homemade_ecdf(reshape(IoN_0,[1,57*n_a]));
%[CDF_SINR_P,CDF_SINR_V]=homemade_ecdf(reshape(SINR,[1,57*n_a]));
%[CDF_SINR_0_P,CDF_SINR_0_V]=homemade_ecdf(reshape(SINR_0,[1,57*n_a]));
%[CDF_throughput_P,CDF_throughput_V]=homemade_ecdf(reshape(throughput,[1,57*n_a]));
%[CDF_throughput_0_P,CDF_throughput_0_V]=homemade_ecdf(reshape(throughput_0,[1,57*n_a]));
%[CDF_loss_sys_P,CDF_loss_sys_V]=homemade_ecdf(reshape(bitloss,[1,57*n_a]));

[CDF_IoN_ref_P, CDF_IoN_ref_V]=homemade_ecdf(reshape(IoN(1,8,:),[1,n_a]));
[CDF_loss_ref_P, CDF_loss_ref_V]=homemade_ecdf(reshape(bitloss(1,8,:),[1,n_a]));

loss_ref_50th=throughput_loss_map(1,8);
IoN_ref_50th=IoN_map(1,8);

loss_ref_95th(d)=CDF_loss_ref_V(round(0.95*(n_a+1)));
IoN_ref_95th=CDF_IoN_ref_V(round(0.95*(n_a+1)));

%loss_sys_95th=CDF_loss_sys_V(round(0.95*(n_a+1)));
end



%Plotting hex outlines and UE positions
cmap=jet;
hexValues=throughput_loss_map;
hexColor=size(cmap,1)*(hexValues-min(hexValues)+1)./(max(hexValues)-min(hexValues));

hold on
for n=1:19
for s=1:3
    
    
    
 plot(outline_b(:,1,n,s),outline_b(:,2,n,s),'b')
 plot(outline_a(:,1,n,s),outline_a(:,2,n,s),'r')
 
 if s==1 && n==8
     fill(outline_b(:,1,n,s),outline_b(:,2,n,s),'r')
 end
 %{
    for i=1:n_a
      if s~=1 || n~=8
     scatter(UE_b(1,s,n,i),UE_b(2,s,n,i),[],10.*log10(I_b(1,8,s,n,i)),'filled')
      end
     scatter(UE_a(1,s,n,i),UE_a(2,s,n,i),[],10.*log10(I_a(1,8,s,n,i)),'filled')
    %plot(UE(1,i),UE(2,i),'k.')
    end
    %}
end
end
%plot(squeeze(UE_b(1,2,8,:)),squeeze(UE_b(2,2,8,:)))
xlabel('x position (m)')
ylabel('y position (m)')
hold off

%{
%Plotting hex outlines and UE positions
hold on
for n=1:19
for s=1:3
 plot(outline(:,2,n,s),outline(:,1,n,s),'b')
 plot(outline(:,2,n,s+3),outline(:,1,n,s+3),'r')
 
    for i=1:n_a
      if s~=1 || n~=8
     scatter(UE_b(1,s,n,i),UE_b(2,s,n,i),[],10.*log10(I_b(1,8,s,n,i)),'filled')
      end
     scatter(UE_a(1,s,n,i),UE_a(2,s,n,i),[],10.*log10(I_a(1,8,s,n,i)),'filled')
    %plot(UE(1,i),UE(2,i),'k.')
    end
    
end
end
xlabel('x position (m)')
ylabel('y position (m)')

hold off
%}
%plot(n_runs(3:n_a),average_loss(3:n_a)-average_loss(2:(n_a-1)))
%xlabel("Number of Runs n_a")
%ylabel("Throughput Loss change per run (%)")

%plot(d_s,mean_throughput_loss)
%xlabel("Seperation Distance (m)")
%ylabel("Network Throughput Loss (%)")

toc
