
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
G_0=15; %BS antenna gain (dBi)
d_sep=18000; % (m) seperation distance between both networks
r_a=4000; %cell radius (m)
n_a=100; %number of monte carlo runs

max_dBmUE=36.95 - UE_gain; %max transmiting power of UEs (dBm)
%according to RSS 130

min_dBmUE=-40; %min transmiting power of UEs (dBm)
nf=5; %noise figure (dB)





%Noise
noisefloor=10.*log10(1.380649e-23*290*bw*1000)+5;  %dBm
N_thermal=10.^(noisefloor./10); %mW


%Matrix for storing positions of hexagon corners
%so we can plot the outline later
%7 because there are 7 corner (first + last overlap)
%2 because x&y
%6 because we have two tri-sectors (A&B) that share a number
outline=NaN(7,2,19,6);


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
[UE_a(1:2,1,n,:),outline(:,:,n,1)]=randHEX(r_a,BS_a(2,n)+r_a,BS_a(1,n),n_a);
[UE_a(1:2,2,n,:),outline(:,:,n,2)]=randHEX(r_a,BS_a(2,n)-0.5*r_a,BS_a(1,n)+0.866*r_a,n_a);
[UE_a(1:2,3,n,:),outline(:,:,n,3)]=randHEX(r_a,BS_a(2,n)-0.5*r_a,BS_a(1,n)-0.866*r_a,n_a);
end

%outline(:,:,1)=[BS_a1(1,2)+r_a,BS_a1(1,1)]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
%outline(:,:,1)=[BS_a1(1,2)-0.5*r_a,BS_a1(1,1)-sqrt(3)*r_a/2]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
%outline(:,:,2)=[BS_a1(1,2)-0.5*r_a,BS_a1(1,1)+sqrt(3)*r_a/2]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];
%outline(:,:,3)=[BS_a1(1,2)+r_a,BS_a1(1,1)]+[r_a,0;r_a/2,sqrt(3)*r_a/2;-r_a/2,sqrt(3)*r_a/2;-r_a,0;-r_a/2,-sqrt(3)*r_a/2;r_a/2,-sqrt(3)*r_a/2;r_a,0];



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
[UE_b(1:2,1,n,:),outline(:,:,n,4)]=randHEX(r_a,BS_b(2,n)-r_a,BS_b(1,n),n_a);
[UE_b(1:2,2,n,:),outline(:,:,n,5)]=randHEX(r_a,BS_b(2,n)+0.5*r_a,BS_b(1,n)+0.866*r_a,n_a);
[UE_b(1:2,3,n,:),outline(:,:,n,6)]=randHEX(r_a,BS_b(2,n)+0.5*r_a,BS_b(1,n)-0.866*r_a,n_a);
end




%count number of UEs over or under tx power range
n_UEmaxB=zeros(3,19);
n_UEminB=zeros(3,19);
n_UEmaxA=zeros(3,19);
n_UEminA=zeros(3,19);

%WANTED SIGNAL
S_const=10.^((SINR_target+noisefloor)./10); %(mW)

%Precalculate Loss
for d=1:300000
P1546(d)=-P1546FieldStrMixed(fc./1000000,50,h_BSb,h_UEb,10,'Rural',d./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', h_BSb);
end

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

I_dBmUEa=NaN(3,19,n_a);
I_dBmUEb=NaN(3,19,n_a);
toc

%Loop to calculate tx power for all UEs using power control
for i=1:n_a %iterate over each monte carlo run
for n=1:19  %iterate over each BS
for s=1:3   %iterate over each sector



UE_a(3,s,n,i)=(rand*270)+30; %generate random drone altitude

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
L_b_b(s,n,i)=P1546(round(vecnorm(d_b(s,n,i,1:2))));

I_dBmUEb(s,n,i)=SINR_target+noisefloor-AntGain_b_b(s,n,i)-L_b_b(s,n,i)-UE_gain-L_feeder-L_body;


if I_dBmUEb(s,n,i)>max_dBmUE
    I_dBmUEb(s,n,i)=max_dBmUE;
    n_UEmaxB(s,n)=n_UEmaxB(s,n)+1;
end
if I_dBmUEb(s,n,i)<min_dBmUE
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
L_a_a(s,n,i)=20.*log10(lambda./(4.*pi.*vecnorm(d_a(s,n,i,:))));




I_dBmUEa(s,n,i)=SINR_target+noisefloor-AntGain_a_a(s,n,i)-L_a_a(s,n,i)-UE_gain-L_feeder;


if I_dBmUEa(s,n,i)>max_dBmUE
    I_dBmUEa(s,n,i)=max_dBmUE;
    n_UEmaxA(s,n)=n_UEmaxA(s,n)+1;
end
if I_dBmUEa(s,n,i)<min_dBmUE
    I_dBmUEa(s,n,i)=min_dBmUE;
    n_UEminA(s,n)=n_UEminA(s,n)+1;
end



end
end
end
toc

%Seperation Distance Sweep
%d_s=linspace(1000,200000,200);
%for t=1:200

 % UE_b(2,:,:,:) = UE_b(2,:,:,:)+d_s(t);
  %BS_b(2,:) = BS_b(2,:)+d_s(t);
  
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


for N=1:19   %Change reference BS
for S=1:3    %Change reference sector
    

for n=1:19  %iterate over each BS in both networks
for s=1:3   %iterate over each sector
for i=1:n_a
    
if s~=S || n~=N %don't count self interference from the reference sectors own UE 
    

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
L_b_r(S,N,s,n,i)=P1546(round(vecnorm(D_b(S,N,s,n,i,1:2))));

I_b(S,N,s,n,i)=10.^((I_dBmUEb(s,n,i)+AntGain_b_r(S,N,s,n,i)+L_b_r(S,N,s,n,i)+UE_gain+L_body+L_feeder)./10); %mW





end %end if statment since all sectors in A interfere to each reference sector

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
L_a_r(S,N,s,n,i)=20.*log10(lambda./(4.*pi.*vecnorm(D_a(S,N,s,n,i,:))));


I_a(S,N,s,n,i)=10.^((I_dBmUEa(s,n,i)+AntGain_a_r(S,N,s,n,i)+L_a_r(S,N,s,n,i)+UE_gain+L_feeder)./10); %mW



end
end

end
end

end



%Aggregate all interference (sum along n,s)
I_0=squeeze(sum(sum(I_b,3),4));
I=squeeze(sum(sum(I_a,3),4));

%Baseline SINR
SINR_0=S_const./(N_thermal+I_0);

%SINR with drone interference
SINR=S_const./(N_thermal+I_0+I);

%throughput (b/s)
throughput_0 = 0.4.*bw.*log2(1+SINR_0);
throughput = 0.4.*bw.*log2(1+SINR);

%throughput loss (%)
bitloss=100.*(log2(1+SINR_0)-log2(1+SINR))./(log2(1+SINR_0));

%average over all n_a (gives value for each sector in B network)
throughput_loss_map=squeeze(mean(bitloss,3));

%average over whole network (sectors & BS)
mean_throughput_loss=squeeze(mean(mean(throughput_loss_map,1),2));
%end


network_loss_i=squeeze(mean(mean(bitloss,1),2));
for i=2:n_a
average_loss(i)=mean(network_loss_i(2:i));
n_runs(i)=i;
end


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
%plot(n_runs(3:n_a),average_loss(3:n_a)-average_loss(2:(n_a-1)))
%xlabel("Number of Runs n_a")
%ylabel("Throughput Loss change per run (%)")

%plot(d_s,mean_throughput_loss)
%xlabel("Seperation Distance (m)")
%ylabel("Network Throughput Loss (%)")

toc
