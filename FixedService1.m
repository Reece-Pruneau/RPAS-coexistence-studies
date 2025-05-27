

%clear all stored values that may interphere
clear all

%start run timer
tic

%for d=1:16
 %  tilt=-1*(d-1)
   
%for a=1:60 




%{
%SPECTRA Low Band URBAN
fc=708e6; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)
bw=5e6; %bandwidth (Hz)
tilt=0;  %mechanical downtilt (degrees)
L_feeder=0; %(dB)
L_body=-4; %(dB)
SINR_target=18.5; %(dB)
UE_gain=-3; %(dBi)
L_entry=-14.1; %(dB)
indoor_ratio=0.7; %ratio of terestrial UEs indoors, between 0 and 1 (0 being all outdoor) 
h_BSa=44.4; %BS antenna height (m)
h_BSb=44.4;  %BS antenna height (m)
h_BSc=44.4;
h_UEb=1.5; %UE antenna height (m)
G_0=15.9; %BS antenna gain (dBi)
%d_sep=50000; % (m) seperation distance between both networks
n_a=1000; %number of monte carlo runs
offset=0; %degrees
r_a=1000; %cell radius (m)
alpha=0.8;
loading_factor=0.5;
aerial_factor=1;
%}


%{
%SPECTRA Low Band RURAL
fc=708e6; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)
bw=5e6; %bandwidth (Hz)
tilt=0;  %mechanical downtilt (degrees)
L_feeder=0; %(dB)
L_body=-4; %(dB)
SINR_target=10; %(dB)
UE_gain=-3; %(dBi)
L_entry=-14.1; %(dB)
indoor_ratio=0.5; %ratio of terestrial UEs indoors, between 0 and 1 (0 being all outdoor) 
h_BSa=76; %BS antenna height (m)
h_BSb=76;  %BS antenna height (m)
h_BSc=76;
h_UEb=1.5; %UE antenna height (m)
G_0=16.1; %BS antenna gain (dBi)
%d_sep=10000; % (m) seperation distance between both networks
n_a=1000; %number of monte carlo runs
offset=0; %degrees
r_a=2500; %cell radius (m)
alpha=1;
loading_factor=0.5;
aerial_factor=0.33;
%}

%{
%Low Band RURAL
fc=708e6; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)
bw=5e6; %bandwidth (Hz)
tilt=-3;  %mechanical downtilt (degrees)
L_feeder=-3; %(dB)
L_body=-4; %(dB)
SINR_target=10; %(dB)
UE_gain=-3; %(dBi)
L_entry=-14.1; %(dB)
indoor_ratio=0.5; %ratio of terestrial UEs indoors, between 0 and 1 (0 being all outdoor) 
h_BSa=30; %BS antenna height (m)
h_BSb=30;  %BS antenna height (m)
h_BSc=30;
h_UEb=1.5; %UE antenna height (m)
G_0=15; %BS antenna gain (dBi)
d_sep=64000; % (m) seperation distance between both networks
n_a=200; %number of monte carlo runs
offset=30; %degrees
r_a=2500; %cell radius (m)
alpha=1;
loading_factor=0.5;
aerial_factor=0.33;
%}

%{
%Low Band URBAN
fc=708e6; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)
bw=5e6; %bandwidth (Hz)
tilt=-3;  %mechanical downtilt (degrees)
L_feeder=-3; %(dB)
L_body=-4; %(dB)
SINR_target=18.5; %(dB)
UE_gain=-3; %(dBi)
L_entry=-14.1; %(dB)
indoor_ratio=0.7; %ratio of terestrial UEs indoors, between 0 and 1 (0 being all outdoor) 
h_BSa=30; %BS antenna height (m)
h_BSb=30;  %BS antenna height (m)
h_BSc=30;
h_UEb=1.5; %UE antenna height (m)
G_0=15; %BS antenna gain (dBi)
d_sep=47000; % (m) seperation distance between both networks
n_a=1000; %number of monte carlo runs
offset=0; %degrees
r_a=1000; %cell radius (m)
alpha=0.8;
loading_factor=0.5;
aerial_factor=0.33;
%}

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
G_0=30; %BS antenna gain (dBi)
d_sep=11000; % (m) seperation distance between both networks
n_a=10000; %number of monte carlo runs
offset=60; %degrees
r_a=2500; %cell radius (m)
alpha=0.8;
loading_factor=0.5;
aerial_factor=0.33;

%{
%Mid Band URBAN
fc=1900e6; %center frequency (Hz)
lambda=299792458/fc; %wavelength (m)
bw=5e6; %bandwidth (Hz)
tilt=-10;  %mechanical downtilt (degrees)
L_feeder=-3; %(dB)
L_body=-4; %(dB)
SINR_target=25.3; %(dB)
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
offset=60; %degrees
r_a=250; %cell radius (m)
%}


OOBE=-8; %dBm/MHz
D=99.*lambda;
ACLR=-30;


%P_max=4.6e-6; %dBm/Hz max power
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

BS_r=[0,0,40];
BS_rV=[0,-1,0];


%count number of UEs over or under tx power range

n_UEmaxA=zeros(3,19);
n_UEminA=zeros(3,19);



%Precalculate Loss
%for d=1:300000
%P1546(d)=-P1546FieldStrMixed(fc./1000000,50,h_UEb,h_BSb,10,'Rural',d./1000,'Land',0, 'q', 50, 'Ptx', 1, 'ha', h_UEb);
%end
%Dat=load("P1546_708MHz_50%_1.5m_30m_10m_rural_1to3e5m.mat");
%P1546=Dat.P1546;

Dat=load("P452_10percent_Below22m.npy.mat");
P452=Dat.P452_10percent_Below22m;

temp=load("P528_0to100km_logspace_10to300m_1900MHz_0_50%");
L528=temp.L528;
d528=temp.d;
h528=temp.ht;


d_a=NaN(3,19,n_a,3);
phi_a=NaN(3,19,n_a);
theta_a=NaN(3,19,n_a);
AntGain_a_a=NaN(3,19,n_a);
L_a_a=NaN(3,19,n_a);
indoor_a=NaN(3,19,n_a);

I_dBmUEa=NaN(3,19,n_a);

%UE_a(3,:,:,:)=h_UEb;
%Loop to calculate tx power for all UEs using power control
for i=1:n_a %iterate over each monte carlo run
for n=1:19  %iterate over each BS
for s=1:3   %iterate over each sector

   
if rand<aerial_factor
UE_a(3,s,n,i)=rand*300; %generate random drone altitude
indoor_a(s,n,i)=0;
else
UE_a(3,s,n,i)=h_UEb;
  if rand<indoor_ratio
    indoor_a(s,n,i)=1;
  else
    indoor_a(s,n,i)=0;
  end
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
[~, Hindex] = min(abs(h528-UE_a(3,s,n,i)));
[~, Dindex] = min(abs(d528-vecnorm(d_a(s,n,i,1:2))));

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

if indoor_a(s,n,i)==0
I_dBmUEa(s,n,i)=min([max_dBmUE, M_pusch + P_0 - alpha.*(AntGain_a_a(s,n,i)+L_a_a(s,n,i)+UE_gain+L_feeder) ]);
elseif indoor_a(s,n,i)==1
I_dBmUEa(s,n,i)=min([max_dBmUE, M_pusch + P_0 - alpha.*(AntGain_a_a(s,n,i)+L_a_a(s,n,i)+UE_gain+L_feeder+L_entry) ]);
end


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

%I_dBmUEa(:,:,:)=OOBE+10*log10(bw./1e6);
I_dBmUEa=I_dBmUEa+ACLR;

% for d=1:100
 %   d_sep=d*1000 
%note there are two sets of displacement vectors and angles
%UPPER CASE is wrt to reference cell BS 
%lower case is wrt to own BS (calculated above already)


D_a=NaN(3,19,n_a,3);
PHI_a=NaN(3,19,n_a);
THETA_a=NaN(3,19,n_a);
AntGain_a_r=NaN(3,19,n_a);
L_a_r=NaN(3,19,n_a);


I_a=zeros(3,19,n_a); %Interference from aerial network



for n=1:19  %iterate over each BS in both networks
for s=1:3   %iterate over each sector
for i=1:n_a
    
    
%displacement vector to ref BS
D_a(s,n,i,:)=[UE_a(1,s,n,i)-BS_r(1),UE_a(2,s,n,i)-BS_r(2),UE_a(3,s,n,i)-BS_r(3)];

%angle between ref BS max gain vector and displacement vector in xy plane
PHI_a(s,n,i)=acosd(dot(squeeze(D_a(s,n,i,:)),BS_rV)./vecnorm(BS_rV)./vecnorm(squeeze(D_a(s,n,i,:))));

%elevation angle
THETA_a(s,n,i)=acosd(dot([vecnorm(D_a(s,n,i,1:2)),D_a(s,n,i,3)],[norm(BS_rV(1:2)),BS_rV(3)])./norm(BS_rV)./vecnorm(D_a(s,n,i,:))).*D_a(s,n,i,3)./abs(D_a(s,n,i,3));

%Antenna gain to ref BS
AntGain_a_r(s,n,i)=F1245(PHI_a(s,n,i,:),G_0,D,fc);

%propagation loss to ref BS
%free space
%L_a_r(S,N,s,n,i)=20.*log10(lambda./(4.*pi.*vecnorm(D_a(S,N,s,n,i,:))));

if UE_a(3,s,n,i)>22
%ITU P528
[~, Hindex] = min(abs(h528-UE_a(3,s,n,i)));
[~, Dindex] = min(abs(d528-vecnorm(D_a(s,n,i,1:2))));

    L_a_r(s,n,i)=-L528(Hindex,Dindex);
elseif UE_a(3,s,n,i)<22 && UE_a(3,s,n,i)>=19
    L_a_r(s,n,i)=-P452(round(vecnorm(D_a(s,n,i,1:2))+1),5);
elseif UE_a(3,s,n,i)<19 && UE_a(3,s,n,i)>=14
    L_a_r(s,n,i)=-P452(round(vecnorm(D_a(s,n,i,1:2))+1),4);
elseif UE_a(3,s,n,i)<14 && UE_a(3,s,n,i)>=9
    L_a_r(s,n,i)=-P452(round(vecnorm(D_a(s,n,i,1:2))+1),3);
elseif UE_a(3,s,n,i)<9 && UE_a(3,s,n,i)>=4
    L_a_r(s,n,i)=-P452(round(vecnorm(D_a(s,n,i,1:2))+1),2);
elseif UE_a(3,s,n,i)<4
    L_a_r(s,n,i)=-P452(round(vecnorm(D_a(s,n,i,1:2))+1),1);
end

if indoor_a(s,n,i)==0
I_a(s,n,i)=10.^((I_dBmUEa(s,n,i)+AntGain_a_r(s,n,i)+L_a_r(s,n,i)+UE_gain+L_feeder)./10); %mW
else
 I_a(s,n,i)=10.^((I_dBmUEa(s,n,i)+AntGain_a_r(s,n,i)+L_a_r(s,n,i)+UE_gain+L_feeder+L_entry)./10); %mW
end   
% 50% LOADING factor
if rand>loading_factor
   I_a(s,n,i)=0;
end

end

end
end

%Aggregate all interference (sum along n,s)
I=squeeze(sum(sum(I_a,1),2));

%I over Thermal
IoT_map=squeeze(mean(10*log10((I+N_thermal)./N_thermal),3));

%I over N
IoN=10*log10(I./N_thermal);

[CDF_IoN_P,CDF_IoN_V]=homemade_ecdf(IoN');
[CDF_I_P,CDF_I_V]=homemade_ecdf(10.*log10(I'));
plot(CDF_IoN_V,CDF_IoN_P)


%{

%Plotting hex outlines and UE positions
cmap=parula;
hexValues=log10(throughput_loss_map);
%hexValues=hexValues-min(hexValues);
hexColor=round((size(cmap,1)-1).*(hexValues-min(min(hexValues)))./(max(max(hexValues))-min(min(hexValues)))+1);

hold on
for n=1:19
for s=1:3
    
    
 plot(outline_c(:,1,n,s),outline_c(:,2,n,s),'g')   
 fill(outline_b(:,1,n,s),outline_b(:,2,n,s),cmap(hexColor(s,n),:))
 plot(outline_a(:,1,n,s),outline_a(:,2,n,s),'r')
 
 text(mean(outline_b(:,1,n,s)), mean(outline_b(:,2,n,s)) , num2str(throughput_loss_map(s,n)), 'HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'FontWeight', 'bold');

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
%}
toc
