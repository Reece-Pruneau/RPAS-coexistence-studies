
Glenn_H=dlmread('F1336H_Pattern.csv');
%Glenn_V=dlmread('F1336V_Pattern.xlsx');


Glenn_H(182:360,1)=Glenn_H(182:360,1)-360;
Glenn_H(:,1)=pi*Glenn_H(:,1)./180;

phi=linspace(-90,90,360);
R=linspace(0,0,360);
rad=phi.*pi./180;

for t=1:length(phi)
%R(t)= AntGain(phi(t),30,"APEREC028V01");
R(t)=F1336(0,phi(t),18);
end

%hold on
polarplot(rad,R)
%rlim([min(R),max(R)]);
%polarplot(Glenn_H(:,1),Glenn_H(:,2))
%rlim([min(Glenn_H(:,2)),max(Glenn_H(:,2))]);
rlim([min(R),max(R)]);
%hold off