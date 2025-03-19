function [HEXphi,HEXtheta,HEXd] = PhiDistribution(h_UE,h_BS,r_a,n_a)


UE_a=NaN(2,n_a);
%generate n_a random xy positions in cell
UE_a(1:2,:)=randHEX(r_a,0,0,n_a);

%Altitude is zero
UE_a(3,:)=h_UE;

HEXphi=NaN(1,n_a);
HEXtheta=NaN(1,n_a);
HEXd=NaN(1,n_a);

for i=1:n_a
    
 D=[UE_a(1,i),UE_a(2,i)-r_a,h_UE-h_BS]; %displacment vector
        HEXphi(i)=acosd(dot(D(1:2),[0,-1])./norm(D(1:2)));
        HEXtheta(i)=acosd(dot([norm(D(1:2)),D(3)],[1,0])./norm(D));      
        HEXd(i)=norm(D(1:2));
end


end
