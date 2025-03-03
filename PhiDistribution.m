
r_a=6000;
n_a=1000000;
UE_a=zeros(3,n_a);
tilt=-3;
h=40;

%generate n_a random xy positions in cell
UE_a(1:2,:)=randHEX(r_a,0,0,n_a);

%Altitude is zero
UE_a(3,:)=0;
global HEXphi
global HEXtheta
global HEXd
for i=1:n_a
    
 D=[UE_a(1,i),UE_a(2,i)-r_a,-h]; %displacment vector
        HEXphi(i)=real(180.*acos(dot(D,[0,-1,sin(pi*tilt/180)])./norm(D))./pi);
        HEXtheta(i)=90+tilt-180.*atan(norm(D)./h)./pi;
        HEXd(i)=norm(D(1:2));
end



