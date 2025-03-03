function [Points,outline] = randHEX(r,lat,lon,ntot)


Points=zeros(2,ntot);


n=0;

while n<ntot
    
    y=2*rand-1;
    x=sqrt(3)*rand-(sqrt(3)./2);
 
    if x<=0 && y>=0
        if y<=(x/sqrt(3)+1)
            n=n+1;
            Points(1,n)=x;
            Points(2,n)=y;
        end
    elseif x>0 && y>0
        if y<=(-x/sqrt(3)+1)
            n=n+1;
            Points(1,n)=x;
            Points(2,n)=y;  
        end
    elseif x<0 && y<0
        if y>=-(x/sqrt(3)+1)
            n=n+1;
            Points(1,n)=x;
            Points(2,n)=y;
        end
    elseif x>0 && y<0
        if y>=-(-x/sqrt(3)+1)
            n=n+1;
            Points(1,n)=x;
            Points(2,n)=y;
        end
    end

end 

Points=Points.*r;
Points(1,:)=Points(1,:)+lon;
Points(2,:)=Points(2,:)+lat;

outline=[lat,lon]+[r,0;r/2,sqrt(3)*r/2;-r/2,sqrt(3)*r/2;-r,0;-r/2,-sqrt(3)*r/2;r/2,-sqrt(3)*r/2;r,0];
%plot(Points(1,:),Points(2,:),'.r');

end
