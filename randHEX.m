function [Points,outline] = randHEX(r,lat,lon,ntot)
%function [Points,outline] = randHEX(r,lat,lon,ntot)
%
%
%This function generates ntot random xy position in a hexagonal area
%it also retuns 7 coodinates that represent the edges of the hexagon
%
%the point of the hexagon is oriented vertically
%
%Points: 2 x ntot
%outline: 2 x 7
%
%Inputs
% 1x1 int, single, double
%
% r : hexagon radius, equivalent to length of side 
% lat : 'y' center position of hexagon
% lon : 'x' center position of hexagon
% ntot : number of random positions you wish to generate with one call of
% this function.
%
% NOTE
% r, lat, lon can be any unit so long as they're all the same



Points=zeros(2,ntot);


n=0;

while n<ntot
    
    %generate a random position within a rectangular area
    % 1 x 0.866 
    y=2*rand-1;
    x=sqrt(3)*rand-(sqrt(3)./2);
 
    %check if point falls within bounds of hexagon
    %if not; keep looping until there are ntot points within bounds.
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

%Rescale to chosen radius
Points=Points.*r;

%Translate to chosen centerposition
Points(1,:)=Points(1,:)+lon;
Points(2,:)=Points(2,:)+lat;

%generate outline (note firt and last coordinate are the same)
outline=[lat,lon]+[r,0;r/2,sqrt(3)*r/2;-r/2,sqrt(3)*r/2;-r,0;-r/2,-sqrt(3)*r/2;r/2,-sqrt(3)*r/2;r,0];
outline=flip(outline,2);
%plot(Points(1,:),Points(2,:),'.r');

end
