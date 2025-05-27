function G = F1245(phi,G0,D,fc)

%function G=F1245(phi,G0,D,fc)
%
%function returns the gain of a sectoral antenna in dBi
%Based on the ITU model F.1336 recomendation 3.1
%INPUTS
%phi (degrees) horizontal plane azimuth angle
%theta (degrees) vertical plane elevation angle
%
%
%All other parameters must be manually set bellow
%
lambda=299792458./fc;

phi_r=39.8.*(D./lambda).^-0.8;
G_a=G0-(2.5e-3*(D.*phi./lambda).^2);
G_b=2+15*log10(D./lambda);
F=10*log10(0.9*sind(3*pi*phi./(2*phi_r)).^2+0.1);

if D/lambda<=100
    if fc>=1e9 && fc<=70e9
        
        if phi>=0 && phi<phi_r
            
            G=max(G_a,G_b);
            
        elseif phi>=phi_r && phi<48
            
            G=42-5*log10(D./lambda)-25.*log10(phi)+F;
            
        elseif phi>=48 && phi<=180
            
            G=-5*log10(D/lambda)+F;
            
        else
            error('Input angle must be between 0 and 180 degrees')
        end
        
    else
        error('Frequency must be between 1 GHz and 70 GHz. Input as Hertz')
    end
else
    error('antenna Diameter must be less that 100 times carier wavelength')
end






end