% 20190528
% convert wavelength to RGB by linearly interpolating the 7 colors in
% rainbow
% input:
% wavelength: in unit of nm, range from 380 to 645
% see: http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm
% created by Pai Peng
function y=wavelength2RGB(Wavelength)
    gamma=0.8;
    if((Wavelength >= 380) && (Wavelength<440))
        Red = -(Wavelength - 440) / (440 - 380);
        Green = 0.0;
        Blue = 1.0;
    
    else
        if((Wavelength >= 440) && (Wavelength<490))
        Red = 0.0;
        Green = (Wavelength - 440) / (490 - 440);
        Blue = 1.0;
    else
        if((Wavelength >= 490) && (Wavelength<510))
        Red = 0.0;
        Green = 1.0;
        Blue = -(Wavelength - 510) / (510 - 490);
    else
        if((Wavelength >= 510) && (Wavelength<580))
        Red = (Wavelength - 510) / (580 - 510);
        Green = 1.0;
        Blue = 0.0;
    else
        if((Wavelength >= 580) && (Wavelength<=645))
        Red = 1.0;
        Green = -(Wavelength - 645) / (645 - 580);
        Blue = 0.0;
    else
        Red = 0.0;
        Green = 0.0;
        Blue = 0.0;
    end
    end
    end
    end
    end

    
if (Wavelength>=380) && (Wavelength<420)
    factor=0.3+0.7*(Wavelength-380)/(420-380);
else
    factor=1;
end
y=[Red,Green,Blue];
y=y*factor;
y=y.^gamma;