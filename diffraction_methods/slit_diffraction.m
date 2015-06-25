function slit_diffraction()
% Need to compare this to diffraction integrals

slit_width = 4;
incident_amplitude = 1;
space_width = 10;

x = -space_width:space_width/80:space_width;
z = 0:space_width/80:space_width * 5;
[xMesh,zMesh] = meshgrid(x,z);

figure();
surface(xMesh,zMesh,abs(Amplitude()),'LineStyle', 'none');

    function amp = Amplitude()
        ampArray = integral(@(kx) Integrand(xMesh(:),zMesh(:),kx), -space_width, space_width, 'ArrayValued', true);
        amp = reshape(ampArray, size(xMesh));
    end

    function val = Integrand(x,z,kx)
        kz = (space_width.^2-kx.^2).^0.5;
        phase = kx*x + kz*z;
        const = 1/(2*pi) * incident_amplitude * slit_width;
        val = const * sinc(kx * slit_width/(2*pi)) .* exp(1i*phase);
    end

end