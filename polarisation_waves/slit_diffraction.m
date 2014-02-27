function slit_diffraction()
% Need to compare this to diffraction integrals

L = 4;
incidentAmp = 1;
k = 10;

x = -k:k/80:k;
z = 0:k/80:k * 5;
[xMesh,zMesh] = meshgrid(x,z);

figure();
surface(xMesh,zMesh,abs(Amplitude()),'LineStyle', 'none');

    function amp = Amplitude()
        ampArray = integral(@(kx) Integrand(xMesh(:),zMesh(:),kx), -k, k, 'ArrayValued', true);
        amp = reshape(ampArray, size(xMesh));
    end

    function val = Integrand(x,z,kx)
        kz = (k.^2-kx.^2).^0.5;
        phase = kx*x + kz*z;
        const = 1/(2*pi) * incidentAmp * L;
        val = const * sinc(kx * L/(2*pi)) .* exp(1i*phase);
    end

end