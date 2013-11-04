function slit_diffraction()
% Need to compare this to diffraction integrals

L = 5;
incidentAmp = 1;
k = 10;

x = -k:k/40:k;
z = 0:k/40:k * 5;
[xMesh,zMesh] = meshgrid(x,z);

figure();
surface(xMesh,zMesh,abs(Amplitude()),'LineStyle', 'none');

    function amp = Amplitude()
        ampArray = integral(@(kx) Integrand(xMesh(:),zMesh(:),kx), -k*2, k*2, 'ArrayValued', true);
        amp = reshape(ampArray, size(xMesh));
    end

    function val = Integrand(x,z,kx)
        kz = (k.^2-kx.^2).^0.5;
        val = 1/(2*pi) * incidentAmp * L * sinc(kx * L/(2*pi)) .* exp(1i*(kx*x+kz*z));
    end

end