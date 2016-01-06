function w = quartic1(C,D,F0,L,V,Z0,lambda,v)
w = [(F0.*lambda)./V;...
    (((1.0./2.0)./Z0+L.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*1.0./(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(3.0./2.0)+(v.*(1.0./2.0))./(V.*Z0)-F0.^2.*L.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*1.0./(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(5.0./2.0).*3.0-(C.*F0.*L.*1.0./V.^2.*lambda.*1.0./(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(3.0./2.0).*(1.0./2.0))./Z0)./(L.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*1.0./sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0)-F0.^2.*L.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*1.0./(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(3.0./2.0)+1.0)+((F0.*lambda)./V+(F0.*L.*lambda.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*1.0./(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(3.0./2.0))./V).*1.0./(L.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*1.0./sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0)-F0.^2.*L.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*1.0./(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(3.0./2.0)+1.0).^2.*((C.*L.*1.0./sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(1.0./2.0))./(V.*Z0)-F0.^3.*L.*1.0./V.^3.*lambda.^3.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*1.0./(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(5.0./2.0).*3.0+(F0.*L.*lambda.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*1.0./(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(3.0./2.0).*3.0)./V-(C.*F0.^2.*L.*1.0./V.^3.*lambda.^2.*1.0./(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(3.0./2.0).*(1.0./2.0))./Z0))./(L.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*1.0./sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0)-F0.^2.*L.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*1.0./(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(3.0./2.0)+1.0);...
    -((F0.^2.*1.0./V.^2.*lambda.^2+1.0).^2.*((C.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(1.0./2.0))./(V.*Z0)+(C.*L.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)))./(V.*Z0)-(F0.*L.*lambda.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*3.0)./V-(F0.*L.^3.*lambda.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*3.0)./V-F0.^3.*L.*1.0./V.^3.*lambda.^3.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*9.0-F0.^5.*L.*1.0./V.^5.*lambda.^5.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*9.0-F0.^7.*L.*1.0./V.^7.*lambda.^7.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*3.0+(C.*L.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(1.0./2.0))./(V.*Z0)-F0.^3.*L.^2.*1.0./V.^3.*lambda.^3.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*6.0+(C.*F0.^2.*1.0./V.^3.*lambda.^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(5.0./2.0))./Z0+(C.*F0.^4.*1.0./V.^5.*lambda.^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*5.0)./Z0+(C.*F0.^6.*1.0./V.^7.*lambda.^6.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*5.0)./Z0+(C.*F0.^8.*1.0./V.^9.*lambda.^8.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(5.0./2.0))./Z0+(C.*F0.^10.*1.0./V.^11.*lambda.^10.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(1.0./2.0))./Z0-(F0.*L.^2.*lambda.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*6.0)./V+(C.*F0.^2.*L.*1.0./V.^3.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*4.0)./Z0+(C.*F0.^4.*L.*1.0./V.^5.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*6.0)./Z0+(C.*F0.^6.*L.*1.0./V.^7.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*4.0)./Z0+(C.*F0.^8.*L.*1.0./V.^9.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)))./Z0+(C.*F0.^2.*L.^2.*1.0./V.^3.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0))./Z0+(C.*F0.^4.*L.^2.*1.0./V.^5.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(1.0./2.0))./Z0))./(L.^3.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*1.0e1+L.^5.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5+(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(1.5e1./2.0)+L.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*5.0+L.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(1.3e1./2.0).*1.0e1+L.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(1.1e1./2.0).*5.0+F0.^2.*L.^3.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*3.0e1+F0.^4.*L.^3.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*3.0e1+F0.^6.*L.^3.*1.0./V.^6.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*1.0e1+F0.^2.*L.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*3.0e1+F0.^4.*L.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*7.5e1+F0.^6.*L.*1.0./V.^6.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*1.0e2+F0.^8.*L.*1.0./V.^8.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*7.5e1+F0.^10.*L.*1.0./V.^10.*lambda.^10.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*3.0e1+F0.^12.*L.*1.0./V.^12.*lambda.^12.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*5.0-F0.^2.*L.^2.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(1.1e1./2.0).*2.0e1-F0.^2.*L.^4.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(9.0./2.0).*2.0e1+F0.^4.*L.^2.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(9.0./2.0).*1.0e1+F0.^4.*L.^4.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(7.0./2.0).*3.0e1-F0.^6.*L.^4.*1.0./V.^6.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(5.0./2.0).*2.0e1+F0.^8.*L.^4.*1.0./V.^8.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*(F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(3.0./2.0).*5.0);...
    ((F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(5.0./2.0).*1.0./(-L.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2+F0.^2.*1.0./V.^2.*lambda.^2.*3.0+F0.^4.*1.0./V.^4.*lambda.^4.*3.0+F0.^6.*1.0./V.^6.*lambda.^6+1.0).^4.*(D.*1.2e1-L.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*9.0+L.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^7.*6.0+L.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^9.*3.0+L.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.0+D.*L.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*2.4e1-D.*L.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*3.6e1+L.^3.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^6.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*6.0-L.^5.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^8.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*9.0+D.*L.^3.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*2.4e1+D.*L.^5.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*1.2e1+D.*F0.^2.*1.0./V.^2.*lambda.^2.*1.32e2+D.*F0.^4.*1.0./V.^4.*lambda.^4.*6.6e2+D.*F0.^6.*1.0./V.^6.*lambda.^6.*1.98e3+D.*F0.^8.*1.0./V.^8.*lambda.^8.*3.96e3+D.*F0.^10.*1.0./V.^10.*lambda.^10.*5.544e3+D.*F0.^12.*1.0./V.^12.*lambda.^12.*5.544e3+D.*F0.^14.*1.0./V.^14.*lambda.^14.*3.96e3+D.*F0.^16.*1.0./V.^16.*lambda.^16.*1.98e3+D.*F0.^18.*1.0./V.^18.*lambda.^18.*6.6e2+D.*F0.^20.*1.0./V.^20.*lambda.^20.*1.32e2+D.*F0.^22.*1.0./V.^22.*lambda.^22.*1.2e1-D.*L.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.6e1-C.^2.*L.*1.0./V.^2.*1.0./Z0.^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(3.0./4.0)+C.^2.*L.^4.*1.0./V.^2.*1.0./Z0.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*3.0+F0.^2.*L.^2.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*9.0+F0.^2.*L.^4.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^7.*1.56e2+F0.^4.*L.^2.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*2.43e2+F0.^2.*L.^6.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^9.*1.5e1+F0.^4.*L.^4.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^7.*4.32e2+F0.^6.*L.^2.*1.0./V.^6.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*7.65e2+F0.^6.*L.^4.*1.0./V.^6.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^7.*4.2e2+F0.^8.*L.^2.*1.0./V.^8.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*1.125e3+F0.^8.*L.^4.*1.0./V.^8.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^7.*1.38e2+F0.^10.*L.^2.*1.0./V.^10.*lambda.^10.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*8.91e2+F0.^12.*L.^2.*1.0./V.^12.*lambda.^12.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*3.69e2+F0.^14.*L.^2.*1.0./V.^14.*lambda.^14.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*6.3e1+C.^2.*L.^2.*1.0./V.^2.*1.0./Z0.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*3.0+F0.^2.*L.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*9.0-F0.^4.*L.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*2.1e1-F0.^6.*L.*1.0./V.^6.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*1.47e2-F0.^8.*L.*1.0./V.^8.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.15e2-F0.^10.*L.*1.0./V.^10.*lambda.^10.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.57e2-F0.^12.*L.*1.0./V.^12.*lambda.^12.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*2.31e2-F0.^14.*L.*1.0./V.^14.*lambda.^14.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*8.1e1-F0.^16.*L.*1.0./V.^16.*lambda.^16.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*1.2e1-C.^2.*L.^3.*1.0./V.^2.*1.0./Z0.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(9.0./2.0)-C.^2.*L.^5.*1.0./V.^2.*1.0./Z0.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(3.0./4.0)+D.*F0.^2.*L.^2.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*1.92e2-D.*F0.^2.*L.^4.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*1.8e2+D.*F0.^4.*L.^2.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*6.72e2-D.*F0.^4.*L.^4.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*3.6e2+D.*F0.^6.*L.^2.*1.0./V.^6.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*1.344e3-D.*F0.^6.*L.^4.*1.0./V.^6.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*3.6e2+D.*F0.^8.*L.^2.*1.0./V.^8.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*1.68e3-D.*F0.^8.*L.^4.*1.0./V.^8.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*1.8e2+D.*F0.^10.*L.^2.*1.0./V.^10.*lambda.^10.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*1.344e3-D.*F0.^10.*L.^4.*1.0./V.^10.*lambda.^10.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*3.6e1+D.*F0.^12.*L.^2.*1.0./V.^12.*lambda.^12.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*6.72e2+D.*F0.^14.*L.^2.*1.0./V.^14.*lambda.^14.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*1.92e2+D.*F0.^16.*L.^2.*1.0./V.^16.*lambda.^16.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*2.4e1-F0.^2.*L.^3.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^6.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*1.08e2-F0.^2.*L.^5.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^8.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*8.1e1-F0.^4.*L.^3.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^6.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*4.92e2-F0.^4.*L.^5.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^8.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*7.2e1-F0.^6.*L.^3.*1.0./V.^6.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^6.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*7.68e2-F0.^8.*L.^3.*1.0./V.^8.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^6.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*5.22e2-F0.^10.*L.^3.*1.0./V.^10.*lambda.^10.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^6.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*1.32e2-D.*F0.^2.*L.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.24e2-D.*F0.^4.*L.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*1.296e3-D.*F0.^6.*L.*1.0./V.^6.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.024e3-D.*F0.^8.*L.*1.0./V.^8.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*4.536e3-D.*F0.^10.*L.*1.0./V.^10.*lambda.^10.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*4.536e3-D.*F0.^12.*L.*1.0./V.^12.*lambda.^12.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.024e3-D.*F0.^14.*L.*1.0./V.^14.*lambda.^14.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*1.296e3-D.*F0.^16.*L.*1.0./V.^16.*lambda.^16.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.24e2-D.*F0.^18.*L.*1.0./V.^18.*lambda.^18.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.6e1+C.^2.*F0.^2.*L.^2.*1.0./V.^4.*1.0./Z0.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*2.4e1+(C.*F0.^3.*L.^2.*1.0./V.^4.*lambda.^3.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*2.52e2)./Z0+(C.*F0.^3.*L.^4.*1.0./V.^4.*lambda.^3.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*1.44e2)./Z0+C.^2.*F0.^4.*L.^2.*1.0./V.^6.*1.0./Z0.^2.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*8.4e1+(C.*F0.^5.*L.^2.*1.0./V.^6.*lambda.^5.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*7.56e2)./Z0+(C.*F0.^5.*L.^4.*1.0./V.^6.*lambda.^5.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*2.16e2)./Z0+C.^2.*F0.^6.*L.^2.*1.0./V.^8.*1.0./Z0.^2.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*1.68e2+(C.*F0.^7.*L.^2.*1.0./V.^8.*lambda.^7.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*1.26e3)./Z0+(C.*F0.^7.*L.^4.*1.0./V.^8.*lambda.^7.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*1.44e2)./Z0+C.^2.*F0.^8.*L.^2.*1.0./V.^10.*1.0./Z0.^2.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*2.1e2+(C.*F0.^9.*L.^2.*1.0./V.^10.*lambda.^9.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*1.26e3)./Z0+(C.*F0.^9.*L.^4.*1.0./V.^10.*lambda.^9.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*3.6e1)./Z0+C.^2.*F0.^10.*L.^2.*1.0./V.^12.*1.0./Z0.^2.*lambda.^10.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*1.68e2+(C.*F0.^11.*L.^2.*1.0./V.^12.*lambda.^11.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*7.56e2)./Z0+C.^2.*F0.^12.*L.^2.*1.0./V.^14.*1.0./Z0.^2.*lambda.^12.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*8.4e1+(C.*F0.^13.*L.^2.*1.0./V.^14.*lambda.^13.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*2.52e2)./Z0+C.^2.*F0.^14.*L.^2.*1.0./V.^16.*1.0./Z0.^2.*lambda.^14.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*2.4e1+(C.*F0.^15.*L.^2.*1.0./V.^16.*lambda.^15.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*3.6e1)./Z0+C.^2.*F0.^16.*L.^2.*1.0./V.^18.*1.0./Z0.^2.*lambda.^16.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).*3.0-C.^2.*F0.^2.*L.*1.0./V.^4.*1.0./Z0.^2.*lambda.^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(2.7e1./4.0)-C.^2.*F0.^4.*L.*1.0./V.^6.*1.0./Z0.^2.*lambda.^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*2.7e1-C.^2.*F0.^6.*L.*1.0./V.^8.*1.0./Z0.^2.*lambda.^6.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*6.3e1-C.^2.*F0.^8.*L.*1.0./V.^10.*1.0./Z0.^2.*lambda.^8.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(1.89e2./2.0)-C.^2.*F0.^10.*L.*1.0./V.^12.*1.0./Z0.^2.*lambda.^10.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(1.89e2./2.0)-C.^2.*F0.^12.*L.*1.0./V.^14.*1.0./Z0.^2.*lambda.^12.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*6.3e1-C.^2.*F0.^14.*L.*1.0./V.^16.*1.0./Z0.^2.*lambda.^14.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*2.7e1-C.^2.*F0.^16.*L.*1.0./V.^18.*1.0./Z0.^2.*lambda.^16.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(2.7e1./4.0)-C.^2.*F0.^18.*L.*1.0./V.^20.*1.0./Z0.^2.*lambda.^18.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(3.0./4.0)+(C.*F0.*L.^2.*1.0./V.^2.*lambda.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*3.6e1)./Z0+(C.*F0.*L.^4.*1.0./V.^2.*lambda.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*3.6e1)./Z0+C.^2.*F0.^2.*L.^4.*1.0./V.^4.*1.0./Z0.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*1.5e1+C.^2.*F0.^4.*L.^4.*1.0./V.^6.*1.0./Z0.^2.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*3.0e1+C.^2.*F0.^6.*L.^4.*1.0./V.^8.*1.0./Z0.^2.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*3.0e1+C.^2.*F0.^8.*L.^4.*1.0./V.^10.*1.0./Z0.^2.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*1.5e1+C.^2.*F0.^10.*L.^4.*1.0./V.^12.*1.0./Z0.^2.*lambda.^10.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*3.0+D.*F0.^2.*L.^3.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*1.44e2+D.*F0.^2.*L.^5.*1.0./V.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.6e1+D.*F0.^4.*L.^3.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.6e2+D.*F0.^4.*L.^5.*1.0./V.^4.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.6e1+D.*F0.^6.*L.^3.*1.0./V.^6.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*4.8e2+D.*F0.^6.*L.^5.*1.0./V.^6.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^5.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*1.2e1+D.*F0.^8.*L.^3.*1.0./V.^8.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*3.6e2+D.*F0.^10.*L.^3.*1.0./V.^10.*lambda.^10.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*1.44e2+D.*F0.^12.*L.^3.*1.0./V.^12.*lambda.^12.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^3.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*2.4e1-(C.*F0.^3.*L.*1.0./V.^4.*lambda.^3.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*7.2e1)./Z0-(C.*F0.^5.*L.*1.0./V.^6.*lambda.^5.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*2.52e2)./Z0-(C.*F0.^7.*L.*1.0./V.^8.*lambda.^7.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*5.04e2)./Z0-(C.*F0.^9.*L.*1.0./V.^10.*lambda.^9.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*6.3e2)./Z0-(C.*F0.^11.*L.*1.0./V.^12.*lambda.^11.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*5.04e2)./Z0-(C.*F0.^13.*L.*1.0./V.^14.*lambda.^13.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*2.52e2)./Z0-(C.*F0.^15.*L.*1.0./V.^16.*lambda.^15.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*7.2e1)./Z0-(C.*F0.^17.*L.*1.0./V.^18.*lambda.^17.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*9.0)./Z0-(C.*F0.*L.*1.0./V.^2.*lambda.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*9.0)./Z0-(C.*F0.^3.*L.^3.*1.0./V.^4.*lambda.^3.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*2.7e2)./Z0-(C.*F0.^3.*L.^5.*1.0./V.^4.*lambda.^3.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^6.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*1.8e1)./Z0-(C.*F0.^5.*L.^3.*1.0./V.^6.*lambda.^5.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*5.4e2)./Z0-(C.*F0.^5.*L.^5.*1.0./V.^6.*lambda.^5.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^6.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*9.0)./Z0-(C.*F0.^7.*L.^3.*1.0./V.^8.*lambda.^7.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*5.4e2)./Z0-(C.*F0.^9.*L.^3.*1.0./V.^10.*lambda.^9.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*2.7e2)./Z0-(C.*F0.^11.*L.^3.*1.0./V.^12.*lambda.^11.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*5.4e1)./Z0-(C.*F0.*L.^3.*1.0./V.^2.*lambda.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*5.4e1)./Z0-(C.*F0.*L.^5.*1.0./V.^2.*lambda.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^6.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*9.0)./Z0-C.^2.*F0.^2.*L.^3.*1.0./V.^4.*1.0./Z0.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*2.7e1-C.^2.*F0.^2.*L.^5.*1.0./V.^4.*1.0./Z0.^2.*lambda.^2.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(9.0./4.0)-C.^2.*F0.^4.*L.^3.*1.0./V.^6.*1.0./Z0.^2.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(1.35e2./2.0)-C.^2.*F0.^4.*L.^5.*1.0./V.^6.*1.0./Z0.^2.*lambda.^4.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(9.0./4.0)-C.^2.*F0.^6.*L.^3.*1.0./V.^8.*1.0./Z0.^2.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*9.0e1-C.^2.*F0.^6.*L.^5.*1.0./V.^8.*1.0./Z0.^2.*lambda.^6.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^4.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(3.0./4.0)-C.^2.*F0.^8.*L.^3.*1.0./V.^10.*1.0./Z0.^2.*lambda.^8.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(1.35e2./2.0)-C.^2.*F0.^10.*L.^3.*1.0./V.^12.*1.0./Z0.^2.*lambda.^10.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*2.7e1-C.^2.*F0.^12.*L.^3.*1.0./V.^14.*1.0./Z0.^2.*lambda.^12.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)).^2.*sqrt(F0.^2.*1.0./V.^2.*lambda.^2+1.0).*(9.0./2.0)))./((F0.^2.*1.0./V.^2.*lambda.^2+1.0).^(3.0./2.0)+L.*((1.0./2.0)./Z0+(v.*(1.0./2.0))./(V.*Z0)))];
end