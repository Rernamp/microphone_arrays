clear;
close all;

N = 16; 
c = 343;

fs = 8000;
d = 0.04;
f = 2000;
phi_sig = 0;
teta_sig = 0;
phi_noise = 0;
teta_noise = 90;
c = 343;
d = c/(8*f);
d = 0.04;
p_loc = gen_place_el(4,4,d,d,1)';

k_sig = 2*pi*f.*[-cosd(teta_sig).*cosd(phi_sig) ; -cosd(teta_sig).*sind(phi_sig) ; -sind(teta_sig)]/c;

k_noise = 2*pi*f.*[-cosd(teta_noise).*cosd(phi_noise) ; -cosd(teta_noise).*sind(phi_noise) ; -sind(teta_noise)]/c;

d_sig = exp(-1i*k_sig'*p_loc);
d_noise = exp(-1i*k_noise'*p_loc);
w = d_sig - (d_noise*d_sig')*d_noise/N; 
w = d_sig;
w = ones(1,N);
% w = [1 zeros(1,N-1)];
phi = 0:0.5:360;
teta = -90:0.5:90;
[PHI, TETA] = meshgrid(phi,teta);

BP = zeros(length(PHI(:,1)),length(PHI(1,:)));
i =1;

for phi_i = 0:0.5:360
    a_phi = [-cosd(teta).*cosd(phi_i) ; -cosd(teta).*sind(phi_i) ; -sind(teta)];
    k = 2*pi*f*a_phi/c;
    
    BP(:,i) = exp(-1i*k'*p_loc)*w';
    i = i +1;
end

BP = abs(BP).^2; 

    
BP_dB = 10*log10(BP);
% BP_dB(BP_dB < -40) = -40;

PHI = pi*PHI/180;
TETA = pi*TETA/180;

[X,Y,BP_dec] = sph2cart(PHI,TETA,BP_dB);

surf(X,Y,BP_dec)
shading interp 
grid on %
colormap gray
xlabel('X');
ylabel('Y');

% [B,BB] = plot_bp_for_place(p_loc,phi_sig,teta_sig,N_f,w_nm,fs,J,N);
