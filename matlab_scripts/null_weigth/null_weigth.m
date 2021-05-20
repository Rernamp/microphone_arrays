clear;
close all;

N = 16; 
J = 40;
c = 343;
N_f = 20;
fs = 8000;
d = 0.04;
phi_sig = 0;
teta_sig = 0;
phi_noise = 50;
teta_noise = 60;
N_f = J;
c = 343;
p_loc = gen_place_el(4,4,d,d,1)';


w_nm = ones(1,N*J);
f = [0:J-1 -J:-1]*fs/(2*J); 
f = f( 1 : length(f)/2 );
% f = 4000;
k_sig = 2*pi*f.*[-cosd(teta_sig).*cosd(phi_sig) ; -cosd(teta_sig).*sind(phi_sig) ; -sind(teta_sig)]/c;

k_noise = 2*pi*f.*[-cosd(teta_noise).*cosd(phi_noise) ; -cosd(teta_noise).*sind(phi_noise) ; -sind(teta_noise)]/c;

d_sig = exp(-1i*k_sig'*p_loc);
d_noise = exp(-1i*k_noise'*p_loc);
% d_sig(2:end,:) = zeros(39,16);
w = d_sig - (d_noise*d_sig')*d_noise/N; 
% w = d_sig -(d_sig'*d_noise*inv(d_noise*d_noise')*d_noise');
% w = [ones(1,16); zeros(39, 16)];
w_nm = reshape(w',1,[]);
%%
% w_nm = meshgrid(w,ones(1,J));
[B,BB] = plot_bp_for_place(p_loc,phi_sig,teta_sig,N_f,w_nm,fs,J,N);
