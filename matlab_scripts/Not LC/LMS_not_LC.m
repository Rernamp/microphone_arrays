clear;
close all;

N = 16; 
c = 343;

fs = 8000;

f = 2000;
phi_sig = 0;
teta_sig = 0;
phi_noise = 60;
teta_noise = 30;
c = 343;
J = 16;
% d = c/(2*f);
t = 2;
time = 0:1/fs:t-1/fs;
d = 0.04;
N_f = 50; 
mu = 0.1;
del_p = 0.005;
p_loc = gen_place_el(sqrt(N),sqrt(N),d,d,1)';

[noise,fs] = audioread('nois.wav');
[signal,fs] = audioread('sound.wav');

t = 2;
a_sig = [-cosd(teta_sig).*cosd(phi_sig)  -cosd(teta_sig).*sind(phi_sig)  -sind(teta_sig)];

tau_sig = a_sig*p_loc/c;

a_noise = [-cosd(teta_noise).*cosd(phi_noise)  -cosd(teta_noise).*sind(phi_noise)  -sind(teta_noise)];

tau_noise = a_noise*p_loc/c;



signal = signal(1:t*fs);
noise = noise(1:t*fs);

signal_shift = shift_plane_tau(signal,fs,tau_sig,p_loc);
noise_shift = shift_plane_tau(noise,fs,tau_noise,p_loc);

sig_in_MR = signal_shift+noise_shift;
mu = 1;

[y,w] = func_LMS(sig_in_MR, J, N, mu );