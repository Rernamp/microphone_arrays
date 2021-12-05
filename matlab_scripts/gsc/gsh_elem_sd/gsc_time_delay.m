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
% p_loc(2,3) = p_loc(2,3) + del_p;
% p_loc(3,3) = p_loc(3,3) - del_p;
% 
% p_loc(2,6) = p_loc(2,6) - del_p;
% p_loc(3,6) = p_loc(3,6) + del_p;
[noise,fs] = audioread('nois.wav');
[signal,fs] = audioread('sound.wav');

t = 2;
a_sig = [-cosd(teta_sig).*cosd(phi_sig)  -cosd(teta_sig).*sind(phi_sig)  -sind(teta_sig)];

tau_sig = a_sig*p_loc/c;

a_noise = [-cosd(teta_noise).*cosd(phi_noise)  -cosd(teta_noise).*sind(phi_noise)  -sind(teta_noise)];

tau_noise = a_noise*p_loc/c;

del_tau = 2.6667e-7;

tau_noise(1) = tau_noise(1) + del_tau;
tau_sig(1) = tau_sig(1) + del_tau;

signal = signal(1:t*fs);
noise = noise(1:t*fs);

signal_shift = shift_plane_tau(signal,fs,tau_sig,p_loc);
noise_shift = shift_plane_tau(noise,fs,tau_noise,p_loc);

sig_in_MR = signal_shift+noise_shift;
%% algoritm
y_GSC = GSC_RLS(sig_in_MR,J,N);
% [y_LMS,W_LMS] = spat_filt_wb_time_lc_lms(sig_in_MR,J,N,mu);
%%
figure()
hold on
grid on
% plot(time,signal)
plot(time,abs(fft(signal)))
plot(time,abs(fft(y_GSC)))

% plot(time,y_LMS)
% plot(time,y_GSC)
legend('input','LC NLMS','GSC LMS')
% legend('input','GSC LMS')
