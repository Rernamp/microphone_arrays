clear;
close all;

N = 9; 
c = 343;

fs = 8000;

f = 2000;
phi_sig = 0;
teta_sig = 0;
phi_noise = 30;
teta_noise = 30;
c = 343;
J = 54;
% d = c/(2*f);
t = 2;
time = 0:1/fs:t-1/fs;
d = 0.04;
N_f = 50; 
mu = 0.1;
p_loc = gen_place_el(sqrt(N),sqrt(N),d,d,1)';

[noise,fs] = audioread('nois.wav');
[signal,fs] = audioread('speech_dft_8kHz.wav');

t = 2;

signal = signal(1:t*fs);
noise = noise(1:t*fs);

signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);

sig_in_MR = signal_shift+noise_shift;
%% algoritm
y_GSC_RLS = GSC_RLS(sig_in_MR,J,N);
y_RLS = spat_filt_wb_time_lc_rls(sig_in_MR,J,N);

%%
figure()
hold on
grid on
% plot(time,sig_in_MR(1,:))

plot(time,y_GSC_RLS)
plot(time,y_RLS)
% legend('input','LC NLMS','GSC LMS')
legend('input','GSC RLS')
