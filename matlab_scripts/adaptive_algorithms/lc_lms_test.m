clear;
close all;

N = 16; 
c = 343;

fs = 8000;

f = 2000;
phi_sig = 0;
teta_sig = 0;
phi_noise = 60;
teta_noise = 0;
c = 343;
d = c/(2*f);
d = 0.04;
L = 60;
N_f = 50; 
mu = 0.01;
p_loc = gen_place_el(4,4,d,d,1)';

[noise,fs] = audioread('nois.wav');
[signal,fs] = audioread('speech_dft_8kHz.wav');

t = 2;

signal = signal(1:t*fs);
noise = noise(1:t*fs);

signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);

sig_in_MR = signal_shift+noise_shift;
u = sig_in_MR(1,:);
p = sig_in_MR(4,:);
[y,W] = func_Frost(sig_in_MR,L,N,mu);
%%
figure()
hold on
plot(sig_in_MR(1,:))

plot(y)
plot(signal)
sound(y,fs)

[B,BB] = plot_bp_for_place(p_loc,phi_sig,teta_sig,N_f,W,fs,L,N);
% plot_wb_bp_freq_3d(W,N,L,p_loc,N_f,phi_sig,teta_sig);
