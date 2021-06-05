
clear;close all;


K = 9; 
c = 343;

fs = 8000;

f = 2000;
phi_sig = 0;
teta_sig = 0;
phi_noise_1 = 50;
teta_noise_1 = 0;
phi_noise_2 = -50;
teta_noise_2 = 0;
c = 343;
d = c/(2*f);
d = 0.04;
L = 32;
N_f = 50; 
mu = 0.01;
c_noise =35;
p_loc = gen_place_el(sqrt(K),sqrt(K),d,d,1)';

[noise,fs] = audioread('nois.wav');
[signal,fs] = audioread('speech_dft_8kHz.wav');
osh_in = mean(signal.^2)/mean(noise.^2);
t = 2;
time = 0:1/fs:t-1/fs;
signal = signal(1:t*fs);
noise = noise(1:t*fs);

signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
noise_shift_1 = shift_plane(noise,phi_noise_1,teta_noise_1,p_loc,fs);
noise_shift_2 = shift_plane(noise,phi_noise_2,teta_noise_2,p_loc,fs);
noise_shift = [noise_shift_1 noise_shift_2];
signal_shift = [signal_shift signal_shift];
sig_in_MR = signal_shift + noise_shift;


time = 0:1/fs:(length(sig_in_MR(1,:))-1)/fs;

%%

[y,W_n] = func_LC_RLS(sig_in_MR,L,K);
figure()
hold on
plot(time,y)
xline(time(length(time)/2))
grid on
ylabel("Амплитуда, усл.ед")
xlabel("Время, с")




figure()
hold on
plot(time, sig_in_MR(1,:))
xline(time(length(time)/2))
grid on
ylabel("Амплитуда, усл.ед")
xlabel("Время, с")

%%
phi_const = 0;
teta_const = 0;
[B,BB] = plot_bp_for_place(p_loc,phi_const,teta_const,N_f,W_n,fs,L,K);
