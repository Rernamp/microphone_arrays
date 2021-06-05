clear;
close all;

N = 16; 
c = 343;

fs = 8000;

f = 2000;
phi_sig = 0;
teta_sig = 0;
phi_noise = 50;
teta_noise = 0;
c = 343;
d = c/(2*f);
d = 0.04;
L = 32;
N_f = 50; 
mu = 1;
c_noise =35;
p_loc = gen_place_el(4,4,d,d,1)';

[noise,fs] = audioread('nois.wav');
[signal,fs] = audioread('speech_dft_8kHz.wav');

t = 2;
time = 0:1/fs:t-1/fs;
signal = signal(1:t*fs);
noise = noise(1:t*fs);

signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);

sig_in_MR = signal_shift+noise_shift;
% for i = 1:length(sig_in_MR(:,1))
%     sig_in_MR(i,:) = awgn(sig_in_MR(i,:),c_noise+rand);
% end

[y_LMS,W_LMS] = func_LC_NLMS(sig_in_MR,L,N,mu);
[y_RLS,W_RLS] = func_LC_RLS(sig_in_MR,L,N);
%%

close all
figure()
hold on
grid on
plot(time,sig_in_MR(1,:))

ylabel("Амплитуда, усл.ед")
xlabel("Время, с")


figure()
hold on
grid on
plot(time,y_LMS)

ylabel("Амплитуда, усл.ед")
xlabel("Время, с")
%%
phi_const = 0;
teta_const = 0;
[B,BB] = plot_bp_for_place(p_loc,phi_const,teta_const,N_f,W_LMS,fs,L,N);
%%

close all
figure()
hold on
grid on
plot(time,sig_in_MR(1,:))

ylabel("Амплитуда, усл.ед")
xlabel("Время, с")
%%

figure()
hold on
grid on
plot(time,y_RLS)

ylabel("Амплитуда, усл.ед")
xlabel("Время, с")
%%
phi_const = 0;
teta_const = 0;
[B,BB] = plot_bp_for_place(p_loc,phi_const,teta_const,N_f,W_RLS,fs,L,N);
