clear;
close all;

N = 16; 
c = 343;
L = 32;
fs = 8000;
f = 2000;
phi_sig = 0;
teta_sig = 0;
phi_noise = 0;
teta_noise = 55;
c = 343;
d = c/(2*f);
d = 0.04;
p_loc = gen_place_el(4,4,d,d,1)';

[noise,fs] = audioread('nois.wav');

[sig,fs] = audioread('speech_dft_8kHz.wav');

t = 2;
m = 1000;
time = 0:1/fs:(length(noise)-1)/fs;
%%
sig = sig(1:t*fs);
noise = noise(1:t*fs);
SNR = -40:2:40;
mu = 0.5;

for i = 1:length(SNR)
    noise(:,i) = awgn(sig,SNR(i)+15);
    noise(:,i) = noise(:,i) - sig;
    osh_in(i) = 10*log10(mean(sig.^2)/mean(noise(:,i).^2));
    signal_shift = shift_plane(sig,phi_sig,teta_sig,p_loc,fs);
    noise_shift = shift_plane(noise(:,i),phi_noise,teta_noise,p_loc,fs);

    [y_signal,W_n] = spat_filt_wb_time_lc_lms(signal_shift, L, N, mu);
    [y_noise,W_n] = spat_filt_wb_time_lc_lms(noise_shift, L, N, mu);
    
    y_signal = y_signal(m:end);
    y_noise = y_noise(m:end);
    plot(noise)
    plot(y_noise)
    osh_out(i) = 10*log10(mean(y_signal.^2)/mean(y_noise.^2));
    osh_out(i) = osh_out(i) - osh_in(i);
end

%%
figure()
plot(osh_in, osh_out)
