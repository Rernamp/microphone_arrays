clear;close all;


K = 9; 
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
noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);

sig_in_MR = signal_shift + noise_shift;
mu = 1;

time = 0:1/fs:(length(noise)-1)/fs;

mu = 1;

time = 0:1/fs:(length(noise)-1)/fs;
%%
phi = -90:2:90;


for l = 1:length(phi)
    noise_shift = shift_plane(noise,phi(l),teta_noise,p_loc,fs);
    sig_in_MR = signal_shift + noise_shift;
    [y_cl,y,W_n] = spat_filt_wb_time_lc_rls_pesq(signal_shift,sig_in_MR,L,K);
    k_phi(l,1:2) = pesqbin(y_cl,y,fs,'nb');
    
   
end
%%
teta = 0:2:90;
for l = 1:length(teta)
    noise_shift = shift_plane(noise,phi_noise,teta(l),p_loc,fs);
    sig_in_MR = signal_shift + noise_shift;
    [y_cl,y,W_n] = spat_filt_wb_time_lc_rls_pesq(signal_shift,sig_in_MR,L,K);
    k_teta(l,1:2) = pesqbin(y_cl,y,fs,'nb');
    
   
end
%%


figure()

plot(phi,k_phi(:,1))
grid on
title("PESQ для LC RLS")
xlabel("J , K = 4")
ylabel("PESQ")



figure()

plot(teta,k_teta(:,1))
grid on
title("PESQ для LC RLS")
xlabel("K, J = 32")
ylabel("PESQ")

save('PESQ_RLS.mat','phi','teta','k_phi','k_teta');
