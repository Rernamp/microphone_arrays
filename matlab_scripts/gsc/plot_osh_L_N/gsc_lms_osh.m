clear;
close all;

N = 9; 
c = 343;

fs = 8000;

f = 2000;
phi_sig = 0;
teta_sig = 0;
phi_noise = 60;
teta_noise = 0;
c = 343;
J = 32;
% d = c/(2*f);
t = 2;
time = 0:1/fs:t-1/fs;
d = 0.04;
N_f = 50; 
mu = 0.01;
p_loc = gen_place_el(sqrt(N),sqrt(N),d,d,1)';

[noise,fs] = audioread('nois.wav');
[signal,fs] = audioread('speech_dft_8kHz.wav');

t = 2;
m = 1;



signal = signal(1:t*fs);
noise = noise(1:t*fs);
osh_in = mean(signal.^2)/mean(noise.^2);

signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);
        
mu_GSC =  0.1;
mu_LMS = 1;
L_max = 64;
L_min = 16;
L_del = 2;
L = L_min:L_del:L_max; 
i = 1;
%%
for l = L_min:L_del:L_max
    [y_noise_GSC] = GSC_NLMS(noise_shift, l,N, mu_GSC);
    [y_sig_GSC] = GSC_NLMS(signal_shift, l, N, mu_GSC);
    
    [y_noise_LMS] = func_LC_NLMS(noise_shift, l,N, mu_LMS);
    [y_sig_LMS] = func_LC_NLMS(signal_shift, l, N, mu_LMS);
    
    osh_l_GSC(i) = mean(y_sig_GSC(m:end).^2)/mean(y_noise_GSC(m:end).^2);
    osh_l_LMS(i) = mean(y_sig_LMS(m:end).^2)/mean(y_noise_LMS(m:end).^2);
    i = i+1;
   
   
end
%%
k = [4 9 16 25];
L_k = 20;

for p = 1:length(k)
    
    
    p_loc = gen_place_el(sqrt(k(p)),sqrt(k(p)),d,d,1)'; 
    signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
    noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);
        
    
    [y_noise_GSC] = GSC_NLMS(noise_shift, L_k,k(p), mu_GSC);
    [y_sig_GSC] = GSC_NLMS(signal_shift, L_k, k(p), mu_GSC);
    
    [y_noise_LMS] = func_LC_NLMS(noise_shift, L_k,k(p), mu_LMS);
    [y_sig_LMS] = func_LC_NLMS(signal_shift, L_k, k(p), mu_LMS);

    osh_k_GSC(p) = mean(y_sig_GSC(m:end).^2)/mean(y_noise_GSC(m:end).^2);
    osh_k_LMS(p) = mean(y_sig_LMS(m:end).^2)/mean(y_noise_LMS(m:end).^2);
    
end
%%

SNR_L_gsc = db(osh_l_GSC-osh_in);
SNR_L_lms = db(osh_l_LMS-osh_in);

figure()
hold on
plot(L,SNR_L_gsc)
plot(L,SNR_L_lms)
legend('GSC','LMS')
xlabel('L')
ylabel('SNR, dB')
grid on
%%

SNR_K_gsc = db(osh_k_GSC-osh_in);
SNR_K_lms = db(osh_k_LMS-osh_in);

figure()
hold on
plot(k,SNR_K_gsc)
plot(k,SNR_K_lms)
legend('GSC','LMS')
xlabel('N')
ylabel('SNR, dB')
grid on
