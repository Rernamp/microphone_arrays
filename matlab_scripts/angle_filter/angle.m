clear;
close all;

N = 9; 
c = 343;

fs = 8000;

f = 2000;
phi_sig = 10;
teta_sig = 30;
phi_noise = 60;
teta_noise = 0;

c = 343;
d = c/(2*f);
d = 0.04;
L = 32;
N_f = 50; 
mu = 1;
c_noise =35;
p_loc = gen_place_el(sqrt(N),sqrt(N),d,d,1)';

[noise,fs] = audioread('nois.wav');
[signal,fs] = audioread('speech_dft_8kHz.wav');

t = 2;
time = 0:1/fs:t-1/fs;
signal = signal(1:t*fs);
noise = noise(1:t*fs);

signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);

sig_in_MR = signal_shift+noise_shift;

a_phi = [-cosd(teta_sig).*cosd(phi_sig) ; -cosd(teta_sig).*sind(phi_sig) ; -sind(teta_sig)];
a_p = a_phi'*p_loc;

f_up = 4000;    
f_down = 0; 
f_den = 0.5;
f_pos = f_down:f_den:f_up;

w_direct = exp(-1i*2*pi*f_pos'.*a_p/c);

w_2zone = conj(w_direct(end:-1:2,:));
% w_direct(end-1,:) = real(w_direct(end,:));

w_w = [ w_direct ; w_2zone];

imp_har1 = ifft(w_w);
imp_har1 = fftshift(imp_har1);
imp_dec = fftshift(imp_har1(8000-(L/2 - 1):8000+(L/2),:));

J = L;

C = zeros(N*J,J);

for j = 1:J
   C(:,j) = [zeros(1,(j-1)*N) ones(1,N) zeros(1,J*N-j*N)].';
end

w = reshape(imp_dec',1,[])';

f = C'*w;
mu = 0.1;
[ y ] = GSC_NLMS(sig_in_MR, J, N, mu,f );
%%
figure()
hold on
plot(sig_in_MR(1,:))
plot(y)

sound(y,fs)

% figure()
% hold on
% for i = 1:length(imp_har1(1,:))
%     plot((imp_dec(:,i)))
% end
