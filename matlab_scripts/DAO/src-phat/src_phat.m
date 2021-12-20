clear;close all;

c = 343;
f_re = 4000;
lamda = c/f_re;
K = 2;
d = 0.04;


thetaan = 0; 

% [noise,fs] = audioread('sound_44.1kHz.wav');
[noise,fs] = audioread('sound.wav');

t = 1;

noise = noise(1:t*fs);



tau = (d*sind(thetaan))/(c)
% tau = 0.8/fs;

noise_shift = zeros(K, length(noise));

%%
for k = 1:K
    [y_shift] = shift( noise, tau*(k-1), fs);

    noise_shift(k,:) = y_shift;
end


t = 0:1/fs:2-1/fs;

a = [0.5 0.5];

X_fft = [fft(noise_shift(1,:)); fft(noise_shift(2,:))];

%%
R = [];
P = [];

P_k = [];
for i = 1:400
    for phi = -30:30
        

        R_i = X_fft(:,i)*X_fft(:,i)';

        v_k = exp(-1i*(2*pi*d*[0:K-1]*sind(phi)*i*fs)/(length(noise_shift(1,:))*c))';
        P_k_i = v_k'*R_i*v_k;
        P_k = [P_k ;P_k_i];

    end

    P = [P P_k];
    P_k = [];
    
end


%%
P_ab = db(abs(P));
figure()
hold on;

plot(db(abs(fft(noise))));

%%
[w,phi] = meshgrid([1:400]*fs/length(noise_shift(1,:)),-30:30);

figure()
surf(w,phi,P_ab);
shading interp 
grid on %
colormap gray
