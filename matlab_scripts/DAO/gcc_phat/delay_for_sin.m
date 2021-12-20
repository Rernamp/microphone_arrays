clear;close all;

c = 343;
f = 1e2;
lamda = c/f;
K = 2;
d = 0.1;

fs = 8e3;
time_sin = 0.5;
N = time_sin*fs*2;
thetaan = -90;
tau = (d*sind(thetaan))/(c);

t = 0:1/fs:time_sin - 1/fs;
y = sin(2*pi*f*t);


tau_f = [];
y_mas_shift = [];

for k = 1:K
    [y_shift] = shift( y, tau*(k-1), fs);
    y_mas_shift = [ y_mas_shift ;y_shift] ;
end

fft_y_mas_shift = [fft(y_mas_shift(1,:), N) ; fft(y_mas_shift(2,:), N)];


R = (fft_y_mas_shift(1,:).*conj(fft_y_mas_shift(2,:)));
R = R./abs(R);
R_t = fftshift(ifft(R));
%     R_t = R_t(end:-1:1);
[max_val max_in] = max(R_t);
tau_est = gccphat(y_mas_shift',y');
%     tau_f = [tau_f tau_est(2)];
tau_f = [tau_f (N/2 - max_in)/fs];

%%
tau_f = tau_f;
error_tau = abs(tau - tau_f);
