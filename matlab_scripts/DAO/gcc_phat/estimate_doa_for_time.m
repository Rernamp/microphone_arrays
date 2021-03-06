clear;close all;

c = 343;
f = 4e3;
lamda = c/f;
K = 2;
d = 0.10;

fs = 16e3;
time_sin = 0.5;
M = time_sin*fs;



t = 0:1/fs:time_sin - 1/fs;
y = sin(2*pi*f*t);
y = awgn(y,20);
N_data_one = 500;
N = N_data_one*2 - 1;
Nfft = 2^nextpow2(N);
N_iter = M/N_data_one;
tau_f = [];

phi_n = 7;
thetaan = 90:-phi_n:90-15*phi_n;
tau = (d*sind(thetaan))/(c);

for i = 1:N_iter
    
    y_mas_shift = [];

    for k = 1:K
       [y_shift] = shift( y, tau(i)*(k-1), fs);
       y_mas_shift = [ y_mas_shift ;y_shift] ;
    end
    
    
    fft_y_mas_shift = [fft(y_mas_shift(1,(i-1)*N_data_one + 1:i*N_data_one),N) ; fft(y_mas_shift(2,(i-1)*N_data_one + 1:i*N_data_one),N)];
    R = (fft_y_mas_shift(1,:).*conj(fft_y_mas_shift(2,:)));
    R = R./(abs(R));
    R_t = fftshift(ifft(R));
    [max_val max_in] = max((R_t));
    tau_find = ((N+1)/2 - max_in)/fs;
    tau_f = [tau_f tau_find];

end

%%

figure()
hold on
plot(thetaan,tau)
plot(thetaan,tau_f)
