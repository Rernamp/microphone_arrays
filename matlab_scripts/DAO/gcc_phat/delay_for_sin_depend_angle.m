clear;close all;

c = 343;
f = 4e3;
lamda = c/f;
K = 2;
d = 0.10;

fs = 16e3;
time_sin = 0.5;
M = time_sin*fs;
N = M*2 - 1;
thetaan = -90:90;
tau = (d*sind(thetaan))/(c);
t = 0:1/fs:time_sin - 1/fs;
y = sin(2*pi*f*t);
y = awgn(y,20);
Nfft = 2^nextpow2(N);

tau_f = [];

for i = 1:length(tau)
    y_mas_shift = [];

    for k = 1:K
       [y_shift] = shift( y, tau(i)*(k-1), fs);
       y_mas_shift = [ y_mas_shift ;y_shift] ;
    end

    fft_y_mas_shift = [fft(y_mas_shift(1,:),N) ; fft(y_mas_shift(2,:),N)];
    % fft_y_mas_shift = [fft(y_mas_shift(1,:)) ; fft(y_mas_shift(2,:))];

%     a = 0.1;
    R = (fft_y_mas_shift(1,:).*conj(fft_y_mas_shift(2,:)));
    R = R./(abs(R));
    % R_t = (ifft(R));
    R_t = fftshift(ifft(R));
    % R_t = R_t(Nfft/2+1-(M-1)/2:Nfft/2+1+(M-1)/2);
    % R_t = R_t(end:-1:1);
    [max_val max_in] = max((R_t));
%   tau_est = gccphat(y_mas_shift',y');

%     [r,tau_two] = gccphat_two( y_mas_shift(2,:)', y_mas_shift(1,:)', fs );

    %     tau_f = [tau_f tau_est(2)];
    tau_f = [ tau_f ((N+1)/2 - max_in)/fs];
    
end

%%

error_tau = abs(tau - tau_f);

thetaan_exp = asind((tau_f*c)/(d));
figure()
hold on;
plot(thetaan,tau)
plot(thetaan,tau_f)
grid on
legend("Orig tau", "TDOA")
xlabel("\phi , град")
ylabel("time , c")

figure()
hold on;
plot(thetaan,thetaan)
plot(thetaan,thetaan_exp)
grid on
legend("Orig %phi", "DOA")
xlabel("\phi , град")
ylabel("\phi , град")
