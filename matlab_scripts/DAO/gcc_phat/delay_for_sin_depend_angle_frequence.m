clear;close all;

c = 343;
f = 4e3;
lamda = c/f;
K = 2;
d = 0.1;

fs = [8 16 32 44.1].*1e3;
tau_ff = [];
thetaan = -76:76;
for ff = 1:length(fs)
    time_sin = 0.5;
    M = time_sin*fs(ff);
    N = M*2 - 1;

    tau = (d*sind(thetaan))/(c);
    t = 0:1/fs(ff):time_sin - 1/fs(ff);
    y = sin(2*pi*f*t);
    y = awgn(y,20);
    Nfft = 2^nextpow2(N);

    tau_f = [];

    for i = 1:length(tau)
        y_mas_shift = [];

        for k = 1:K
           [y_shift] = shift( y, tau(i)*(k-1), fs(ff));
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

    %     [r,tau_two] = gccphat_two( y_mas_shift(2,:)', y_mas_shift(1,:)', fs(ff) );

        %     tau_f = [tau_f tau_est(2)];
        tau_f = [ tau_f ((N+1)/2 - max_in)/fs(ff)];

    end
    tau_ff = [tau_ff ; tau_f];
end
%%

figure()
hold on;
plot(thetaan,tau)
for ff = 1:length(fs)
    
    
    
    plot(thetaan,tau_ff(ff,:))
    
    
end
legend("Orig tau", ["f = " + num2str(fs(1))], ["f = " + num2str(fs(2))], ["f = " + num2str(fs(3))], ["f = " + num2str(fs(4))])
xlabel("\phi , град")
ylabel("time , c")
grid on

figure()
hold on;
plot(thetaan,thetaan)
for ff = 1:length(fs)
    thetaan_exp = asind((tau_ff(ff,:)*c)/(d));
    plot(thetaan,thetaan_exp)
    
end

grid on
legend("Orig \phi", ["f = " + num2str(fs(1))], ["f = " + num2str(fs(2))], ["f = " + num2str(fs(3))], ["f = " + num2str(fs(4))])
xlabel("\phi , град")
ylabel("\phi , град")

%%




for ff = 1:length(fs)
    
    figure()
    hold on;
    plot(thetaan,tau)
    plot(thetaan,tau_ff(ff,:))
    legend("Orig tau", ["f = " + num2str(fs(ff))])
    xlabel("\phi , град")
    ylabel("time , c")
    grid on
    
end



%%

for ff = 1:length(fs)
    thetaan_exp = asind((tau_ff(ff,:)*c)/(d));
    figure()
    hold on;
    plot(thetaan,thetaan)
    plot(thetaan,thetaan_exp)
    grid on
    legend("Orig \phi", ["f = " + num2str(fs(ff))])
    xlabel("\phi , град")
    ylabel("\phi , град")

end

