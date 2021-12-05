function [signal_shift] = shift_plane_tau(signal,fs,tau,p_el)
c = 343;

signal_fft = fft(signal);
signal_fft(1) = 0;
signal = ifft(signal_fft);


signal_shift = zeros(length(p_el(1,:)),length(signal));

for i = 1:length(p_el(1,:))
    signal_shift(i,:) = shift(signal,tau(i),fs);
end

end

