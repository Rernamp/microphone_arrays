function [signal_shift] = shift_plane(signal,phi,theta,p_el,fs)
c = 343;

signal_fft = fft(signal);
signal_fft(1) = 0;
signal = ifft(signal_fft);
a = [-sind(theta)*cosd(phi) -sind(theta)*sind(phi) -cosd(theta)];

tau = a*p_el/c;

signal_shift = zeros(length(p_el(1,:)),length(signal));

for i = 1:length(p_el(1,:))
    signal_shift(i,:) = shift(signal,tau(i),fs);
end

end

