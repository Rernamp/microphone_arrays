function [ y,W ] = func_LMS(s_in, J, M, mu )

x = zeros(M*J,1);
h = zeros(M*J,1);
y = zeros(1,length(s_in));
for i = 1:length(s_in) 

    x = [s_in(:,i); x(1:end-M)];
    y(i) = h' * x;
    
    alpha = - h' * x;

    h = h + 2*mu*x*alpha;
    
    
end
W = h;


end

