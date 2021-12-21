function [PCM] = pdm2pcm_matlab(PDM)
%% CIC decimator x32
R = 32;         %коэффициент децимации
D = 2;          %задержка
N = 4;          %порядок фильтра

cicdec = dsp.CICDecimator(R,D,N);

out_CIC = cicdec(PDM);  
out_CIC = floor(out_CIC/(16));

%% FIR compensator - decimator x2
load('ImResp_CompFIR.mat', 'bhi');
out_CIC_FIR = filter(bhi,1,out_CIC);
out_CIC_FIR = out_CIC_FIR(1:2:end);
out_CIC_FIR = floor(out_CIC_FIR/(2^15));


%% Hight-Pass FIR filter 
HP = fir1(256,300/8000,'high');
PCM = filter(HP,1,out_CIC_FIR);
end

