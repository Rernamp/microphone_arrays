function [PCM] = pdm2pcm_mcu(PDM)

%% CIC decimator x32
R = 32;         %коэффициент децимации
out_CIC = func_cic_n4(PDM,R,32);  
out_CIC = floor(out_CIC/(16));

%% FIR compensator - decimator x2
load('ImResp_CompFIR.mat', 'bhi');
out_CIC_FIR = filter(bhi,1,out_CIC);
out_CIC_FIR = out_CIC_FIR(1:2:end);
out_CIC_FIR = floor(out_CIC_FIR/(2^15));


%% Hight-Pass IIR biquad filter 
[b,a] = butter(2,150/8000,'high'); 
SOS = [b a];
SOS = round(32767*SOS)/32767;

biquad = dsp.BiquadFilter('Structure','Direct form I',...
    'SOSMatrix',SOS,'ScaleValues',1);
PCM = biquad(out_CIC_FIR.');
end

