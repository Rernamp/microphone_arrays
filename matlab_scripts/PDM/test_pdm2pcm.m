clc; clear; close all;
% load('pcm_exp.mat','PCM')
load('mic1.txt')
load('mic2.txt')

% pcm_exp = PCM;
% pdm_exp = PDM;
pdm_exp = mic1' + 128;
pcm_exp = mic2' + 128;


% _exp      - experimental data
% _matlab   - matlab functions result
% _mcu      - microcontroller simulation result

%% cut off the first 50 ms
t_cut = 40;

pcm_err_ms = 8*t_cut; 
pdm_err_ms = 64*t_cut;

pcm_exp(1:pcm_err_ms) = [];
pdm_exp(1:pdm_err_ms) = [];

pdm_exp_bi = fliplr( de2bi(pdm_exp) );

pdm_exp_bi_rshp = reshape(pdm_exp_bi.',[],1);

for ii =1:length(pdm_exp_bi_rshp)
    if(pdm_exp_bi_rshp(ii) == 0)
        pdm_exp_bi_rshp(ii) = -1;
    end
end

pcm_matlab = pdm2pcm_matlab(pdm_exp_bi_rshp);
pcm_mcu = pdm2pcm_mcu(pdm_exp_bi_rshp);

pcm_exp = pcm_exp/max(pcm_exp);
pcm_matlab = pcm_matlab/max(pcm_matlab);
pcm_mcu = pcm_mcu/max(pcm_mcu);

figure(555)
hold on
plot(pcm_matlab,'DisplayName','matlab')
% plot(pcm_exp,'DisplayName','exp')
plot(pcm_mcu,'DisplayName','mcu')
legend()


