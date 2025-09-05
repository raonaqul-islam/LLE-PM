%% POWER SPECTRAL DENSITY OF NO-MODULATION CASE
% This subroutine loads the centroid data from a .mat file and calculates
% the power spectral density for the case where there is no-modulation to
% the pump. It also calculates it analytically and plots them on top of
% each other.

clear
close all

% ---- Load centroid data ---- %
load ./data/mc_nopm_090125.mat                 

% ---- Parameters ---- %
R      = 2.5e-3;                                   % Radius of the microresonator
N      = 5e6;                                      % No. of steps
tph    = 2.5e-6;                                   % Photon lifetime
dt     = 1e-4;                                     % Time step in LLE
dt     = dt*tph;                                   % De-normalizing with photo-lifetime
theta0 = nopmcircmean;                               % Average
fs     = 1/dt;                                     % Sampling frequency
theta1 = unwrap(theta0,[],1);                      % Unwrapping the centroids   
psd1   = zeros(N-1,size(theta0,2));                % Initiate array for repetition-rate psd

% ---- Calculate velocity ---- %
v     = diff(theta1,[],1)/(dt);                  % Calculate velocity
v_fft = fft(v,[],1);                               % Take fourier transform of velocity

% ----- Calculate PSD ---- %
Nv    = size(v_fft,1);                             % Array size

for kk = 1:size(v_fft,2)
    psd1(:,kk)  = abs(v_fft(:,kk)).^2/(Nv*fs);     % PSD without windowing    
end
psd1  = sum(psd1,2)/sqrt(size(psd1,2));            % Take average

% ----- Take only positive frequencies ---- %
if rem(Nv,2)==0
    % If array size is even
    psd1_onesided = psd1(1:Nv/2+1);
    psd1_onesided(2:end-1) = 2*psd1_onesided(2:end-1);
    f_one_sided = (0:Nv/2)'*(fs/Nv);
else
    % If array size is odd
    psd1_onesided = psd1(1:(Nv+1)/2);
    psd1_onesided(2:end) = 2*psd1_onesided(2:end);
    f_one_sided = (0:(Nv-1)/2)'*(fs/Nv);
end

% ---- Plot LLE-PSD ---- %
fig = figure('Units','inches','Position',[5 5 12 5],'Color','White');
set(gcf,'DefaultLineLineWidth',1);
loglog(f_one_sided,psd1_onesided,'r-');
hold on;

% ---- Calculate and plot PSD analytically ---- %
nu       = f_one_sided;    % Frequency in Hz
kappa    = 1/tph;
gamma_n  = 2e5;     % Dissipation factor
eta      = 1.65e6;   % Variation in repetition rate
A        = 10;             % Noise amplitude
PSDnopm  = 2*eta^2*A./(gamma_n^2+(2*pi*nu).^2);
loglog(nu,PSDnopm,'k--');

% ---- Plot settings ---- %
set(gca,'FontSize',14,'FontName','Times New Roman','LineWidth',1.2,'TickDir','Out');
yticklabels(strrep(yticklabels,'-',char(0x2212)))
xticklabels(strrep(xticklabels,'-',char(0x2212)))
title('No-modulation (Mone-Carlo)')
box off;
ylabel('Power Spectral Density');
xlabel('Frequency (Hz)');
xline(gamma_n/(2*pi),'Color','k','alpha',0.2,'LineStyle','--','LineWidth',3);
legend('LLE','Analytical')

% ---- Save data ---- %
psdnopm.LLE = psd1_onesided;
psdnopm.analytic = PSDnopm;
psdnopm.nu = f_one_sided;
save('./data/psdnopm2.mat','psdnopm');