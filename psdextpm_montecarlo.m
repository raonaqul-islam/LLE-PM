%% POWER SPECTRAL DENSITY OF EXTERNAL-MODULATION CASE
% This subroutine loads the centroid data from a .mat file and calculates
% the power spectral density for the case where there is no-modulation to
% the pump. It also calculates it analytically and plots them on top of
% each other.

clear
close all

% ---- Load centroid data ---- %
load ./data/mc_extpm_020925.mat                  

% ---- Parameters ---- %
R      = 2.5e-3;
tph    = 2.5e-6;                                   % Photon lifetime
N      = 5e6;                                      % No. of steps
dt     = 1e-4;                                     % Time step
dt     = dt*tph;                                   % Time-step de-normalized
theta0 = extpmcircmean;                              % Average
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
set(gcf,'DefaultLineLineWidth',3);
loglog(f_one_sided,psd1_onesided,'r-');
hold on;

% ---- Calculate and plot PSD analytically ---- %
kappa     = 1/tph;
nu        = f_one_sided;
eta       = 1.65e6;
gamma_n   = 2e5;                           % Normalized gamma_n 
A         = 10;                            % Normalized noise amplitude
lda       = 0.2;                           % Position shifting eigenvalue
PSDextpm  = (2*eta^2*A.*(2*pi*nu).^2)./((gamma_n^2+(2*pi*nu).^2).*(kappa^2*lda^2+(2*pi*nu).^2));
loglog(nu,PSDextpm,'k--');

% ---- Plot settings ---- %
legend('LLE','Analytical')
set(gca,'FontSize',14,'FontName','Times New Roman','LineWidth',1.2,'TickDir','Out');
yticklabels(strrep(yticklabels,'-',char(0x2212)))
xticklabels(strrep(xticklabels,'-',char(0x2212)))
title('External-modulation (Mone-Carlo)')
box off;
ylabel('Power Spectral Density');
xlabel('Frequency (Hz)');
xline(gamma_n/(2*pi),'Color','k','alpha',0.2,'LineStyle','--','LineWidth',3);

% ---- Save data ---- %
psdextpm.LLE = psd1_onesided;
psdextpm.analytic = PSDextpm;
psdextpm.nu = f_one_sided;
save('./data/psdextpm2.mat','psdextpm');