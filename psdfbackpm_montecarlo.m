%% POWER SPECTRAL DENSITY OF FEEDBACK-MODULATION CASE
% This subroutine loads the centroid data from a .mat file and calculates
% the power spectral density for the case where there is no-modulation to
% the pump. It also calculates it analytically and plots them on top of
% each other.

clear
close all

% ---- Load centroid data ---- %
load ./data/mc_fbackpm_090425.mat              

% ---- Parameters ---- %
R      = 2.5e-3;
tph    = 2.5e-6;                                   % Photon lifetime
delayedsteps = 1e6;
N      = 6e6-delayedsteps;                         % No. of steps
dt     = 1e-4;                                     % Time step in normalized form   
dt     = tph*dt;                                   % De-normalized time-step
theta0 = fbackpmcircmean(delayedsteps+1:end,:);    % Skipped delayed steps
fs     = 1/dt;                                     % Sampling frequency
theta1 = unwrap(theta0,[],1);                      % Unwrapping the centroids   
psd1   = zeros(N-1,size(theta0,2));                % Initiate array for repetition-rate psd

% ---- Calculate velocity ---- %
v     = diff(theta1,[],1)/(dt);                      % Calculate velocity
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
% set(gcf,'DefaultLineLineWidth',3);
loglog(f_one_sided,psd1_onesided,'r-','LineWidth',3);
hold on;

% ---- Calculate and plot PSD analytically ---- %
kappa     = 1/tph;
nu        = f_one_sided;
omega_d   = 2*pi*nu;
eta       = 1.65e6;
gamma_n   = 2e5;                           % Normalized gamma_n 
A         = 10;                            % Normalized noise amplitude
delay     = 1e4;                           % No. of time-steps used as feedback delay
tf        = delay*dt;                      % Delay de-normalized
lda       = 0.2;                           % Position shifting eigenvalue
% lda       = lda/delay;                        % Average stiffness reduced due to delay
mod_term   = (exp(-1i.*omega_d.*tf)-1);    % Numerator term 
PSDfbackpm = (2*eta^2.*(omega_d.^2)*A).*abs(mod_term).^2./((gamma_n^2+(omega_d).^2).*(kappa^2*lda^2+(omega_d).^2));
loglog(f_one_sided,PSDfbackpm,'k-','LineWidth',1);
hold on;
xline(gamma_n/(2*pi),'Color','k','alpha',0.2,'LineStyle','--','LineWidth',3);

% ---- Plot settings ---- %
legend('LLE','Analytical')
set(gca,'FontSize',14,'FontName','Times New Roman','LineWidth',1.2,'TickDir','Out');
yticklabels(strrep(yticklabels,'-',char(0x2212)))
xticklabels(strrep(xticklabels,'-',char(0x2212)))
title('Feedback-modulation (Mone-Carlo)')
box off
ylabel('Power Spectral Density');
xlabel('Frequency (Hz)');

% ---- Save data ---- %
psdfbackpm.LLE = psd1_onesided;
psdfbackpm.analytic = PSDfbackpm;
psdfbackpm.nu = f_one_sided;
save('./data/psdfbackpm2.mat','psdfbackpm');
