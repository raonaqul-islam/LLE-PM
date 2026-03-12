%% POWER SPECTRAL DENSITY OF FEEDBACK-MODULATION CASE

clc
clear
close all

% *************************************************************************
%                        SIMULATION PARAMETERS                            %
% *************************************************************************
issave       = true;                              % If you want to save data
N_sim        = 120;                                % No. of simulations
N            = 1e6;                                % No. of samples
R            = 225e-6;                             % Radius of the microresonator
dt           = 1e-2;                               % Discretization
tph          = 8.23e-10;                           % Photon lifetime
kappa        = 1/tph;                              % Cavity decay rate
dt           = tph*dt;                             % De-normalized time-step
fs           = 1/(100*dt);                               % Sampling frequency
skip         = 0;                                  % No. of samples to skip

% *************************************************************************
%                     LOAD DATA AND CALCULATE PSD                         %
% *************************************************************************
Nw = length(skip+1:N)-1;
Swrep_all = zeros(Nw,N_sim);
for ik = 1:N_sim
    load(sprintf('../data/mc_extpm_112325/mc_extpm_112325_%d.mat',ik));    
    circmean  = unwrap(circmean);                  % Unwrap to avoid boundary jumps
    dwrep     = diff(circmean(skip+1:N))/(dt)*R;   % Calculate f_rep noise
    dwrep_fft = fft(dwrep);                        % Fourier transform of f_rep noise   
    Swrep     = abs(dwrep_fft).^2/(Nw*fs);         % Power spectral density of f_rep noise
    Swrep_all(:,ik) = Swrep;                       % Store current iteration    
    clear circmean theta0 dwrep dwrep_fft Swrep    % Clear variables to save memory 
end

% Averaged frequency noise
Swrep  = sum(Swrep_all,2)/sqrt(N_sim);

% Take only positive frequencies
if rem(Nw,2)==0
    % If array size is even
    Swrep_onesided = Swrep(1:Nw/2+1);
    Swrep_onesided(2:end-1) = 2*Swrep_onesided(2:end-1);
    f_one_sided = (0:Nw/2)'*(fs/(Nw));    
else
    % If array size is odd
    Swrep_onesided = Swrep(1:(Nw+1)/2);
    Swrep_onesided(2:end) = 2*Swrep_onesided(2:end);
    f_one_sided = (0:(Nw-1)/2)'*(fs/(Nw));    
end

% *************************************************************************
%                            PLOT NUMERICAL                               %
% *************************************************************************
fig = figure('Units','inches','Position',[5 5 12 5],'Color','White');
set(gcf,'DefaultLineLineWidth',2);

loglog(f_one_sided,Swrep_onesided,'r-');
xline(kappa*0.13,'Color','Black','LineStyle','--','Linewidth',1);
ylabel('Frequency Noise (Hz^2/Hz)');
xlabel('Offset Frequency (Hz)');
title('External PM')
% xlim([1e3 1e11]);
hold on

set(gca,'FontSize',14,'FontName','Arial','LineWidth',1,'TickDir','Out');
yticklabels(strrep(yticklabels,'-',char(0x2212)))
xticklabels(strrep(xticklabels,'-',char(0x2212)))
box off

% Save data
if issave == true
    num_extpm.LLE      = Swrep_onesided;    
    num_extpm.nu       = f_one_sided;    
    save('../data/psdextpm_112925.mat','num_extpm');
end
