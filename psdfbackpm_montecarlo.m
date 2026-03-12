% POWER SPECTRAL DENSITY OF FEEDBACK-MODULATION CASE

tic
clc
clear
close all

% *************************************************************************
%                        SIMULATION PARAMETERS                            %
% *************************************************************************
issave       = true;                               % If you want to save data
N_sim        = 60;                                 % No. of simulations
N            = 5e6;                                % No. of samples
R            = 225e-6;                             % Radius of the microresonator
dt           = 1e-2;                               % Discretization
tph          = 8.23e-10;                           % Photon lifetime
dt           = tph*dt;                             % De-normalized time-step
fs           = 1/dt;                               % Sampling frequency
skip         = 1e6;                                % No. of samples to skip

% *************************************************************************
%                     LOAD DATA AND CALCULATE PSD                         %
% *************************************************************************
Nw = N-skip-1;
Swrep_all = zeros(Nw,N_sim);

for ik = 1:N_sim
    data = load(sprintf('../data/mc_fbackpm_121025/mc_fbackpm_121025_%d.mat',ik),'circmean');    
    circmean  = unwrap(data.circmean(1:end));           % Unwrap to avoid boundary jumps
    dwrep     = diff(circmean(skip+1:N))/(dt)*R;   % Calculate f_rep noise
    dwrep_fft = fft(dwrep);                        % Fourier transform of f_rep noise   
    Swrep     = abs(dwrep_fft).^2/(Nw*fs);         % Power spectral density of f_rep noise
    Swrep_all(:,ik) = Swrep;                       % Store current iteration    
    % clear circmean theta0 dwrep dwrep_fft Swrep    % Clear variables to save memory 
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
ylabel('Frequency Noise (Hz^2/Hz)');
xlabel('Offset Frequency (Hz)');
title('Feedback-modulation')
xlim([1e3 1e11]);
hold on

set(gca,'FontSize',14,'FontName','Arial','LineWidth',1.2,'TickDir','Out');
yticklabels(strrep(yticklabels,'-',char(0x2212)))
xticklabels(strrep(xticklabels,'-',char(0x2212)))
box off

% Save data
if issave == true
    num_fbackpm.LLE      = Swrep_onesided;    
    num_fbackpm.nu       = f_one_sided;    
    save('../data/psdfbackpm_121025.mat','num_fbackpm');
end

elapsed_time = toc;
