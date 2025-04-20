% LLE time evolution including phase modulation
% This version uses generic microresonator parameters, again, and shifts
% the initial soliton to -pi/2. It is observed that after 100 photon
% lifetimes the soliton moves from -pi/2 to the stable position pi/2.
% -----------------------------+
% Raonaqul Islam, UMBC         |
% Date started: March 19, 2025 |
% Last updated: April 18, 2025 |
% -----------------------------+

clear
clc
close all

% ************
% Parameters
% ************
% Parameters captured from Pradyoth

alpha   = 3.5;                            % Detuning
beta    = -0.1549;                        % Dispersion, beta2 only in this case
F       = 2.3;                            % Pump amplitude
gamma   = 1;                              % Normalized nonlinear coefficient
deltaM  = 0.65;                           % Phase modulation

% ******************************
% Spatial domain discretization
% ******************************
N       = 512;                            % No. of points on axis
naxis   = (-N/2:N/2-1).';                 % General axis
dtheta  = 2*pi/N;                         % Spatial step-size
theta   = naxis*dtheta;                   % Spatial (Azimuthal) domain
dmu     = 1;                              % Mode number domain step-size
mu      = fftshift(dmu*naxis);            % Mode number domain

% ****************************
% Power with input modulation
% ****************************
F_tilde = fft(F*exp(1i*deltaM*sin(theta))); % Fourier transform of power 
                                            % including the modulation term

% ***************************                                            
% Time domain discretization
% ***************************
dt      = 1e-4;                             % Temporal step size
N_tau_p = 200;                              % No. of photon lifetimes
N_steps = N_tau_p/dt;                       % No. of steps

% Generate figure to plot on
figure('Units','inches','Position',[2 2 8 4]); 

% ****************
% Time evolution 
% ****************

% Shifting previous solution to -pi/2
load ./data/sol0.mat                        % Load previous solution
psi_0 = circshift(psi_0,-N/2);              % Using previous solution shifted to -pi/2

for ix = 1:N_steps  
    
    % Nonlinear part (half-step)        
    psi_nl      = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));

    % Linear part (full-step)
    psi_0_tilde = fft(psi_nl);    
    A_tilde     = -(1+1i*alpha)+1i.*beta.*(mu).^2;      
    psi_l_tilde = (psi_0_tilde+F_tilde./A_tilde).*exp(A_tilde*dt)-F_tilde./A_tilde;
    psi_l       = ifft(psi_l_tilde);
    
    % Nonlinear part (half-step)
    psi_0       = psi_l;
    psi_nl2     = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));
    psi_0       = psi_nl2; 

    % Output after one-step
    psi_out     = psi_nl2;
    
    % Plot solution for every photon lifetime (1/dt)
    if rem(ix,1/dt) == 0
        % transparency = ix/(2*N_steps);
        % plot(theta/pi,abs(psi_out).^2,'LineWidth',2.5,'Color',[0.07 0.62 1 transparency])
        plot(theta/pi,abs(psi_out).^2,'LineWidth',2.5,'Color',[0.07 0.62 1]);
        xlabel('\theta/\pi')
        ylabel('|\psi|^2')
        title(['No. of \tau_{p} = ' num2str(ix*dt)])
        xlim([-1 1]);
        set(gca,'FontSize',16,'FontName','Times New Roman','LineWidth',1.2)  
        xticklabels(strrep(xticklabels,'-','−'))
        xline(0.5,'LineWidth',2,'Color','k','LineStyle','--')
        grid on
        box off
        % hold on
        drawnow
    end

end

% save('./data/sol0.mat','psi_out');

% ************************************
% Plot soliton in mode number domain
% ************************************
% Prepare the solution
psi_fft  = abs(fftshift(fft(psi_out))).^2;  % Obtain fourier transform of the solution
psi_norm = psi_fft/max(psi_fft);            % Normalize the output
psi_db   = 10*log10(psi_norm);              % Calculate normalized output power in dB

% Plot in stem
figure('Units','inches','Position',[2 2 15 4]);
% plot(fftshift(mu),psi_db,'LineWidth',3,'Color','Blue');
soliton = stem(fftshift(mu),psi_db,'LineWidth',2,'Marker','None','Color','Red');
soliton.BaseValue = -400; % Turn the plot upside down
xlabel('Relative mode number, \mu')
ylabel('Power (dB)')
set(gca,'FontSize',16,'FontName','Times New Roman','LineWidth',1.2)  
xticklabels(strrep(xticklabels,'-','−'))
yticklabels(strrep(yticklabels,'-','−'))
xlim([-100 100])
ylim([-300 0])
% xticks([-100 -50 0 50 100])
% xticklabels('-100','-50','0','50','100')
grid on
% hold on
% load cold_cavity.mat
