%% MONTE-CARLO SIMULATION FOR NO PHASE-MODULATION OF PUMP
% This routine simulates the system where there is no phase modulation but
% there is "random" noise in the system.

% ------------------------------+
% Raonaqul Islam, UMBC          |
% Date started: March 19, 2025  |
% Last updated: Sep 5, 2025     |
% ------------------------------+

clear
close all

tic

% Load previous solution
load ./data/intnopm2.mat

% Monte-Carlo parameters
N_sim        = 300;                    % No. of Monte-Carlo simulations
N_steps      = 5e6;                    % No. of steps in each iteration
nopmtheta0   = zeros(N_steps,N_sim);   % Store peaks
nopmcircmean = zeros(N_steps,N_sim);   % Store circular means

parfor k = 1:N_sim

    % Use different random seeds for each simulation
    rng(k);

    disp(k);

    % *********************************************************************
    %                        SIMULATION PARAMETERS                        %
    % *********************************************************************
    N          = 512;                             % No. of points on axis | Spatial resolution
    dt         = 1e-4;                            % Temporal step size
    psi_0      = psi_prev;                        % Initial condition from previous solution
    theta0     = zeros(N_steps,1);                % Initiate peak position array
    circmean   = zeros(N_steps,1);                % Initiate circular mean array

    % *********************************************************************
    %                          NOISE PARAMETERS                           %
    % *********************************************************************
    A          = 10;                               % Noise amplitude    
    noise      = zeros(N_steps+1,1);               % Noise array

    % Time steps
    for ix = 1:N_steps

        % *****************************************************************
        %                           PARAMETERS                            %
        % *****************************************************************
        alpha   = 2.0;                            % Detuning
        beta    = -0.2;                           % Dispersion, beta2 only in this case
        F       = 1.4103;                         % Pump amplitude
        gamma   = 1.0;                            % Normalized nonlinear coefficient
        deltaM  = 0;                              % Phase modulation        
        gamma_n = 2e5;                            % Noise disspation factor
        tph     = 2.5e-6;                         % Photon lifetime

        % *****************************************************************
        %                         DISCRETIZATION                          %
        % *****************************************************************
        naxis   = (-N/2:N/2-1).';                 % General axis
        dtheta  = 2*pi/N;                         % Spatial step-size
        theta   = naxis*dtheta;                   % Spatial (Azimuthal) domain
        dmu     = 1;                              % Mode number domain step-size
        mu      = fftshift(dmu*naxis);            % Mode number domain

        % *****************************************************************
        %                            PUMP POWER                           %
        % *****************************************************************
        F_tilde = fft(F*ones(1,N)).';             % Fourier transform of power
                                                  % including the modulation term

        % *****************************************************************
        %                            SPLIT STEP                           %
        % *****************************************************************

        % Set random noise        
        % NOTE: time-step is normalized to photon-lifetime in order to scale
        % noise with dissipation
        noise(ix+1) = noise(ix) - dt*tph*gamma_n*noise(ix) + sqrt(dt*tph*A)*randn;

        % Nonlinear part (half-step)
        psi_nl      = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));

        % Linear part (full-step)
        psi_0_tilde = fft(psi_nl);
        A_tilde     = -(1+1i*alpha) - 1i.*noise(ix+1).*mu + 1i.*beta.*(mu).^2;
        psi_l_tilde = (psi_0_tilde+F_tilde./A_tilde).*exp(A_tilde.*dt)-F_tilde./A_tilde;
        psi_l       = ifft(psi_l_tilde);

        % Nonlinear part (half-step)
        psi_0       = psi_l;
        psi_nl2     = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));
        psi_0       = psi_nl2;

        % Output after one-step
        psi_out     = psi_nl2;

        % Calculate centroid        
        I           = abs(psi_out).^2;               % Intensity
        theta0(ix)  = sum(I.*theta)/sum(I);          % Centroid
        circmean(ix) = angle(sum(I.*exp(1i*theta))); % Weighted circular mean
    end

    % Store this iteration
    nopmtheta0(:,k) = theta0;
    nopmcircmean(:,k) = circmean;

end

time_elapsed = toc;

% Save the results in a mat file (v7.3 for larger files)
save('./data/mc_nopm_090125.mat','nopmtheta0','nopmcircmean','time_elapsed','-v7.3');