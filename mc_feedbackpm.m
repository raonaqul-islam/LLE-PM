%% MONTE-CARLO SIMULATION FOR FEEDBACK PHASE-MODULATION OF PUMP
% This routine simulates when there is noise in the system that adds up in
% every time step, but there is also external phase modulation using
% external RF source.

% ------------------------------+
% Raonaqul Islam, UMBC          |
% Date started: March 19, 2025  |
% Last updated: Sep 04, 2025    |
% ------------------------------+

tic

clear
close all

% Initial condition from previous solution
load ./data/intextpm2.mat

% Monte-Carlo parameters
N_sim         = 300;                       % No. of Monte-Carlo simulations
% ------------------------
% No. of steps in each iteration
% NOTE: Since there is a delay in the feedback reaching the modulator,
% during that period, there is no modulation. Hence, the single soliton
% undergoes some chaotic behavior and takes around 15 photon lifetimes to
% settle down. That is why we simulate for extra 100 photon lifetimes and
% ignore the first 100 in our PSD calculation. 100 photon lifetimes is 100*1e4
% = 1e6 steps. So the total steps become 5e6+1e6=6e6 steps.
% (~Raonaq, 8/27/25)
N_steps       = 6e6;                     
% -------------------------
fbackpmtheta0 = zeros(N_steps,N_sim);      % Store peaks

parfor k = 1:N_sim

    disp(k)

    % Use different random seeds for each iteration
    rng(k);

    % *********************************************************************
    %                        SIMULATION PARAMETERS                        %
    % *********************************************************************
    N          = 512;                             % No. of points on axis | Spatial resolution
    dt         = 1e-4;                            % Temporal step size
    psi_0      = psi_prev;                        % Initial condition from previous solution
    theta0     = zeros(N_steps,1);                % Initiate peak position array
    circmean   = zeros(N_steps,1);                % Initiate circular mean array    

    % *********************************************************************
    %                            NOISE PARAMETERS                         %
    % *********************************************************************
    A          = 10;                             % Noise amplitude    
    noise      = zeros(N_steps+1,1);              % Noise array
    nprev      = 0;                               % Initialize

    % Time steps
    for ix = 1:N_steps

        % *****************************************************************
        %                           PARAMETERS                            %
        % *****************************************************************
        tph     = 2.5e-6;                         % Photon-lifetime
        alpha   = 2.0;                            % Detuning
        beta    = -0.2;                           % Dispersion, beta2 only in this case
        F       = 1.4103;                         % Pump amplitude
        gamma   = 1.0;                            % Normalized nonlinear coefficient
        deltaM  = 0.65;                           % Phase modulation
        gamma_n = 2e5;                              % Noise disspation factor


        % *****************************************************************
        %                         DISCRETIZATION                          %
        % *****************************************************************
        naxis   = (-N/2:N/2-1).';                 % General axis
        dtheta  = 2*pi/N;                         % Spatial step-size
        theta   = naxis*dtheta;                   % Spatial (Azimuthal) domain
        dmu     = 1;                              % Mode number domain step-size
        mu      = fftshift(dmu*naxis);            % Mode number domain

        % *****************************************************************
        %                        NOISE CALCULATION                        %
        % *****************************************************************
    
        % Add random noise
        noise(ix+1) = noise(ix) - dt*tph*gamma_n*noise(ix) + sqrt(dt*tph*A)*randn;
        
        % Add time delay for feedback
        delay = 1/dt;                                      % For now, delay added for 1 photon-lifetime
        if ix < delay 
            % Wait for some time to get the feedback        
            alphaM  = noise(ix+1);                         % Feedback off
            F_tilde = fft(F*ones(1,N)).';                  % No modulation            
        else 
            % Starts feedback and modulation
            alphaM  = noise(ix+1) - noise(ix+1-delay);     % Feedback on
            F_tilde = fft(F*exp(1i*deltaM*sin(theta)));    % Modulation on
        end

        % *****************************************************************
        %                            SPLIT STEP                           %
        % *****************************************************************

        % Nonlinear part (half-step)
        psi_nl      = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));

        % Linear part (full-step)
        psi_0_tilde = fft(psi_nl);
        A_tilde     = -(1+1i*alpha) - 1i*alphaM*mu + 1i.*beta.*(mu).^2;
        psi_l_tilde = (psi_0_tilde+F_tilde./A_tilde).*exp(A_tilde*dt)-F_tilde./A_tilde;
        psi_l       = ifft(psi_l_tilde);

        % Nonlinear part (half-step)
        psi_0       = psi_l;
        psi_nl2     = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));
        psi_0       = psi_nl2;

        % Output after one-step
        psi_out     = psi_nl2;

        % Calculate centroid        
        I           = abs(psi_out).^2;             % Intensity
        theta0(ix)  = sum(I.*theta)/sum(I);        % Center of mass
        circmean(ix) = angle(sum(I.*exp(1i*theta))); % Weighted circular mean

    end

    % Store this iteration
    fbackpmtheta0(:,k) = theta0;
    fbackpmcircmean(:,k) = circmean;

end

time_elapsed = toc;

% Save the results in a mat file (v7.3 for larger files)
save('./data/mc_fbackpm_090425.mat','fbackpmtheta0','fbackpmcircmean','time_elapsed','-v7.3');