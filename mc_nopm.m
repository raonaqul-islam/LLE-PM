%% MONTE-CARLO SIMULATION FOR NO PHASE-MODULATION OF PUMP
% This routine simulates the system where there is no phase modulation but
% there is "random" noise in the system.

% ------------------------------+
% Raonaqul Islam, UMBC          |
% Date started: March 19, 2025  |
% Last updated: Nov 01, 2025    |
% ------------------------------+

clear
close all

tic

% Load previous solution
load ./data/intnopm100ghz_3.mat

% Monte-Carlo parameters
plot_evo     = false;
N_sim        = 60;                     % No. of Monte-Carlo simulations (multiples of 64 because 64 parallel processes will be used)
N_steps      = 6e6;                    % No. of steps in each iteration
savestep     = 1;                      % Save every "S" steps to reduce memory usage

for k = 1:N_sim

    % Use different random seeds for each simulation
    rng(k);

    disp(k);

    % *********************************************************************
    %                        SIMULATION PARAMETERS                        %
    % *********************************************************************
    N          = 512;                             % No. of points on axis | Spatial resolution
    dt         = 1e-2;                            % Temporal step size
    psi_0      = psi_prev;                        % Initial condition from previous solution
    theta0     = zeros(N_steps/savestep,1);       % Initiate peak position array
    circmean   = zeros(N_steps/savestep,1);       % Initiate circular mean array
    dw_rep = zeros(N_steps,1);
    count      = 1;                               % Sets the index of saving-to-file

    % *********************************************************************
    %                          NOISE PARAMETERS                           %
    % *********************************************************************
    A          = 10;                              % Noise amplitude    
    noise      = zeros(N_steps+1,1);              % Noise array

    % Generate figure to plot on with a white background
    if plot_evo == true
        figure('Units','inches','Position',[6 4 8 4],'Color','White');
    end

    % Time steps
    for ix = 1:N_steps

        % disp(ix)

        % *****************************************************************
        %                           PARAMETERS                            %
        % *****************************************************************
        alpha   = 3.2;                            % Detuning
        beta    = -0.2;                           % Dispersion, beta2 only in this case
        F       = 2.2;                            % Pump amplitude
        gamma   = 1.0;                            % Normalized nonlinear coefficient
        deltaM  = 0;                              % Phase modulation        
        gamma_n = 1e6;                            % Noise disspation factor
        tph     = 8.23e-6;                        % Photon lifetime
        R       = 225e-6;                         % Cavity radius

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

        if plot_evo==true && rem(ix,1/dt)==0 % Plot every photon lifetime
            plot(theta/pi,abs(psi_out).^2,'LineWidth',2.5,'Color',[0.07 0.62 1]);
            xlabel('\theta/\pi');
            ylabel('|\psi|^2');
            title(sprintf('No. of {\\it{\\tau}}_p = %0.3f/%d',ix*dt,dt*N_steps));
            xlim([-1 1]);
            ax = gca;
            set(ax,'FontSize',16,'FontName','Times New Roman','LineWidth',1.2)
            xticklabels(strrep(xticklabels,'-',char(0x2212)))
            xline(0,'LineWidth',2,'Color','k','LineStyle','--')
            grid on
            ax.GridLineWidth = 3;
            ax.GridAlpha = 0.02;
            ax.GridLineStyle = '-';
            drawnow
        end

        % Calculate centroid  
        if rem(ix,savestep)==0
            I               = abs(psi_out).^2;               % Intensity
            theta0(count)   = sum(I.*theta)/sum(I);          % Centroid
            circmean(count) = angle(sum(I.*exp(1i*theta)));  % Weighted circular mean          
            if ix>1
                dw_rep(count) = (theta0(count)-theta0(count-1))./dt;
            end
            count           = count+1;                       % Index count
        end
    end

    % *********************************************************************
    %                             SAVE DATA                               %
    % *********************************************************************        
    % Form a structure to store data (because of the parallel pool)
    data_structure = struct('theta0',theta0,'circmean',circmean);                           
    % Save the structure   
    save(sprintf('./data/mc_nopm_110325/mc_nopm_110325_%d.mat',k),'-fromstruct',data_structure); 

end

time_elapsed = toc;

% Show elapsed time
disp(['Elapsed time = ' num2str(time_elapsed)]);
