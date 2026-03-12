% SIMULATION FOR FEEDBACK PHASE-MODULATION OF PUMP
%
% ------------------------------+
% Raonaqul Islam, UMBC          |
% Date started: March 19, 2025  |
% Last updated: Nov 21, 2025    |
% ------------------------------+
% Feb 15, 26: Reduce A to check whether the ripples are visible or not
% That would confirm that when the noise is too low, we cannot numerically
% observe the ripples

tic

clc
clear
close all

% Initial condition from previous solution
load ./data/intextpm100ghz_3.mat

% Monte-Carlo parameters
N_sim         = 120;                        % No. of Monte-Carlo simulations
N_steps       = 1e8;                       % No. of time steps
savestep      = 1;                         % Save every "S" steps to reduce memory usage
delay         = 1e4;                       % Delay in time-steps
skip          = 1e4;                       % Skip a few time-steps in the beginning

% Do you want to plot the time evolution? (true/false)
plot_evo      = false;                  

parfor k = 1:N_sim
    
    disp(k);
    
    % Use different random seeds for each iteration
    rng(k);

    % *********************************************************************
    %                        SIMULATION PARAMETERS                        %
    % *********************************************************************
    N          = 512;                             % No. of points on axis | Spatial resolution
    dt         = 1e-2;                            % Temporal step size
    psi_0      = psi_prev;                        % Initial condition from previous solution
    theta0     = zeros(N_steps/savestep,1);       % Initiate peak position array
    circmean   = zeros(N_steps/savestep,1);       % Initiate circular mean array
    phi_un     = zeros(N_steps/savestep,1);
    dw_rep     = zeros(N_steps/savestep,1);       % Initiate w_rep noise array
    count      = 1;                               % Sets the index of saving-to-file
    R          = 225e-6;                          % Cavity radius
    deltaM     = 0.65;                            % Phase modulation
    
    % *********************************************************************
    %                         NOISE PARAMETERS                            %
    % *********************************************************************
    A          = 10;                         % Noise amplitude
    noise      = zeros(N_steps+1,1);              % Noise array
    
    % Generate figure to plot on with a white background
    if plot_evo == true
        figure('Units','inches','Position',[6 4 8 4],'Color','White');
    end

    % Time steps
    for ix = 1:N_steps

        % *****************************************************************
        %                           PARAMETERS                            %
        % *****************************************************************
        tph     = 8.23e-10;                       % Photon-lifetime
        kappa   = 1/tph;                          % Cavity loss
        alpha   = 3.2;                            % Detuning
        beta    = -0.2;                           % Dispersion, beta2 only in this case
        F       = 2.2;                            % Pump amplitude
        gamma   = 1.0;                            % Normalized nonlinear coefficient        
        gamma_n = 1e6;                            % Noise disspation factor

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

        % Add gaussian white noise
        noise(ix+1) = noise(ix) - dt*tph*gamma_n*noise(ix) + sqrt(dt*tph*A)*randn;
        
        % Add time delay for feedback        
        if ix < 2*skip
            % Wait for some time to get the feedback
            alphaM  = noise(ix+1);                           % Feedback off
            F_tilde = fft(F*ones(1,N)).';                    % No modulation
        else            
            % Modulate pump
            F_tilde = fft(F*exp(1i*deltaM*sin(theta)));
            % Starts feedback and modulation
            alphaM  = dw_rep(ix-delay) + noise(ix+1);
        end

        % *****************************************************************
        %                            SPLIT STEP                           %
        % *****************************************************************

        % Nonlinear part (half-step)
        psi_nl      = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));

        % Linear part (full-step)
        psi_0_tilde = fft(psi_nl);
        A_tilde     = -(1+1i*alpha) + 1i*alphaM*mu + 1i.*beta.*(mu).^2;
        psi_l_tilde = (psi_0_tilde+F_tilde./A_tilde).*exp(A_tilde*dt)-F_tilde./A_tilde;
        psi_l       = ifft(psi_l_tilde);

        % Nonlinear part (half-step)
        psi_0       = psi_l;
        psi_nl2     = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));
        psi_0       = psi_nl2;

        % Output after one-step
        psi_out     = psi_nl2;
        
        % *****************************************************************
        %                     CENTROID AND WREP NOISE
        % *****************************************************************
        
        % Skip few steps before starting to store centroid (let soliton settle down)    
        
        I              = abs(psi_out).^2;                   % Intensity
        theta0(count)     = sum(I.*theta)/sum(I);              % Center of mass
        z              = sum(I.*exp(1i*theta));             % Weighted circular mean
        circmean(count)   = angle(z);                          % Weighted mean of angle
        if ix>2*skip
            dw_rep(count) = (theta0(count)-theta0(count-1))./dt;
        end
        count = count + 1;
        

        if plot_evo==true && rem(ix,100/dt)==0 % Plot every 100 photon lifetime
            plot(theta/pi,abs(psi_out).^2,'LineWidth',2.5,'Color',[0.07 0.62 1]);
            xlabel('\theta/\pi');
            ylabel('|\psi|^2');
            title(sprintf('No. of {\\it{\\tau}}_p = %0.3f/%d',ix*dt,dt*N_steps));
            xlim([-1 1]);
            ax = gca;
            set(ax,'FontSize',16,'FontName','Times New Roman','LineWidth',1.2)
            xticklabels(strrep(xticklabels,'-',char(0x2212)))
            xline(1/2,'LineWidth',2,'Color','k','LineStyle','--')
            grid on
            ax.GridLineWidth = 3;
            ax.GridAlpha = 0.02;
            ax.GridLineStyle = '-';
            drawnow
        end

    end

    % *********************************************************************
    %                             SAVE DATA                               %
    % *********************************************************************
    % Form a structure to store data (because of the parallel pool)
    data_structure = struct('circmean',circmean);
    % Save the structure
    save(sprintf('./data/mc_fbackpm_021526/mc_fbackpm_021526_%d.mat',k),'-fromstruct',data_structure);

end

time_elapsed = toc;

% Show elapsed time
disp(['Elapsed time = ' num2str(time_elapsed)]);
