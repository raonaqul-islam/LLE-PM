% Simulating thermorefractive noise (TRN)
clear
close all
load ./data/intextpm100ghz_3.mat
% Simulation parameters
N       = 512;                            % No. of points on axis | Spatial resolution
dt      = 1e-2;                           % Temporal step size
count   = 1;                              % Sets the index of saving-to-file
% Device parameters
tph     = 8.23e-10;                       % Photon-lifetime
kappa   = 1/tph;                          % Cavity loss
alpha0  = 3.2;                            % Detuning
beta    = -0.2;                           % Dispersion, beta2 only in this case
F       = 2.2;                            % Pump amplitude
gamma   = 1.0;                            % Normalized nonlinear coefficient
gamma_T = 2*pi*1e5;                       % Noise disspation factor
R       = 225e-6;                         % Radius of the cavity
L       = 2*pi*R;                         % Length of the cavity
c       = 2.998e8;                        % Speed of light
naxis   = (-N/2:N/2-1).';                 % General axis
dtheta  = 2*pi/N;                         % Spatial step-size
theta   = naxis*dtheta;                   % Spatial (Azimuthal) domain
dmu     = 1;                              % Mode number domain step-size
mu      = fftshift(dmu*naxis);            % Mode number domain
deltaM  = 0;                               % Forcing no-modulation
plot_evo= false;
N_sim   = 120;                       
N_steps = 1e8;                      
parfor k = 1:N_sim
    disp(k);    
    rng(k);      
    if plot_evo == true
        figure('Units','inches','Position',[6 4 8 4],'Color','White');
    end  
    % Initialize    
    circmean      = zeros(N_steps,1);               
    noise         = zeros(N_steps+1,1);   
    dw_rep        = zeros(N_steps,1);
    psi_0         = psi_prev;               
    for ix = 1:N_steps        
        T0         = 300;                                        
        A          = 7.39e-5;   
        eta_nu     = 2.45e-5;                                    
        m          = 2000;
        noise(ix+1)= noise(ix) - dt*tph*gamma_T*noise(ix) + sqrt(dt*tph*A)*randn; 
        dT         = noise(ix+1); 
        dT         = -1e-3;
        T          = T0 + dT;                                    
        wu0        = 2*pi*193e12;
        wpmp       = wu0 + (kappa/2)*alpha0;
        D1         = 2*pi*100e9;                                 
        D2         = -(kappa/2)*beta;                             
        wu_T0      = wu0 + D1.*mu + D2*mu.^2;   
        wu_T       = wu_T0 - dT*eta_nu*wu_T0.^2*L./((mu+m).*c);  
        Dint_T     = -(2/kappa)*(wu_T - wu0 - D1.*mu);  
        alpha      = -(2/kappa)*(wu_T(1) - wpmp);                 
        F_tilde    = fft(F*exp(1i*deltaM*sin(theta)));    
        % Nonlinear part (half-step)
        psi_nl      = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));
        % Linear part (full-step)
        psi_0_tilde = fft(psi_nl);
        % A_tilde     = -(1+1i*alpha) + 1i*alphaM*mu + 1i*beta.*mu.^2;
        A_tilde     = -(1+1i*alpha) + 1i*Dint_T;
        psi_l_tilde = (psi_0_tilde+F_tilde./A_tilde).*exp(A_tilde*dt)-F_tilde./A_tilde;
        psi_l       = ifft(psi_l_tilde);
        % Nonlinear part (half-step)
        psi_0       = psi_l;
        psi_nl2     = psi_0.*exp(1i.*(gamma*abs(psi_0).^2.*dt/2));
        psi_0       = psi_nl2;
        % Output after one-step
        psi_out     = psi_nl2;        
        % Centroid calculation
        I              = abs(psi_out).^2;                   
        % theta0(ix)     = sum(I.*theta)/sum(I);              
        z              = sum(I.*exp(1i*theta));              
        circmean(ix)   = angle(z);        
        if ix>1
            dw_rep(ix+1) = (circmean(ix)-circmean(ix-1))/dt;
        end
    end    
    % Save data
    data_structure = struct('circmean',circmean);    
    save(sprintf('./data/nopmTRN030126/nopmTRN030126_%d.mat',k),'-fromstruct',data_structure);
end
