% MONTE-CARLO SIMULATION OF EXTERNAL PHASE MODULATION WITH THERMOREFRACTIVE NOISE
clear
close all
load ./data/intextpm100ghz_3.mat
% Parameters
N       = 512;                 
dt      = 1e-2;                         
count   = 1;      
tph     = 8.23e-10;                   
kappa   = 1/tph;                      
alpha0  = 3.2;                         
beta    = -0.2;                       
F       = 2.2;                           
gamma   = 1.0;                           
deltaM  = 0.65;                          
gamma_T = 2*pi*1e5;                      
R       = 225e-6;                       
L       = 2*pi*R;                        
c       = 2.998e8;                       
naxis   = (-N/2:N/2-1).';                 
dtheta  = 2*pi/N;                     
theta   = naxis*dtheta;                  
dmu     = 1;                             
mu      = fftshift(dmu*naxis); 
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
    theta0        = circmean;
    psi_0         = psi_prev;               
    for ix = 1:N_steps        
        T0         = 300;                                        
        A          = 7.39e-5;   
        eta_nu     = 2.45e-5;                                    
        m          = 2000;
        noise(ix+1)= noise(ix) - dt*tph*gamma_T*noise(ix) + sqrt(dt*tph*A)*randn; 
        % dT         = noise(ix+1); 
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
        theta0(ix)     = sum(I.*theta)/sum(I);              
        z              = sum(I.*exp(1i*theta));              
        circmean(ix)   = angle(z);     
        if ix>1
            dw_rep(ix+1) = (circmean(ix)-circmean(ix-1))/dt;
        end
    end    
    % Save data
    data_structure = struct('circmean',circmean);    
    save(sprintf('./data/extpmTRN030126/extpmTRN030126_%d.mat',k),'-fromstruct',data_structure);
end
