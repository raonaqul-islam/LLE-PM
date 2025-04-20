% LLE stability analysis with fsolve. This version considers input phase
% modulation. Important to note that it uses the solution from split-step
% time evolution solver. Therefore, to obtain solution from this code, you
% need to run the time-evolution first, save the "psi_out" as "sol0.mat" in a
% folder named "data" in the code folder. 

% -----------------------------+
% Raonaqul Islam, UMBC         |
% Date started: March 22, 2025 |
% Last updated: April 20, 2025 |
% -----------------------------+

clear
close
clc

% ******************
% Input Parameters
% ******************
alpha   = 3.5;                                  % Detuning
beta    = -0.1549;                              % Dispersion, beta2 only in this case
F       = 2.3;                                  % Power
gamma   = 1;                                    % Normalized nonlinear coefficient
method  = 'fft';                                % Method for differentiation, 'fdm' or 'fft'
deltaM  = 0.65;                                 % Modulation depth

% ******************************
% Spatial domain discretization
% ******************************
N       = 512;                                  % No. of points on axis
naxis   = (-N/2:N/2-1).';                       % General axis
dtheta  = 2*pi/N;                               % Spatial step-size
theta   = naxis*dtheta;                         % Spatial domain
dmu     = 1;                                    % Mode number domain step-size
mu      = fftshift(dmu*naxis);                  % Mode number domain

% *****************
% Phase modulation
% *****************
F       = F*exp(1i*deltaM*sin(theta));          % Input power after phase modulation

% *****************
% fsolve operation
% *****************
% Initial guess
load ./data/sol0.mat                            % Using solution from time-evolution
psi_0   = psi_out;
psi_r0  = real(psi_0);                          % Real part
psi_m0  = imag(psi_0);                          % Imaginary part
psi_0   = [psi_r0;psi_m0];                      % Complete guess

% fsolve function
switch method
    case 'fdm'
        fun     = LLE_fdm(alpha,beta,gamma,F,dtheta,mu,N); 
    case 'fft'
        fun     = LLE_fft(alpha,beta,gamma,F,dtheta,mu,N);
end

% Verify the jacobian matrix
fun.verify_jacobian(psi_0,0); % 1: Plot verification; 0: Do not plot verification

% fsolve options
options = optimset('Display','iter','Jacobian','on','TolFun',1e-10,'TolX',...
          1e-10,'Algorithm','levenberg-marquardt','ScaleProblem','Jacobian',...
          'MaxIter',100);

% find roots using fsolve
[psi,fval,exitflag,output,Jacobian] = fsolve(@fun.findroots,psi_0,options);    

% **************************************
% Extract and form the complex solution
% **************************************
psi_r   = psi(1:N);
psi_m   = psi(N+1:end);
psi_out = psi_r + 1i*psi_m;

% *************
% Plot Soliton
% *************
figure('Units','Inches','Position',[2 2 12 8]);
subplot(211)
plt1    = plot(theta./pi,abs(psi_out),'Color',[0.07 0.62 1]);
xtext   = '\theta/\pi';
ytext   = '|\Psi|';
xline(0.5,'LineWidth',2,'Color','Black','LineStyle','--');
customplot(plt1,xtext,ytext);

% *******************
% Stability Analysis
% *******************
eigenvalues = eig(Jacobian); % Extracting the eigenvalues of the Jacobian matrix

% Plot eigenvalues
subplot(212)
plt2     = scatter(real(eigenvalues),imag(eigenvalues),'MarkerEdgeColor',[0.07 0.62 1],'MarkerFaceColor','None');
xtext    = 'Real';
ytext    = 'Imaginary';
xline(0,'LineWidth',3,'LineStyle','--','Color','k');
customplot(plt2,xtext,ytext);