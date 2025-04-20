% This solves the LLE with fsolve while using finite difference method to
% calculate the second derivatives.

classdef LLE_fdm
    
    % #####################################################################
    properties
        alpha
        beta
        gamma
        F
        dx
        mu 
        N
    end
    
    % #####################################################################
    methods
        
        % #################################################################
        function self = LLE_fdm(alpha,beta,gamma,F,dx,mu,N)
            self.alpha = alpha;
            self.beta  = beta;
            self.gamma = gamma;
            self.F     = F;
            self.dx    = dx;
            self.mu    = mu;            
            self.N     = N;
        end

        % #################################################################
        function y = der2(self)                        
            a = spdiags([1 -2 1],-1:1,self.N,self.N);                        
            a(1) = -1;
            a(end) = -1;
            y = a./(self.dx^2);            
        end

        % #################################################################
        function [psi,J] = findroots(self,psi_0)

            % Extract real and imaginary parts from initial guess            
            psi_r   = psi_0(1:self.N);
            psi_m   = psi_0(self.N+1:2*self.N);
            
            % Equation-1                     
            F1      = - psi_r + self.alpha.*psi_m + (self.beta/2).*self.der2*psi_m ... 
                      - psi_r.^2.*psi_m.*self.gamma - psi_m.^3.*self.gamma + self.F;
            
            % Equation-2
            F2       = - psi_m - self.alpha.*psi_r - (self.beta/2).*self.der2*psi_r ...
                      + psi_r.^3.*self.gamma + psi_m.^2.*psi_r.*self.gamma;
            
            % Final solution containing both real and imaginary parts
            psi     = [F1;
                       F2];

            % Find Jacobian
            J       = self.jacobian(psi_0);
        end

        % #################################################################
        function [J] = jacobian(self,psi_0)           
            
            psi_r = psi_0(1:self.N,1);
            psi_m = psi_0(self.N+1:2*self.N,1); 

            e   = eye(self.N);
           
            J11 = - e - 2*diag(psi_r)*diag(psi_m)*self.gamma;

            J12 = self.alpha*e + diag(self.beta/2)*self.der2 ...
                  - diag(psi_r.^2)*self.gamma - 3*diag(psi_m.^2)*self.gamma;
            
            J21 = - self.alpha*e - diag(self.beta/2)*self.der2 ...
                  + 3.*diag(psi_r.^2)*self.gamma + diag(psi_m.^2)*self.gamma;
                        
            J22 = - e + 2*diag(psi_m)*diag(psi_r)*self.gamma;

            J   = [J11 J12;
                   J21 J22];
        end
    end
end