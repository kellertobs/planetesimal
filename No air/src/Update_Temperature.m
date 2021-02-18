% Update Temperature
function [Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,T_diff,dT] =...
    Update_Temperature(Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,nx1,nz1,dx,dz,...
    vx_out,vz_out,Eta_out,Eta_mid,Rho_vz,Alpha_mid,T_mid,gz,...
    Nx,Nz,k_vx,k_vz,Rho_mid,Cp_mid,Hr,T_top,dt)

Hxz = Sigxx;
%% calculate stresses and shear/adiabatic heating
for j = 1:1:nx1
for i = 1:1:nz1
    Epsxz(i,j)  = 1/2 * ((vx_out(i+1,j) - vx_out(i,j))/dz...
                       + (vz_out(i,j+1) - vz_out(i,j))/dx);   %shear strain on normal grid
    Sigxz(i,j)  =       2*Eta_out(i,j)*Epsxz(i,j);          %deviatoric stress on normal grid
end 
end
for j = 2:1:nx1 %calculate normal stress/strain components on the middle of grid/pressure nodes
for i = 2:1:nz1
    Epsxx(i,j)  = (vx_out(i,j) - vx_out(i,j-1))/dx; %Normal strain rate
    Sigxx(i,j)  = 2*Eta_mid(i,j)*Epsxx(i,j); %deviatoric stress
end
end
    
%compute shear heating and adiabatic heating
for j = 2:1:nx1 
for i = 2:1:nz1
    Hxz(i,j)         = (Epsxz(i-1,j-1)*Sigxz(i-1,j-1)    + Epsxz(i-1,j)*Sigxz(i-1,j) +...
                   Epsxz(i-1,j)*Sigxz(i-1,j)        + Epsxz(i,j)*Sigxz(i,j))/4; %average of xz products
    Hs(i,j)     = 2*Hxz(i,j) + 2*Epsxx(i,j)*Sigxx(i,j); %shear heating
    % now compute adiabatic heating
    Ha(i,j)     = (vz_out(i,j) + vz_out(i-1,j))/2 * (Rho_vz(i,j) + Rho_vz(i-1,j))/2 ...
        *T_mid(i,j)*gz*Alpha_mid(i,j);
end
end
        
%% Solve temperature diffusion
%setip matrix and vector for implicit solution
NP          = Nx*Nz; % total number of P nodes to solve + ghost nodes
T_diff      = T_mid;
ind         = reshape(1:NP,Nz,Nx);


II  = [];
JJ  = [];
AA  = [];
IR  = [];
RR  = [];

%% fill implicit matrices



% internal points
%auxilary array
ii = ind(2:nz1,2:nx1); % current node
%left of current
jj = ind(2:nz1,1:nx1-1); kk = ind(2:nz1,1:nx1-1);
II = [II, ii(:)'];   JJ = [JJ, jj(:)'];   AA = [AA, -k_vx(kk(:)')/dx^2];
%right of current
jj = ind(2:nz1,3:nx1+1); kk = ind(2:nz1,2:nx1);
II = [II, ii(:)'];   JJ = [JJ, jj(:)'];   AA = [AA, -k_vx(kk(:)')/dx^2];
%above current
jj = ind(1:nz1-1,2:nx1); kk = ind(1:nz1-1,2:nx1);
II = [II, ii(:)'];   JJ = [JJ, jj(:)'];   AA = [AA, -k_vz(kk(:)')/dz^2];
%below current
jj = ind(3:nz1+1,2:nx1); kk = ind(2:nz1,2:nx1);
II = [II, ii(:)'];   JJ = [JJ, jj(:)'];   AA = [AA, -k_vz(kk(:)')/dz^2];
%current node
RhoCp_mid = Rho_mid(2:nz1,2:nx1).*Cp_mid(2:nz1,2:nx1);
Asum = (RhoCp_mid./dt)...
    +((k_vx(2:nz1,2:nx1) + k_vx(2:nz1,1:nx1-1))/dx^2)...
    +((k_vz(2:nz1,2:nx1) + k_vz(1:nz1-1,2:nx1))/dz^2);
jj = ii;
II = [II, ii(:)'];   JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'];
%RHS
RRSUM = (Hs(2:nz1,2:nx1)+Ha(2:nz1,2:nx1)+Hr(2:nz1,2:nx1))+RhoCp_mid.*T_mid(2:nz1,2:nx1)./dt;
IR = [IR, ii(:)']; RR = [RR, RRSUM(:)'];

% boundary points
%top constant temp
ii  = ind(1,2:nx1); jj  = ind(1,2:nx1); jj1 = ind(2,2:nx1); %below top boundary
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];  AA = [AA, Asum(:)'+1];
%RHS
IR = [IR, ii(:)']; RR = [RR, Asum(:)'*T_top*2];
    
%bottom constant temp
ii  = ind(Nz,2:nx1); jj  = ind(Nz,2:nx1); jj1 = ind(nz1,2:nx1); %below top boundary
Asum = zeros(size(ii)); Asum1 = Asum.*T_mid(Nz,2:nx1).*2;
II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];  AA = [AA, Asum(:)'+1];
%RHS
IR = [IR, ii(:)']; RR = [RR, Asum1(:)'];

%left insulating boundary
ii  = ind(1:Nz,1); jj  = ind(1:Nz,1); jj1 = ind(1:Nz,2); %below top boundary
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];  AA = [AA, Asum(:)'-1];
%RHS
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

%right insulating boundary
ii  = ind(1:Nz,Nx); jj  = ind(1:Nz,Nx); jj1 = ind(1:Nz,nx1); %below top boundary
Asum = zeros(size(ii));
II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
II = [II, ii(:)']; JJ = [JJ, jj1(:)'];  AA = [AA, Asum(:)'-1];
%RHS
IR = [IR, ii(:)']; RR = [RR, Asum(:)'];

%         if (j==1)
%             II = [II, ind(i,j)]; JJ = [JJ, ind(i,j)];   AA = [AA, 1];
%             II = [II, ind(i,j)]; JJ = [JJ, ind(i,j+1)]; AA = [AA, -1];
%             IR = [IR, ind(i,j)]; RR = [RR, 0];

A       = sparse(II,JJ,AA,NP,NP);
RHS     = sparse(IR,ones(size(IR)),RR,NP,1);

X           =  sqrt(abs(diag(A)));
X           =  diag(sparse(1./X));

A           =  X*A*X;
RHS         =  X*RHS;

%% output T grid
Tvec = X*(A\RHS); % backslash operator to output diffused temperature vector
for j=1:1:Nx
for i=1:1:Nz    
    % Reload solution
    T_diff(i,j)=Tvec(ind(i,j));
end
end

% compute dT
dT          = T_diff-T_mid;

% thermal timestepping conditions
maxdT       = max(max(abs(dT)));

% legacy loop for A matrix filler
% if titer<2 && maxdT>dTmax
%     dt      = dt/maxdT*dTmax;
% else
%     break
% end
% for i = 1:1:Nz
% for j = 1:1:Nx    
%     % setup boundary conditions
%     if (i == 1 || i == Nz || j==1 || j==Nx)
%         % Top BC, constant temperature
%         if (i==1 && j>1 && j<Nx)
%             II = [II, ind(i,j)]; JJ = [JJ, ind(i,j)];   AA = [AA, 1];
%             II = [II, ind(i,j)]; JJ = [JJ, ind(i+1,j)]; AA = [AA, 1];
%             IR = [IR, ind(i,j)]; RR = [RR, T_top*2];
%         end
%         % Bottom BC, constant temperature    
%         if (i==Nz && j>1 && j<Nx)
%             II = [II, ind(i,j)]; JJ = [JJ, ind(i,j)];   AA = [AA, 1];
%             II = [II, ind(i,j)]; JJ = [JJ, ind(i-1,j)]; AA = [AA, 1];
%             IR = [IR, ind(i,j)]; RR = [RR, T_mid(i,j)*2];
%             
%         end
%         % Left BC, insulating boundary
%         if (j==1)
%             II = [II, ind(i,j)]; JJ = [JJ, ind(i,j)];   AA = [AA, 1];
%             II = [II, ind(i,j)]; JJ = [JJ, ind(i,j+1)]; AA = [AA, -1];
%             IR = [IR, ind(i,j)]; RR = [RR, 0];
%             
%         end
%         % Right BC, Insulating boundary
%         if (j==Nx)
%             II = [II, ind(i,j)]; JJ = [JJ, ind(i,j)];   AA = [AA, 1];
%             II = [II, ind(i,j)]; JJ = [JJ, ind(i,j-1)]; AA = [AA, -1];
%             IR = [IR, ind(i,j)]; RR = [RR, 0];
%             
%         end    
%     else   
%     % internal points
%     % fill A matrix
%     RhoCp_mid = Rho_mid(i,j)*Cp_mid(i,j);
%     Asum = (RhoCp_mid/dt)+((k_vx(i,j) + k_vx(i,j-1))/dx^2)+((k_vz(i,j) + k_vz(i-1,j))/dz^2);
%     
% %     II = [II, ind(i,j)]; JJ = [JJ, ind(i,j-1)];   AA = [AA, -k_vx(i,j-1)/dx^2]; % left of current node
% %     II = [II, ind(i,j)]; JJ = [JJ, ind(i,j+1)];   AA = [AA, -k_vx(i,j)  /dx^2]; % right of current node
% %     II = [II, ind(i,j)]; JJ = [JJ, ind(i-1,j)];   AA = [AA, -k_vz(i-1,j)/dz^2]; % above current node
% %     II = [II, ind(i,j)]; JJ = [JJ, ind(i+1,j)];   AA = [AA, -k_vz(i,j)  /dz^2]; % below current node
%     II = [II, ind(i,j)]; JJ = [JJ, ind(i,j)];     AA = [AA, Asum];
% 
%                                                               % current node
%     %RHS
%     IR = [IR, ind(i,j)]; RR = [RR, (Hs(i,j)+Ha(i,j)+Hr(i,j))+RhoCp_mid*T_mid(i,j)/dt];
%                             
%     end
% end
% end


            
        