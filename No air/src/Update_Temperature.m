% Update Temperature
function [Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,T_diff,dT,dt] =...
    Update_Temperature(Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,nx1,nz1,dx,dz,...
    vx_out,vz_out,Eta_out,Eta_mid,Rho_vz,Alpha_mid,T_mid,gz,...
    Nx,Nz,k_vx,k_vz,RhoCp_mid,Hr,Number,T_top,dt,dTmax)

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
A           = sparse(Nz*Nx,Nz*Nx);
RHS         = zeros(Nz*Nx,1);


for titer=1:1:2
%% fill implicit matrices
for i = 1:1:Nz
for j = 1:1:Nx    
    % setup boundary conditions
    if (i == 1 || i == Nz || j==1 || j==Nx)
        % Top BC, constant temperature
        if (i==1 && j>1 && j<Nx)
            A(ind(i,j),ind(i,j))    = 1;
            A(ind(i,j),ind(i+1,j))  = 1;
            RHS(ind(i,j))           = T_top*2;
        end
        % Bottom BC, constant temperature    
        if (i==Nz && j>1 && j<Nx)
            A(ind(i,j),ind(i,j))    = 1;
            A(ind(i,j),ind(i-1,j))  = 1;
            RHS(ind(i,j))           = T_mid(i,j)*2;
        end
        % Left BC, insulating boundary
        if (j==1)
            A(ind(i,j),ind(i,j))    = 1;
            A(ind(i,j),ind(i,j+1))  = -1;
            RHS(ind(i,j))           = 0;
        end
        % Right BC, Insulating boundary
        if (j==Nx)
            A(ind(i,j),ind(i,j))    = 1;
            A(ind(i,j),ind(i,j-1))  = -1;
            RHS(ind(i,j))           = 0;
        end    
    else   
    % internal points
    % fill A matrix
    A(ind(i,j),ind(i,j-1))      = -k_vx(i,j-1)/dx^2;     % left of current node
    A(ind(i,j),ind(i,j+1))      = -k_vx(i,j)  /dx^2;     % right of current node
    A(ind(i,j),ind(i-1,j))      = -k_vz(i-1,j)/dz^2;     % above current node
    A(ind(i,j),ind(i+1,j))      = -k_vz(i,j)  /dz^2;     % below current node
    A(ind(i,j),ind(i,j))        = (RhoCp_mid(i,j)/dt) + ...
                                ((k_vx(i,j) + k_vx(i,j-1))/dx^2) + ...
                                ((k_vz(i,j) + k_vz(i-1,j))/dz^2); %current node
                            
    % fill RHS vector
    RHS(ind(i,j))               = (Hs(i,j) + Ha(i,j) + Hr(i,j)) +...
                                   RhoCp_mid(i,j)*T_mid(i,j)/dt;
    end
end
end

%% output T grid
Tvec = A\RHS; % backslash operator to output diffused temperature vector
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
if titer<2 && maxdT>dTmax
    dt      = dt/maxdT*dTmax;
else
    break
end
end
titer


            
        