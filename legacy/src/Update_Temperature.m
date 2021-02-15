% Update Temperature
function [Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,T_diff,dT,dt] =...
    Update_Temperature(Epsxz,Sigxz,Epsxx,Sigxx,Hs,Ha,nx,nz,dx,dz,...
    vx_out,vz_out,Eta_out,Eta_mid,Rho_vz,Alpha_mid,T_mid,gz,...
    nx1,nz1,k_vx,k_vz,RhoCp_mid,Hr,...
    Number,T_top,T_air,dt,dTmax,Material)

%% calculate stresses and shear/adiabatic heating
for j = 1:1:nx
for i = 1:1:nz
    Epsxz(i,j)  = 1/2 * ((vx_out(i+1,j) - vx_out(i,j))/dz...
                       + (vz_out(i,j+1) - vz_out(i,j))/dx);   %shear strain on normal grid
    Sigxz(i,j)  =       2*Eta_out(i,j)*Epsxz(i,j);          %deviatoric stress on normal grid
end 
end
for j = 2:1:nx %calculate normal stress/strain components on the middle of grid/pressure nodes
for i = 2:1:nz
    Epsxx(i,j)  = (vx_out(i,j) - vx_out(i,j-1))/dx; %Normal strain rate
    Sigxx(i,j)  = 2*Eta_mid(i,j)*Epsxx(i,j); %deviatoric stress
end
end
    
%compute shear heating and adiabatic heating
for j = 2:1:nx 
for i = 2:1:nz
    Hxz         = (Epsxz(i-1,j-1)*Sigxz(i-1,j-1)    + Epsxz(i-1,j)*Sigxz(i-1,j) +...
                   Epsxz(i-1,j)*Sigxz(i-1,j)        + Epsxz(i,j)*Sigxz(i,j))/4; %average of xz products
    Hs(i,j)     = 2*Hxz + 2*Epsxx(i,j)*Sigxx(i,j); %shear heating
    % now compute adiabatic heating
    Ha(i,j)     = (vz_out(i,j) + vz_out(i-1,j))/2 * (Rho_vz(i,j) + Rho_vz(i-1,j))/2 ...
        *T_mid(i,j)*gz*Alpha_mid(i,j);
end
end
        
%% Solve temperature diffusion
%setip matrix and vector for implicit solution
T_diff      = T_mid;
Ind         = Number;
A           = sparse(nz1*nx1,nz1*nx1);
RHS         = zeros(nz1*nx1,1);


for titer=1:1:2
%% fill implicit matrices
for i = 1:1:nz1
for j = 1:1:nx1    
    % setup boundary conditions
    if (i == 1 || i == nz1 || j==1 || j==nx1)
        % Top BC, constant temperature
        if (i==1 && j>1 && j<nx1)
            A(Ind(i,j),Ind(i,j))    = 1;
            A(Ind(i,j),Ind(i+1,j))  = 1;
            RHS(Ind(i,j))           = T_top*2;
        end
        % Bottom BC, constant temperature    
        if (i==nz1 && j>1 && j<nx1)
            A(Ind(i,j),Ind(i,j))    = 1;
            A(Ind(i,j),Ind(i-1,j))  = 1;
            RHS(Ind(i,j))           = T_mid(i,j)*2;
        end
        % Left BC, far field boundary
        if (j==1)
            A(Ind(i,j),Ind(i,j))    = 1;
            A(Ind(i,j),Ind(i,j+1))  = -1;
            RHS(Ind(i,j))           = 0;
        end
        % Right BC, far field boundary
        if (j==nx1)
            A(Ind(i,j),Ind(i,j))    = 1;
            A(Ind(i,j),Ind(i,j-1))  = -1;
            RHS(Ind(i,j))           = 0;
        end    
    else   
    % internal points
    % fill A matrix
    A(Ind(i,j),Ind(i,j-1))      = -k_vx(i,j-1)/dx^2;     % left of current node
    A(Ind(i,j),Ind(i,j+1))      = -k_vx(i,j)  /dx^2;     % right of current node
    A(Ind(i,j),Ind(i-1,j))      = -k_vz(i-1,j)/dz^2;     % above current node
    A(Ind(i,j),Ind(i+1,j))      = -k_vz(i,j)  /dz^2;     % below current node
    A(Ind(i,j),Ind(i,j))        = (RhoCp_mid(i,j)/dt) + ...
                                ((k_vx(i,j) + k_vx(i,j-1))/dx^2) + ...
                                ((k_vz(i,j) + k_vz(i-1,j))/dz^2); %current node
                            
    % fill RHS vector
    RHS(Ind(i,j))               = (Hs(i,j) + Ha(i,j) + Hr(i,j)) +...
                                   RhoCp_mid(i,j)*T_mid(i,j)/dt;
    end
end
end

%% output T grid
Tvec = A\RHS; % backslash operator to output diffused temperature vector
for j=1:1:nx1
for i=1:1:nz1    
    % Reload solution
    T_diff(i,j)=Tvec(Ind(i,j));
end
end

% Overwrite air temperature
% indair = find(Material == 2);
T_diff(find(Material == 2)) = T_air;

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


            
        