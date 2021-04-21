function [T_diff,dT,dt] = Thermal_diffusion(nx1,nz1,dx,dz,T_mid,k_vx,k_vz,RhoCp_mid,Number,Hs,Ha,Hr,T_top,T_air,dt,dTmax,Material)
% Thermal iterations


%setip matrix and vector for implicit solution
T_diff      = T_mid;
Ind         = Number;
A           = sparse(nz1*nx1,nz1*nx1);
RHS         = zeros(nz1*nx1,1);
maxdT       = 1000;
titer       = 1;


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
    
        
    % internal points
    else
    
%     % overwrite air
%     if Material(i,j) == 2
%     A(Ind(i,j),Ind(i,j))    = 1;
%     A(Ind(i,j),Ind(i+1,j))  = 1;
%     RHS(Ind(i,j))           = T_top*2;  
%     end
    
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
maxdT       = max(max(abs(dT)));

% thermal timestepping conditions
maxdT       = max(max(abs(dT)));
if titer<2 && maxdT>dTmax
    dt      = dt/maxdT*dTmax;
else
    break
end
end
titer

                        