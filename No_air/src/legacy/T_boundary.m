function [T_mid,T_vx,T_vz] = T_boundary(T_mid,T_vx,T_vz,Plumetype,T_top,T_bot...
    L,xp,nx,nz,nx1,nz1)
    % Apply thermal boundary condition to interpolated nodes
    % upper boundary, constant temperature 
    T_mid(1,2:nx)   = 2*T_top   - T_mid(2,2:nx);
    T_vx(1,2:nx)   = 2*T_top   - T_vx(2,2:nx);
    T_vz(1,2:nx)   = 2*T_top   - T_vz(2,2:nx);
    % lower boundary, constant temperature
    T_mid(Nz,2:nx)  = 2*(T_bot)   - T_mid(nz,2:nx);
    T_vx(Nz,2:nx)   = 2*(T_bot)   - T_vx(nz,2:nx);
    T_vz(Nz-1,2:nx) = 2*(T_bot)   - T_vz(nz,2:nx);    

    
    % left boundary, insulating boundary
    T_mid(:,1)      = T_mid(:,2);
    T_vx(:,1)       = T_vx(:,2);
    T_vz(:,1)       = T_vz(:,2);
    % right boundary, insulating boundary
    T_mid(:,Nx)     = T_mid(:,nx);
    T_vx(:,Nx-1)    = T_vx(:,nx);
    T_vz(:,Nx)      = T_vz(:,nx);
    % plume beneath lithosphere
%     midx = fix(L/2);
    
%     indexP = find(abs(xp(1,:)) <= midx+ (L*0.23) & abs(xp(1,:)) >= midx - (L*0.23));
%     T_mid(nz1:Nz,indexP) = T_bot+300;
%     indexvx = find(abs(xvx(1,:)) <= midx+ (L*0.23) & abs(xvx(1,:)) >= midx - (L*0.23));
%     T_vx(nz:Nz,indexvx) = T_bot+300;
%     indexvz = find(abs(xvz(1,:)) <= midx+ (L*0.23) & abs(xvz(1,:)) >= midx - (L*0.23));
%     T_vz(nz:Nz,indexvz) = T_bot+300;