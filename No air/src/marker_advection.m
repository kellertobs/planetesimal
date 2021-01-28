function [xm,zm] = marker_advection(marknum,xm,zm,Nx,Nz,dx,dz,xvx,zvx,xvz,zvz,xp,zp,...
    vx_out,vz_out,vx_mid,vz_mid,vpratio,dt)
%advect the marker coordinates. Eq. 8.19
    % interpolate grid velocities into markers
    for m=1:1:marknum
        %compute 4th order Runge-kutta velocities
        %create temporary variables for markers
        xmRK    = xm(m);
        zmRK    = zm(m);
        for rk = 1:1:4
        %interpolate vx_mid and vz_mid
        j       = fix((xmRK-xp(1))/dx)+1;
        i       = fix((zmRK-zp(1))/dz)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Nz)
            i=Nz;
        end
        %compute distances
        dxm1    = abs(xmRK-xp(j));
        dzm1    = abs(zmRK-zp(i));
        %compute weights
        wtij    = (1-dxm1/dx)*(1-dzm1/dz);
        wtij1   = (dxm1/dx)*(1-dzm1/dz);
        wti1j   = (1-dxm1/dx)*(dzm1/dz);
        wti1j1  = dxm1*dzm1/dx/dz;
        %compute velocity on node i,j
        vxm(rk) = vx_mid(i,j)*wtij + vx_mid(i,j+1)*wtij1 + vx_mid(i+1,j)*wti1j +...
            vx_mid(i+1,j+1)*wti1j1;
        vzm(rk) = vz_mid(i,j)*wtij + vz_mid(i,j+1)*wtij1 + vz_mid(i+1,j)*wti1j +...
            vz_mid(i+1,j+1)*wti1j1;
        
        % Interpolate vx
        % Define i,j indexes for the upper left node
        j=fix((xmRK-xvx(1))/dx)+1;
        i=fix((zmRK-zvx(1))/dz)+1;
        if(j<1)
            j=1;
        elseif(j>Nx-1)
            j=Nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>Nz)
            i=Nz;
        end
        %calculate distances
        dxm1    = xmRK-xvx(j);
        dzm1    = zmRK-zvx(i);
        %compute weights of distances from real nodes
        wtij    = (1-dxm1/dx)*(1-dzm1/dz);
        wtij1   = (dxm1/dx)*(1-dzm1/dz);
        wti1j   = (1-dxm1/dx)*(dzm1/dz);
        wti1j1  = dxm1*dzm1/dx/dz;
        %compute velocity on node i,j
        vxm(rk) = vxm(rk)*vpratio + (1-vpratio)*vx_out(i,j)*wtij...
            +vx_out(i,j+1)*wtij1+vx_out(i+1,j)*wti1j+...
            vx_out(i+1,j+1)*wti1j1;
        
        % Interpolate vz
        % Define i,j indexes for the upper left node
        j       = fix((xmRK-xvz(1))/dx)+1;
        i       = fix((zmRK-zvz(1))/dz)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Nz-1)
            i=Nz-1;
        end
        %calculate distances
        dxm1    = xmRK-xvz(j);
        dzm1    = zmRK-zvz(i);
        %compute weights of distances from real nodes
        wtij    = (1-dxm1/dx)*(1-dzm1/dz);
        wtij1   = (dxm1/dx)*(1-dzm1/dz);
        wti1j   = (1-dxm1/dx)*(dzm1/dz);
        wti1j1  = dxm1*dzm1/dx/dz;
        %compute velocity on node i,j
        vzm(rk) = vzm(rk)*vpratio + (1-vpratio)*vz_out(i,j)*wtij+vz_out(i,j+1)*wtij1+vz_out(i+1,j)*wti1j+...
            vz_out(i+1,j+1)*wti1j1;
        
        %change coordinates of xm & zm for upper orders
        if (rk==1 || rk==2)
            xmRK = xm(m) + vxm(rk)*dt/2;
            zmRK = zm(m) + vzm(rk)*dt/2;
        elseif (rk==3)
            xmRK = xm(m) + vxm(rk)*dt;
            zmRK = zm(m) + vzm(rk)*dt;
        end
        end
        %compute 4th order RK effective velocity
        vxm_eff = (vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4))/6;
        vzm_eff = (vzm(1)+2*vzm(2)+2*vzm(3)+vzm(4))/6;
        %move markers
        xm(m) = xm(m)+ dt*vxm_eff;   
        zm(m) = zm(m)+ dt*vzm_eff;
    end