function Tm = Update_T_markers(Tm,T_diff,dT,marknum,nx1,nz1,dx,dz,xm,zm,xp,zp,ti)

for m = 1:1:marknum
    % define i,j indeces for upper-left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((zm(m)-zp(1))/dz)+1;
    if(j<2)
        j   =2;
    elseif(j>nx1)
        j   =nx1;
    end
    if(i<2)
        i   =2;
    elseif(i>nz1)
        i   =nz1;
    end
    
    % compute distances
    dxm1    = xm(m) - xp(j);
    dzm1    = zm(m) - zp(i);
    
    % compute weights
    wtij    = (1 - dxm1/dx) * (1 - dzm1/dz);    % topleft node
    wtij1   =      dxm1/dx  * (1 - dzm1/dz);    % topright node
    wti1j   = (1 - dxm1/dx) *      dzm1/dz ;    % bottomleft node
    wti1j1  =      dxm1/dx  *      dzm1/dz;     % bottomright node
    
    % update marker properties
    if ti == 1
        % interpolate the diffused temperature into markers for first
        % timestep
        Tm(m)   =   T_diff(i,j)*wtij    + T_diff(i+1,j)*wti1j +...
                    T_diff(i,j+1)*wtij1 + T_diff(i+1,j+1)*wti1j1;  
    else
        % interpolate dT into markers and add onto old temperature markers
        dTm     =   dT(i,j)*wtij        + dT(i+1,j)*wti1j     +...
                    dT(i,j+1)*wtij1     + dT(i+1,j+1)*wti1j1;
        Tm(m)   =   Tm(m) + dTm;
    end
end
