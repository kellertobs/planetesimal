% Update Temperature
adv_T0 = adv_T;     lapl_T0 = lapl_T;
T0 = T_mid;        
Hxz = Sigxx;
%% calculate stresses and shear/adiabatic heating
Epsxz(1:nz1,1:nx1) = 1/2 .* ((vx_out(2:nz2,1:nx1) - vx_out(1:nz1,1:nx1))./dz...
    +  (vz_out(1:nz1,2:nx2) - vz_out(1:nz1,1:nx1))./dx);
Sigxz(1:nz1,1:nx1) =   2 .* Eta_out(1:nz1,1:nx1).*Epsxz(1:nz1,1:nx1);

%calculate normal stress/strain components on the middle of grid/pressure nodes
Epsxx(2:nz1,2:nx1) = (vx_out(2:nz1,2:nx1) - vx_out(2:nz1,1:nx1-1))./dx; %Normal strain rate
Sigxx(2:nz1,2:nx1) = 2.*Eta_mid(2:nz1,2:nx1).*Epsxx(2:nz1,2:nx1); %deviatoric stress

% compute shear heating and adiabatic heating
Hxz(2:nz1,2:nx1) =    (Epsxz(1:nz1-1,1:nx1-1)  .*Sigxz(1:nz1-1,1:nx1-1)...
    + Epsxz(1:nz1-1,2:nx1)    .*Sigxz(1:nz1-1,2:nx1) +...
    Epsxz(1:nz1-1,2:nx1)    .*Sigxz(1:nz1-1,2:nx1)...
    + Epsxz(2:nz1,2:nx1)      .*Sigxz(2:nz1,2:nx1))./4; %average of xz products
Hs(2:nz1,2:nx1)  =  2.*Epsxx(2:nz1,2:nx1)      .*Sigxx(2:nz1,2:nx1) +...
    2.*Hxz(2:nz1,2:nx1); %shear heating
% now compute adiabatic heating
Ha(2:nz1,2:nx1)         = (vz_out(2:nz1,2:nx1) + vz_out(1:nz1-1,2:nx1))./2 .* (Rho_vz(2:nz1,2:nx1) + Rho_vz(1:nz1-1,2:nx1))./2 ...
    .*T_mid(2:nz1,2:nx1).*gz.*Alpha_mid(2:nz1,2:nx1);

%% Solve temperature diffusion
switch Tsolver
    case 'explicit'
        %assume average thermal conductivity for now
        Kappa0 = mode(k_vx(:)) / mean(Rho_mid(:)) / mode(Cp_mid(:));
        icx = 2:nx2-1; icz = 2:nz2-1;
        if ti > 1
            for it = 1:10
                lapl_T = Kappa0*(diff(T_mid(:,icx),2,1) + diff(T_mid(icz,:),2,2))/dz/dx;... % finite differences
                    + (Hr(2:end-1,2:end-1) + Hs(2:end-1,2:end-1) + Ha(2:end-1,2:end-1))./mean(Rho_mid(:))./ Cp_mantle; % + internal sources
                
%                 adv_T = advection2(vx_out,vz_out,T_mid,dx,dz,AdvRegime);
%                 dMdt  = advection2(vx_out,vz_out,Material,dx,dz,AdvRegime);
                [adv_T] = advection2(vx_out,vz_out,T_mid,dx,dz,AdvRegime);
                T_mid(2:end-1,2:end-1)  =  T0(2:end-1,2:end-1)...
                    + (theta.*(lapl_T-adv_T)+(1-theta).*(lapl_T-adv_T))*dt;
                
            end
            
        else
            lapl_T  = Kappa0*(diff(T_mid(:,icx),2,1) + diff(T_mid(icz,:),2,2))/dz/dx;... % finite differences
                + (Hr(2:end-1,2:end-1) + Hs(2:end-1,2:end-1) +...
                Ha(2:end-1,2:end-1))./mean(Rho_mid(:))./ Cp_mantle; % + internal sources 
            [adv_T] = advection2(vx_out,vz_out,T_mid,dx,dz,AdvRegime);
        end
        T_mid(1,:) = T_top;
        T_mid(end,:) = T_mid(end,:);
        T_mid(:,1) = T_mid(:,2);
        T_mid(:,end) = T_mid(:,end-1);
        
        T_out = T_mid;
        dT = T_out-T0;
    case 'implicit'
        %         dT = T_out-T_mid;
        %     case 'implicit'
        %         %setip matrix and vector for implicit solution
        %         NP          = nx2*nz2; % total number of P nodes to solve + ghost nodes
        %         T_out      = T_mid;
        %         ind         = reshape(1:NP,nz2,nx2);
        %
        %
        %         II  = [];
        %         JJ  = [];
        %         AA  = [];
        %         IR  = [];
        %         RR  = [];
        %
        %         %% fill implicit matrices
        %         % internal points
        %         %auxilary array
        %         ii = ind(2:nz1,2:nx1); % current node
        %         %left of current
        %         jj = ind(2:nz1,1:nx1-1); kk = ind(2:nz1,1:nx1-1);
        %         II = [II, ii(:)'];   JJ = [JJ, jj(:)'];   AA = [AA, -k_vx(kk(:)')/dx^2];
        %         %right of current
        %         jj = ind(2:nz1,3:nx1+1); kk = ind(2:nz1,2:nx1);
        %         II = [II, ii(:)'];   JJ = [JJ, jj(:)'];   AA = [AA, -k_vx(kk(:)')/dx^2];
        %         %above current
        %         jj = ind(1:nz1-1,2:nx1); kk = ind(1:nz1-1,2:nx1);
        %         II = [II, ii(:)'];   JJ = [JJ, jj(:)'];   AA = [AA, -k_vz(kk(:)')/dz^2];
        %         %below current
        %         jj = ind(3:nz1+1,2:nx1); kk = ind(2:nz1,2:nx1);
        %         II = [II, ii(:)'];   JJ = [JJ, jj(:)'];   AA = [AA, -k_vz(kk(:)')/dz^2];
        %         %current node
        %         RhoCp_mid = Rho_mid(2:nz1,2:nx1).*Cp_mid(2:nz1,2:nx1);
        %         Asum = (RhoCp_mid./dt)...
        %             +((k_vx(2:nz1,2:nx1) + k_vx(2:nz1,1:nx1-1))/dx^2)...
        %             +((k_vz(2:nz1,2:nx1) + k_vz(1:nz1-1,2:nx1))/dz^2);
        %         jj = ii;
        %         II = [II, ii(:)'];   JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'];
        %         %RHS
        %         RRSUM = (Hs(2:nz1,2:nx1)+Ha(2:nz1,2:nx1)+Hr(2:nz1,2:nx1))+RhoCp_mid.*T_mid(2:nz1,2:nx1)./dt;
        %         IR = [IR, ii(:)']; RR = [RR, RRSUM(:)'];
        %
        %         % boundary points
        %         %top constant temp
        %         ii  = ind(1,2:nx1); jj  = ind(1,2:nx1); jj1 = ind(2,2:nx1); %below top boundary
        %         Asum = zeros(size(ii));
        %         II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
        %         II = [II, ii(:)']; JJ = [JJ, jj1(:)'];  AA = [AA, Asum(:)'+1];
        %         %RHS
        %         IR = [IR, ii(:)']; RR = [RR, Asum(:)'*T_top*2];
        %
        %         %bottom constant temp
        %         ii  = ind(nz2,2:nx1); jj  = ind(nz2,2:nx1); jj1 = ind(nz1,2:nx1); %below top boundary
        %         Asum = zeros(size(ii)); Asum1 = Asum.*T_mid(nz2,2:nx1).*2;
        %         II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
        %         II = [II, ii(:)']; JJ = [JJ, jj1(:)'];  AA = [AA, Asum(:)'+1];
        %         %RHS
        %         IR = [IR, ii(:)']; RR = [RR, Asum1(:)'];
        %
        %         %left insulating boundary
        %         ii  = ind(1:nz2,1); jj  = ind(1:nz2,1); jj1 = ind(1:nz2,2); %below top boundary
        %         Asum = zeros(size(ii));
        %         II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
        %         II = [II, ii(:)']; JJ = [JJ, jj1(:)'];  AA = [AA, Asum(:)'-1];
        %         %RHS
        %         IR = [IR, ii(:)']; RR = [RR, Asum(:)'];
        %
        %         %right insulating boundary
        %         ii  = ind(1:nz2,nx2); jj  = ind(1:nz2,nx2); jj1 = ind(1:nz2,nx1); %below top boundary
        %         Asum = zeros(size(ii));
        %         II = [II, ii(:)']; JJ = [JJ, jj(:)'];   AA = [AA, Asum(:)'+1];
        %         II = [II, ii(:)']; JJ = [JJ, jj1(:)'];  AA = [AA, Asum(:)'-1];
        %         %RHS
        %         IR = [IR, ii(:)']; RR = [RR, Asum(:)'];
        %
        %         A       = sparse(II,JJ,AA,NP,NP);
        %         RHS     = sparse(IR,ones(size(IR)),RR,NP,1);
        %
        %         X           =  sqrt(abs(diag(A)));
        %         X           =  diag(sparse(1./X));
        %
        %         A           =  X*A*X;
        %         RHS         =  X*RHS;
        %         %% output T grid
        %         Tvec = X*(A\RHS); % backslash operator to output diffused temperature vector
        %         T_out = reshape(Tvec,size(T_out));
        %         % compute dT
        %         dT          = T_out-T_mid;
        %         % thermal timestepping conditions
        %         maxdT       = max(max(abs(dT)));
end




