%numerical setup file
%preferably not user editable

%% setup numerical grid
% setup vectors
[x,z,xvx,zvx,xvz,zvz,xp,zp] = vectorsetup(D,L,dx,dz);

% setup meshes

%% setup grid

%setup temperature grid
T0      = T_mantle; %(placeholder)
T       = zeros(nz1,nx1) + T0;
for j = 1:1:nx1
for i = 1:1:nz1
    w = L/5;
    T(m) = T(m) + 300*exp(- (xp(m)-L/2)^2/w^2 - (zp(m)-D/2)^2/w^2 );
end
end
% %loop to set markers
% for jm=1:1:nxm_all
% for im=1:1:nzm_all
%     % Define marker coordinates
%     xm(m)   = dxm/2+(jm-1)*dxm+(rand-0.5)*dxm; % randomly distributed
%     zm(m)   = dzm/2+(im-1)*dzm+(rand-0.5)*dzm;
%     % Input marker properties
%     switch Ambtype
%         case 'linear'           
%             Tm(m)       = T_top + abs(zm(m))/D*(T_bot-T_top); 
%         case 'constant'
%             Tm(m)       = T_top;
%         case 'gaussian'
%             Tm(m)       = T_top;
%             
%     end
% %     rmark=((xm(m)-L/2)^2+(zm(m)-D/2)^2)^0.5;
% %     if rmark<(L/5)
% %         Tm(m) = 2000;
% %     end
%     Etam(m)     = Eta_mantle;    % Viscosity
%     Alpham(m)   = Alpha_mantle;  % Thermal expansion
%     Cpm(m)      = Cp_mantle;  % Heat capacity
%     Kappam(m)   = Kappa_mantle;  % Thermal conductivity
%     Hrm(m)      = Hr_mantle;     % Radiogenic heating
%     if zm(m)>D*0.95
%         Mtype(m) = 2;
%     elseif zm(m)>D*0.9
%         Mtype(m) = 1;
%     elseif zm(m)>D*0.8
%         Mtype(m) = 2;
%     elseif zm(m)>D*0.7
%         Mtype(m) = 1;
%     elseif zm(m)>D*0.6
%         Mtype(m) = 2;
%     elseif zm(m)>D*0.5
%         Mtype(m) = 1;
%     elseif zm(m)>D*0.3
%         Mtype(m) = 2;
%     elseif zm(m)>D*0.1
%         Mtype(m) = 1;
%     else
%         Mtype(m) = 2;
%     end
%         
% 
%     % Update marker counter
%     m=m+1;
% end
% end

%% initialise arrays
Epsxz           = zeros(Nz,Nx);     % strain rate on the ordinary grid
Sigxz           = zeros(Nz,Nx);     % deviatoric stress on the ordinary grid
Epsxx           = Epsxz;            % strain rate in 
Sigxx           = Sigxz;            % deviatoric stress the middle of grid/pressure nodes
Hs              = Sigxx;            % shear heating, on the pressure nodes
Ha              = Hs;               % adiabatic heating, on pressure nodes
Hr              = Hs;               % radiogenic heating


%% create indexing system
Number  = numsetup(Nz,Nx);  % ordinary grid
