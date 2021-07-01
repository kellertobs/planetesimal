% planetesimal: solve thermo-chemical equations

% print solver header
fprintf(1,'  ---  solve thermo-chemical equations \n');
tic;


%% Solve energy equation explicitly
To    = SOL.T;     % store previous temperature solution
dTdto = SOL.dTdt;  % store previous rate of change

% get advection rate
advn_T = advection(SOL.U,SOL.W,SOL.T,NUM.dx,NUM.dz,NUM.AdvnScheme);

% get diffusion rate
diff_T = (diff(SOL.T(:,2:end-1),2,1)./NUM.dz^2  ...
       +  diff(SOL.T(2:end-1,:),2,2)./NUM.dx^2) ...
       .* MAT.kT(2:end-1,2:end-1);

% heat capacity density
RhoCp = MAT.Rho.*MAT.Cp;

% get total rate of change
SOL.dTdt = - advn_T + (diff_T + MAT.Hr(2:end-1,2:end-1) + SOL.Hs(2:end-1,2:end-1) + SOL.Ha(2:end-1,2:end-1))./RhoCp(2:end-1,2:end-1);

% update temperature solution
SOL.T(2:end-1,2:end-1) = SOL.T(2:end-1,2:end-1) + (  NUM.theta .*SOL.dTdt   ...
                                                + (1-NUM.theta).*    dTdto) .* NUM.dt;

% apply top boundary conditions
switch SOL.BCTempTop
    case 'isothermal'
        SOL.T(1,:) =     To(1,:);
    case 'insulating'
        SOL.T(1,:) = SOL.T (2,:);
end

% apply bottom boundary conditions
switch SOL.BCTempBot
    case 'isothermal'
        SOL.T(end,:) =     To(end  ,:);
    case 'insulating'
        SOL.T(end,:) = SOL.T (end-1,:);
end

% apply side boundary conditions
switch SOL.BCTempSides
    case 'isothermal'
        SOL.T(:,[1 end]) =     To(:,[1 end  ]);
    case 'insulating'
        SOL.T(:,[1 end]) = SOL.T (:,[2 end-1]);
end

toc_solve = toc;  % stop clock temperature solution
fprintf(1,'       solution time %1.4f s \n\n',toc_solve);