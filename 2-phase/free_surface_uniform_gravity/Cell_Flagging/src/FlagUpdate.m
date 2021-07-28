% setup free surface
NUM.Fs = TopoFun( NUM.xP, NUM.nzP, NUM.nxP, NUM.D/10, 'horizontal');

%% P-cells 
NUM.flag.P = zeros(NUM.nzP,NUM.nxP);
% deactivate cells above the free surface
NUM.flag.P(NUM.ZP<NUM.D/10) = 3;
% activate Free surface cells (need to figure out a way to vextorize)
% tic
for j = 1:1:NUM.nxP
    
    for i = 1:1:NUM.nzP
        if (NUM.flag.P(i,j) == 3 && NUM.flag.P(i+1,j) == 0)
        NUM.flag.P(i,j) = 6;
        end
    end
end
% toc

%% Vx-nodes; 0 = active 
NUM.flag.U = ones(NUM.nzU,NUM.nxU);
NUM.flag.U = NUM.flag.P(:,1:end-1).*NUM.flag.P(:,2:end);

%% Vz-nodes; 0 = active
NUM.flag.W = zeros(NUM.nzW,NUM.nxW);
NUM.flag.W = NUM.flag.P(1:end-1,:).*NUM.flag.P(2:end,:);

%% corner nodes
NUM.flag.C = zeros(NUM.nzC,NUM.nxC);
NUM.flag.C = (NUM.flag.P(1:end-1,1:end-1) ...
    +  NUM.flag.P(2:end  ,1:end-1) ...
    +  NUM.flag.P(1:end-1,2:end  ) ...
    +  NUM.flag.P(2:end  ,2:end  ));
NUM.flag.C(NUM.flag.C<=6) = 0;

% figure()
% imagesc(NUM.xP,NUM.zP,NUM.flag.P)
% colorbar
% hold on
% plot(NUM.xP,NUM.Fs,'--k', 'linewidth', 3)
% figure()
% imagesc(NUM.xU,NUM.zU,NUM.flag.U)
% colorbar
% hold on
% plot(NUM.xP,NUM.Fs,'--k', 'linewidth', 3)
% figure()
% imagesc(NUM.xW,NUM.zW,NUM.flag.W)
% colorbar
% hold on
% plot(NUM.xP,NUM.Fs,'--k', 'linewidth', 3)
% figure()
% imagesc(NUM.xC,NUM.zC,NUM.flag.C)
% colorbar
% hold on
% plot(NUM.xP,NUM.Fs,'--k', 'linewidth', 3)

% Topographic function
function Y = TopoFun( X, z, x, Ds, TopoType)
switch TopoType
    case 'horizontal'
        Y = zeros(z, x) + Ds;
    case 'depression'
        
    case 'bulge'
        
    case 'square'
        
    case 'circle'
        
end
end