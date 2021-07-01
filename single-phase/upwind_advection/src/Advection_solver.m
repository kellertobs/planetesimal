% advection solvers, currently only fromm is vectorised
function[dTdt,dMdt] = Advection_solver(vx_out,vz_out,vx_mid,vz_mid,T,Material,dz,dx,nx,nz,nx1,nz1,AdvRegime)
w = vz_out(1:nz1,:); % vz
u = vx_out(:,1:nx1); % vx

u1 = vx_mid(2:nz1,2:nx1); w1 = vz_mid(2:nz1,2:nx1);

wp = w(2:end,2:end-1);
wm = w(1:end-1,2:end-1);
up = u(2:end-1,2:end);
um = u(2:end-1,1:end-1);

a = T; %advected quantity a (Temperature
b = Material;

agh                    = zeros(size(a)+2);
agh(2:end-1,2:end-1)   = a;
bgh                    = zeros(size(b)+2);
bgh(2:end-1,2:end-1)   = b;


% agh([1 2 end-1 end],:) = agh([3 3 end-2 end-2],:);
% agh(:,[1 2 end-1 end]) = agh(:,[3 3 end-2 end-2]);
agh([1 2 end-1 end],:) = agh([4 3 end-2 end-3],:);
agh(:,[1 2 end-1 end]) = agh(:,[4 3 end-2 end-3]);


% bgh([1 2 end-1 end],:) = bgh([3 3 end-2 end-2],:);
% bgh(:,[1 2 end-1 end]) = bgh(:,[3 3 end-2 end-2]);
bgh([1 2 end-1 end],:) = bgh([4 3 end-2 end-3],:);
bgh(:,[1 2 end-1 end]) = bgh(:,[4 3 end-2 end-3]);

acc = agh(3:end-2,3:end-2);
ajp = agh(4:end-1,3:end-2);  ajpp = agh(5:end-0,3:end-2);
ajm = agh(2:end-3,3:end-2);  ajmm = agh(1:end-4,3:end-2);
aip = agh(3:end-2,4:end-1);  aipp = agh(3:end-2,5:end-0);
aim = agh(3:end-2,2:end-3);  aimm = agh(3:end-2,1:end-4);

bcc = bgh(3:end-2,3:end-2);
bjp = bgh(4:end-1,3:end-2);  bjpp = bgh(5:end-0,3:end-2);
bjm = bgh(2:end-3,3:end-2);  bjmm = bgh(1:end-4,3:end-2);
bip = bgh(3:end-2,4:end-1);  bipp = bgh(3:end-2,5:end-0);
bim = bgh(3:end-2,2:end-3);  bimm = bgh(3:end-2,1:end-4);


switch AdvRegime
    case 'fromm'

        dTdt  =     up .*(-(aipp-aip)./dx./8 + (aip + acc)./dx./2 + (acc-aim )./dx./8) ...
              - abs(up).*(-(aipp-aip)./dx./8 + (aip - acc)./dx./4 - (acc-aim )./dx./8) ...
              -     um .*(-(aip -acc)./dx./8 + (acc + aim)./dx./2 + (aim-aimm)./dx./8) ...
              + abs(um).*(-(aip -acc)./dx./8 + (acc - aim)./dx./4 - (aim-aimm)./dx./8) ...
              +     wp .*(-(ajpp-ajp)./dz./8 + (ajp + acc)./dz./2 + (acc-ajm )./dz./8) ...
              - abs(wp).*(-(ajpp-ajp)./dz./8 + (ajp - acc)./dz./4 - (acc-ajm )./dz./8) ...
              -     wm .*(-(ajp -acc)./dz./8 + (acc + ajm)./dz./2 + (ajm-ajmm)./dz./8) ...
              + abs(wm).*(-(ajp -acc)./dz./8 + (acc - ajm)./dz./4 - (ajm-ajmm)./dz./8);
          
        dMdt  =     up .*(-(bipp-bip)./dx./8 + (bip + bcc)./dx./2 + (bcc-bim )./dx./8) ...
              - abs(up).*(-(bipp-bip)./dx./8 + (bip - bcc)./dx./4 - (bcc-bim )./dx./8) ...
              -     um .*(-(bip -bcc)./dx./8 + (bcc + bim)./dx./2 + (bim-bimm)./dx./8) ...
              + abs(um).*(-(bip -bcc)./dx./8 + (bcc - bim)./dx./4 - (bim-bimm)./dx./8) ...
              +     wp .*(-(bjpp-bjp)./dz./8 + (bjp + bcc)./dz./2 + (bcc-bjm )./dz./8) ...
              - abs(wp).*(-(bjpp-bjp)./dz./8 + (bjp - bcc)./dz./4 - (bcc-bjm )./dz./8) ...
              -     wm .*(-(bjp -bcc)./dz./8 + (bcc + bjm)./dz./2 + (bjm-bjmm)./dz./8) ...
              + abs(wm).*(-(bjp -bcc)./dz./8 + (bcc - bjm)./dz./4 - (bjm-bjmm)./dz./8);
    case 'first upwind'              
        vxp = max(u1,0); vxm = min(u1,0);
        vzp = max(w1,0); vzm = min(w1,0);
        axp = (agh(3:end-2,4:end-1)-agh(3:end-2,3:end-2))./dx;
        axm = (agh(3:end-2,3:end-2)-agh(3:end-2,2:end-3))./dx;
        azp = (agh(4:end-1,3:end-2)-agh(3:end-2,3:end-2))./dz;
        azm = (agh(3:end-2,3:end-2)-agh(2:end-3,3:end-2))./dz;
        daxdt = vxp.*axm + vxm.*axp;
        dazdt = vzp.*azm + vzm.*azp;
        dTdt = daxdt + dazdt;
       
%         vxp = max(u1,0); vxm = min(u1,0);
%         vzp = max(w1,0); vzm = min(w2,0);
        bxp = (bgh(3:end-2,4:end-1)-bgh(3:end-2,3:end-2))./dx;
        bxm = (bgh(3:end-2,3:end-2)-bgh(3:end-2,2:end-3))./dx;
        bzp = (bgh(4:end-1,3:end-2)-bgh(3:end-2,3:end-2))./dz;
        bzm = (bgh(3:end-2,3:end-2)-bgh(2:end-3,3:end-2))./dz;
        dbxdt = vxp.*bxm + vxm.*bxp;
        dbzdt = vzp.*bzm + vzm.*bzp;
        dMdt = dbxdt + dbzdt;
        
    case 'second upwind'
        vxp = max(u1,0); vxm = min(u1,0);
        vzp = max(w1,0); vzm = min(w1,0);
        axp = (-3*agh(3:end-2,3:end-2)+4*agh(3:end-2,4:end-1)-agh(3:end-2,5:end))/2/dx;
        axm = (3*agh(3:end-2,3:end-2)-4*agh(3:end-2,2:end-3)+agh(3:end-2,1:end-4))/2/dx;
        azp = (-3*agh(3:end-2,3:end-2)+4*agh(4:end-1,3:end-2)-agh(5:end,3:end-2))/2/dz;
        azm = (3*agh(3:end-2,3:end-2)-4*agh(2:end-3,3:end-2)+agh(1:end-4,3:end-2))/2/dz;
        daxdt = vxp.*axm + vxm.*axp;
        dazdt = vzp.*azm + vzm.*azp;
        dTdt = daxdt + dazdt;   
        
        bxp = (-3*bgh(3:end-2,3:end-2)+4*bgh(3:end-2,4:end-1)-bgh(3:end-2,5:end))/2/dx;
        bxm = ( 3*bgh(3:end-2,3:end-2)-4*bgh(3:end-2,2:end-3)+bgh(3:end-2,1:end-4))/2/dx;
        bzp = (-3*bgh(3:end-2,3:end-2)+4*bgh(4:end-1,3:end-2)-bgh(5:end,3:end-2))/2/dz;
        bzm = ( 3*bgh(3:end-2,3:end-2)-4*bgh(2:end-3,3:end-2)+bgh(1:end-4,3:end-2))/2/dz;
        dbxdt = vxp.*bxm + vxm.*bxp;
        dbzdt = vzp.*bzm + vzm.*bzp;
        dMdt = dbxdt + dbzdt; 
        
    case 'third upwind'
        vxp = max(u1,0); vxm = min(u1,0);
        vzp = max(w1,0); vzm = min(w1,0);
        axp = (-2*agh(3:end-2,2:end-3)-3*agh(3:end-2,3:end-2)+6*agh(3:end-2,4:end-1)-agh(3:end-2,5:end))/6  /dx;
        axm = ( 2*agh(3:end-2,4:end-1)+3*agh(3:end-2,3:end-2)-6*agh(3:end-2,2:end-3)+agh(3:end-2,1:end-4))/6/dx;
        azp = (-2*agh(2:end-3,3:end-2)-3*agh(3:end-2,3:end-2)+6*agh(4:end-1,3:end-2)-agh(5:end,3:end-2))  /6/dz;
        azm = ( 2*agh(4:end-1,3:end-2)+3*agh(3:end-2,3:end-2)-6*agh(2:end-3,3:end-2)+agh(1:end-4,3:end-2))/6/dz;
        daxdt = vxp.*axm + vxm.*axp;
        dazdt = vzp.*azm + vzm.*azp;
        dTdt = daxdt + dazdt;
        
        bxp = (-2*bgh(3:end-2,2:end-3)-3*bgh(3:end-2,3:end-2)+6*bgh(3:end-2,4:end-1)-bgh(3:end-2,5:end))/6  /dx;
        bxm = ( 2*bgh(3:end-2,4:end-1)+3*bgh(3:end-2,3:end-2)-6*bgh(3:end-2,2:end-3)+bgh(3:end-2,1:end-4))/6/dx;
        bzp = (-2*bgh(2:end-3,3:end-2)-3*bgh(3:end-2,3:end-2)+6*bgh(4:end-1,3:end-2)-bgh(5:end,3:end-2))  /6/dz;
        bzm = ( 2*bgh(4:end-1,3:end-2)+3*bgh(3:end-2,3:end-2)-6*bgh(2:end-3,3:end-2)+bgh(1:end-4,3:end-2))/6/dz;
        dbxdt = vxp.*bxm + vxm.*bxp;
        dbzdt = vzp.*bzm + vzm.*bzp;
        dMdt = dbxdt + dbzdt;
end
  

  
