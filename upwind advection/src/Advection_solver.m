% advection solvers, currently only fromm is vectorised
function[dTdt,dMdt] = Advection_solver(vx_out,vz_out,vx_mid,vz_mid,T,Material,dz,dx,nx,nz,nx1,nz1,output)
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


switch output
    case 'fromm'

        dTdt  =     up .*(-(aipp-aip)./dx./8 + (aip + acc)./dx./2 + (acc-aim )./dx./8) ...
              - abs(up).*(-(aipp-aip)./dx./8 + (aip - acc)./dx./4 - (acc-aim )./dx./8) ...
              -     um .*(-(aip -acc)./dx./8 + (acc + aim)./dx./2 + (aim-aimm)./dx./8) ...
              + abs(um).*(-(aip -acc)./dx./8 + (acc - aim)./dx./4 - (aim-aimm)./dx./8) ...
              +     wp .*(-(ajpp-ajp)./dx./8 + (ajp + acc)./dx./2 + (acc-ajm )./dx./8) ...
              - abs(wp).*(-(ajpp-ajp)./dx./8 + (ajp - acc)./dx./4 - (acc-ajm )./dx./8) ...
              -     wm .*(-(ajp -acc)./dx./8 + (acc + ajm)./dx./2 + (ajm-ajmm)./dx./8) ...
              + abs(wm).*(-(ajp -acc)./dx./8 + (acc - ajm)./dx./4 - (ajm-ajmm)./dx./8);
          
        dMdt  =     up .*(-(bipp-bip)./dx./8 + (bip + bcc)./dx./2 + (bcc-bim )./dx./8) ...
              - abs(up).*(-(bipp-bip)./dx./8 + (bip - bcc)./dx./4 - (bcc-bim )./dx./8) ...
              -     um .*(-(bip -bcc)./dx./8 + (bcc + bim)./dx./2 + (bim-bimm)./dx./8) ...
              + abs(um).*(-(bip -bcc)./dx./8 + (bcc - bim)./dx./4 - (bim-bimm)./dx./8) ...
              +     wp .*(-(bjpp-bjp)./dz./8 + (bjp + bcc)./dz./2 + (bcc-bjm )./dz./8) ...
              - abs(wp).*(-(bjpp-bjp)./dz./8 + (bjp - bcc)./dz./4 - (bcc-bjm )./dz./8) ...
              -     wm .*(-(bjp -bcc)./dz./8 + (bcc + bjm)./dz./2 + (bjm-bjmm)./dz./8) ...
              + abs(wm).*(-(bjp -bcc)./dz./8 + (bcc - bjm)./dz./4 - (bjm-bjmm)./dz./8);
    case 'first upwind'
        adv = zeros(size(a)-2);        
        vxp = max(u1,0); vxm = min(u1,0);
        vzp = max(w1,0); vzm = min(w2,0);
        axp = (agh(3:end-2,4:end-1)-agh(3:end-2,3:end-2))./h;
        axm = (agh(3:end-2,3:end-2)-agh(3:end-2,2:end-3))./h;
        azp = (agh(4:end-1,3:end-2)-agh(3:end-2,3:end-2))./h;
        azm = (agh(3:end-2,3:end-2)-agh(2:end-3,3:end-2))./h;
        daxdt = vxp.*axm + vxm.*axp;
        dazdt = vzp.*azm + vzm.*azp;
        adv = daxdt + dazdt;
        
    case 'second upwind'
        adv = zeros(size(a)-2);
        for i = 1:1:nz
        for j = 1:1:nx
            vx = (up(i,j)+um(i,j))/2;
            vz = (wp(i,j)+wm(i,j))/2;

            if vx>0
                dadx = vx*(3*agh(i+2,j+2)-4*agh(i+2,j+1)+agh(i+2,j))/2/dx;
            else
                dadx = vx*(-3*agh(i+2,j+2)+4*agh(i+2,j+3)-agh(i+2,j+4))/2/dx;
            end
            if vz>0
                dadz = vz*(3*agh(i+2,j+2)-4*agh(i+1,j+2)+agh(i,j+2))/2/dz;
            else
                dadz = vz*(-3*agh(i+2,j+2)+4*agh(i+3,j+2)-agh(i+4,j+2))/2/dz;
            end
            adv(i,j) = dadx+dadz;
        end
        end
        
    case 'third upwind'
        adv = zeros(size(a)-2);
        for i = 1:1:nz
        for j = 1:1:nx
            vx = (up(i,j)+um(i,j))/2;
            vz = (wp(i,j)+wm(i,j))/2;

            if vx>0
                dadx = vx*(2*agh(i+2,j+3)+3*agh(i+2,j+2)-6*agh(i+2,j+1)+agh(i+2,j))/6/dx;
            else
                dadx = vx*(-2*agh(i+2,j+1)-3*agh(i+2,j+2)+6*agh(i+2,j+3)-agh(i+2,j+4))/6/dx;
            end
            if vz>0
                dadz = vz*(2*agh(i+3,j+2)+3*agh(i+2,j+2)-6*agh(i+1,j+2)+agh(i,j+2))/6/dz;
            else
                dadz = vz*(-2*agh(i+1,j+2)-3*agh(i+2,j+2)+6*agh(i+3,j+2)-agh(i+4,j+2))/6/dz;
            end
            adv(i,j) = dadx+dadz;
        end
        end        
end
  

  
