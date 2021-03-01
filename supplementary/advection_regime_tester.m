% test and compare advection schemes

% prepare workspace
clear all
close all

% set grid parameter
L   = 1e3;
N   = 400;
h   = L/N;

% set run options
RunID       = 'central gaussian rotation';
movement    = 'rotation';
a0          = 1;     % initialbackground value of advected quantity
a1          = 1;     % initial amplitude of advected quantity
radius      = L/20;  % initial radius of gaussian in advected quantity
CFL         = 0.5;   % Courant number
nt          = 10;

% get coordinate grids
% interior cell nodes
x   = -L/2+h/2:h:(L/2 - h/2);
z   = -L/2+h/2:h:(L/2 - h/2);
[x2d,z2d] = meshgrid(x,z);

% ghosted cell nodes
xa  = -L/2-h/2:h:L/2+h/2;
za  = -L/2-h/2:h:L/2+h/2;
[xa2d,za2d] = meshgrid(xa,za);

% staggered x-face nodes
xu  = -L/2:h:L/2;
zu  = -L/2-h/2:h:L/2+h/2;
[xu2d,zu2d] = meshgrid(xu,zu);

% staggered z-face nodes
xw  = -L/2-h/2:h:L/2+h/2;
zw  = -L/2:h:L/2;
[xw2d,zw2d] = meshgrid(xw,zw);

% initialise arrays for velocity
u   = zeros(N+2,N+1);  % staggered vx
w   = zeros(N+1,N+2);  % staggered vz

% set initial condition for advected quantity
x0  = L/6; z0 = L/6;
gsn = exp(- (xa2d-x0).^2./radius.^2 - (za2d-z0).^2./radius.^2 );
a   = a0 + a1.*gsn; 

% plot initial condition
% figure(1)
% subplot(2,2,1)
% imagesc(x,z,a(2:end-1,2:end-1)); 
% colorbar
% title('initial')

% set velocity field
switch movement
    case 'rotation'
        u = -zu2d;
        w =  xw2d;  
        umid = -za2d;
        wmid =  xa2d;
end

afr  = a;  dadtfr = 0.*a(2:end-1,2:end-1);
au1  = a;  dadtu1 = 0.*a(2:end-1,2:end-1);
au2  = a;  dadtu2 = 0.*a(2:end-1,2:end-1);
au3  = a;  dadtu3 = 0.*a(2:end-1,2:end-1);

time = 0;
dt   = CFL * min( h/2/max(abs(u(:))) , h/2/max(abs(w(:))) );
ti   = 0;

while max([max(afr(:))-min(afr(:)),max(au2(:))-min(au2(:)),max(au3(:))-min(au3(:))]) < 10*a1 && ti <= nt

    % store previous advection rates
    dadtfro = dadtfr;
    dadtu1o = dadtu1;
    dadtu2o = dadtu2;
    dadtu3o = dadtu3;

    dadtfr  = advection(u,w,afr,h,'fromm');
    dadtu1  = advection(u,w,au1,h,'first upwind');
    dadtu2  = advection(u,w,au2,h,'second upwind');
    dadtu3  = advection(u,w,au3,h,'third upwind');
    
    % advance advection time step
    afr(2:end-1,2:end-1)    = afr(2:end-1,2:end-1)  - (dadtfr+dadtfro)/2*dt;
    au1(2:end-1,2:end-1)    = au1(2:end-1,2:end-1)  - (dadtu1+dadtu1o)/2*dt;
    au2(2:end-1,2:end-1)    = au2(2:end-1,2:end-1)  - (dadtu2+dadtu2o)/2*dt;
    au3(2:end-1,2:end-1)    = au3(2:end-1,2:end-1)  - (dadtu3+dadtu3o)/2*dt;

    % advance reference gaussian
    x0   = x0 + umid*dt; z0 = z0 + wmid*dt;
    gsn  = exp(- (xa2d-x0).^2./radius.^2 - (za2d-z0).^2./radius.^2 );
    res_afr = afr - (a0 + a1.*gsn); 
    res_au1 = au1 - (a0 + a1.*gsn); 
    res_au2 = au2 - (a0 + a1.*gsn); 
    res_au3 = au3 - (a0 + a1.*gsn); 

    resnorm_afr = norm(res_afr(:),2)./norm(afr,2);
    resnorm_au1 = norm(res_au1(:),2)./norm(afr,2);
    resnorm_au2 = norm(res_au2(:),2)./norm(afr,2);
    resnorm_au3 = norm(res_au3(:),2)./norm(afr,2);
    
    fprintf(1,'   --- res fromm = %1.3e; res upw1 = %1.3e; res upw2 = %1.3e; res upw3 = %1.3e;\n',resnorm_afr,resnorm_au1,resnorm_au2,resnorm_au3);
        
    % plot results
    if ~mod(ti,10)  
        
        figure(1)
        subplot(2,2,1)
        imagesc(x,z,afr);
        colorbar
        title('fromm scheme')
        
        subplot(2,2,2)
        imagesc(x,z,au1);
        colorbar
        title('1st order upwind')
        
        subplot(2,2,3)
        imagesc(x,z,au2);
        colorbar
        title('2nd order upwind')
        
        subplot(2,2,4)
        imagesc(x,z,au3);
        colorbar
        title('3rd order upwind')
        drawnow;
        
        figure(2)
        subplot(2,2,1)
        imagesc(x,z,res_afr);
        colorbar
        title('fromm scheme')
        
        subplot(2,2,2)
        imagesc(x,z,res_au1);
        colorbar
        title('1st order upwind')
        
        subplot(2,2,3)
        imagesc(x,z,res_au2);
        colorbar
        title('2nd order upwind')
        
        subplot(2,2,4)
        imagesc(x,z,res_au3);
        colorbar
        title('3rd order upwind')
        drawnow;
        
%     figure(2)
%     subplot(2,2,1)
%     quiver(x(1:10:end),z(1:10:end),umid(2:10:end-1,1:10:end),wmid(1:10:end,2:10:end-1))
%     axis ij image
%     
%     subplot(2,2,2)
%     imagesc(x,z,dadtf(2:end-1,2:end-1))
%     colorbar
%     title('fromm scheme da/dt')
%     
%     subplot(2,2,3)
%     imagesc(x,z,dadtu3(2:end-1,2:end-1))
%     colorbar
%     title('3rd order upwind da/dt')
%     
%     subplot(2,2,4)
%     imagesc(x,z,dadtu2(2:end-1,2:end-1))
%     colorbar
%     title('2nd order upwind da/dt')
    
%     pause(0.1)
%     saveas(figure(1),[pwd '/out/', RunID, ' ', num2str(ti), 'iterations.jpg'])
%     saveas(figure(2),[pwd '/out/dadt ', RunID, ' ', num2str(ti), 'iterations.jpg'])
    end
    % update time step
    ti = ti+1;
    time = time +dt;
end

function adv = advection(u,w,a,h,output)

wp  = w(2:end  ,2:end-1);
wm  = w(1:end-1,2:end-1);
up  = u(2:end-1,2:end  );
um  = u(2:end-1,1:end-1);

vx  = (up+um)./2;
vz  = (wp+wm)./2;

vxp = max(vx,0); vxm = min(vx,0);
vzp = max(vz,0); vzm = min(vz,0);
        
agh                    = zeros(size(a)+2);
agh(2:end-1,2:end-1)   = a;

agh([1 2 end-1 end],:) = agh([4 3 end-2 end-3],:);
agh(:,[1 2 end-1 end]) = agh(:,[4 3 end-2 end-3]);

acc = agh(3:end-2,3:end-2);
ajp = agh(4:end-1,3:end-2);  ajpp = agh(5:end-0,3:end-2);
ajm = agh(2:end-3,3:end-2);  ajmm = agh(1:end-4,3:end-2);
aip = agh(3:end-2,4:end-1);  aipp = agh(3:end-2,5:end-0);
aim = agh(3:end-2,2:end-3);  aimm = agh(3:end-2,1:end-4);

switch output
    case 'fromm'

        adv   =     up .*(-(aipp-aip)./h./8 + (aip + acc)./h./2 + (acc-aim )./h./8) ...
              - abs(up).*(-(aipp-aip)./h./8 + (aip - acc)./h./4 - (acc-aim )./h./8) ...
              -     um .*(-(aip -acc)./h./8 + (acc + aim)./h./2 + (aim-aimm)./h./8) ...
              + abs(um).*(-(aip -acc)./h./8 + (acc - aim)./h./4 - (aim-aimm)./h./8) ...
              +     wp .*(-(ajpp-ajp)./h./8 + (ajp + acc)./h./2 + (acc-ajm )./h./8) ...
              - abs(wp).*(-(ajpp-ajp)./h./8 + (ajp - acc)./h./4 - (acc-ajm )./h./8) ...
              -     wm .*(-(ajp -acc)./h./8 + (acc + ajm)./h./2 + (ajm-ajmm)./h./8) ...
              + abs(wm).*(-(ajp -acc)./h./8 + (acc - ajm)./h./4 - (ajm-ajmm)./h./8);
          
    case 'first upwind'
        
        axp   = (aip-acc)./h;
        axm   = (acc-aim)./h;
        azp   = (ajp-acc)./h;
        azm   = (acc-ajm)./h;
        
        daxdt = vxp.*axm + vxm.*axp;
        dazdt = vzp.*azm + vzm.*azp;
        
        adv   = daxdt + dazdt;
        
    case 'second upwind'

        axp   = (-3*acc+4*aip-aipp)/2/h;
        axm   = ( 3*acc-4*aim+aimm)/2/h;
        azp   = (-3*acc+4*ajp-ajpp)/2/h;
        azm   = ( 3*acc-4*ajm+ajmm)/2/h;
        
        daxdt = vxp.*axm + vxm.*axp;
        dazdt = vzp.*azm + vzm.*azp;
        
        adv   = daxdt + dazdt;        
        
    case 'third upwind'

        axp   = (-2*aim-3*acc+6*aip-aipp)/6/h;
        axm   = ( 2*aip+3*acc-6*aim+aimm)/6/h;
        azp   = (-2*ajm-3*acc+6*ajp-ajpp)/6/h;
        azm   = ( 2*ajp+3*acc-6*ajm+ajmm)/6/h;
        
        daxdt = vxp.*axm + vxm.*axp;
        dazdt = vzp.*azm + vzm.*azp;
        
        adv   = daxdt + dazdt;
            
end  
end
