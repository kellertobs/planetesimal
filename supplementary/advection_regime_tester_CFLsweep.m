% prepare workspace
clear all
close all
% profile off
% profile on
% set run options
% RunID       = 'gaussian_rotation';
movement    = 'rotation';
L           = 1;     % domain length
N           = 200;   % number of grid points
h           = L/N;   % grid spacing
CFL         = 1;  % Courant number
theta       = 0.5;  % time stepping mode (0.5 = Crank-Nicolson, 0 = backward Euler, 1 = forward Euler)
M           = 1;   % number of time steps to take
a0          = 1;     % initialbackground value of advected quantity
a1          = 1;     % initial amplitude of advected quantity
radius      = L/20;  % initial radius of gaussian in advected quantity

% RunID = ['CFL' num2str(CFL) ' dx' num2str(h)]

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

% if ~exist([pwd '/out/',RunID],'dir'); mkdir([pwd '/out/',RunID]); end

% set velocity field
switch movement
    case 'rotation'
        u = -zu2d;
        w =  xw2d;
        umid = -za2d;
        wmid =  xa2d;
end

acd  = a;  dadtcd = 0.*a(2:end-1,2:end-1);
afd  = a;  dadtfd = 0.*a(2:end-1,2:end-1);
afr  = a;  dadtfr = 0.*a(2:end-1,2:end-1);
au1  = a;  dadtu1 = 0.*a(2:end-1,2:end-1);
au2  = a;  dadtu2 = 0.*a(2:end-1,2:end-1);
au3  = a;  dadtu3 = 0.*a(2:end-1,2:end-1);
CFL_all = linspace(0.1,1,10);
resend_acd=[];
resend_afd=[];
resend_afr=[];
resend_au1=[];
resend_au2=[];
resend_au3=[];
for i = 1:1:length(CFL_all)
    CFL = CFL_all(i);
    M = 3;
    x0  = L/6; z0 = L/6;
    acd  = a;  dadtcd = 0.*a(2:end-1,2:end-1);
    afd  = a;  dadtfd = 0.*a(2:end-1,2:end-1);
    afr  = a;  dadtfr = 0.*a(2:end-1,2:end-1);
    au1  = a;  dadtu1 = 0.*a(2:end-1,2:end-1);
    au2  = a;  dadtu2 = 0.*a(2:end-1,2:end-1);
    au3  = a;  dadtu3 = 0.*a(2:end-1,2:end-1);
    fprintf(1,'\n\n*****  Advection Scheme Tester: N = %d;  CFL = %1.3f;  M = %d;\n\n',N,CFL,M);
    
    time = 0;
    dt   = CFL * min( h/2/max(abs(u(:))) , h/2/max(abs(w(:))) );
    ti   = 1;
    
    while max([max(afr(:))-min(afr(:)),max(au3(:))-min(au3(:))]) < 1.5*a1 && ti <= M
        
        % store previous advection rates
        dadtcdo = dadtcd;
        dadtfdo = dadtfd;
        dadtfro = dadtfr;
        dadtu1o = dadtu1;
        dadtu2o = dadtu2;
        dadtu3o = dadtu3;
        
        dadtcd  = advection(u,w,acd,h,a0,'centred FD');
        dadtfd  = advection(u,w,afd,h,a0,'flxdiv');
        dadtfr  = advection(u,w,afr,h,a0,'fromm');
        dadtu1  = advection(u,w,au1,h,a0,'first upwind');
        dadtu2  = advection(u,w,au2,h,a0,'second upwind');
        dadtu3  = advection(u,w,au3,h,a0,'third upwind');
        
        % advance advection time step
        acd(2:end-1,2:end-1)    = acd(2:end-1,2:end-1)  - (theta.*dadtcd+(1-theta).*dadtcdo)*dt;
        afd(2:end-1,2:end-1)    = afd(2:end-1,2:end-1)  - (theta.*dadtfd+(1-theta).*dadtfdo)*dt;
        afr(2:end-1,2:end-1)    = afr(2:end-1,2:end-1)  - (theta.*dadtfr+(1-theta).*dadtfro)*dt;
        au1(2:end-1,2:end-1)    = au1(2:end-1,2:end-1)  - (theta.*dadtu1+(1-theta).*dadtu1o)*dt;
        au2(2:end-1,2:end-1)    = au2(2:end-1,2:end-1)  - (theta.*dadtu2+(1-theta).*dadtu2o)*dt;
        au3(2:end-1,2:end-1)    = au3(2:end-1,2:end-1)  - (theta.*dadtu3+(1-theta).*dadtu3o)*dt;
        
        % advance reference gaussian
        x0   = x0 + interp2(xa2d,za2d,umid,x0,z0)*dt;
        z0 = z0 + interp2(xa2d,za2d,wmid,x0,z0)*dt;
        gsn  = exp(- (xa2d-x0).^2./radius.^2 - (za2d-z0).^2./radius.^2 );
        
        res_acd = acd - (a0 + a1.*gsn);
        res_afd = afd - (a0 + a1.*gsn);
        res_afr = afr - (a0 + a1.*gsn);
        res_au1 = au1 - (a0 + a1.*gsn);
        res_au2 = au2 - (a0 + a1.*gsn);
        res_au3 = au3 - (a0 + a1.*gsn);
        
        resnorm_acd = norm(res_acd(:),2)./norm(a0+a1.*gsn,2);
        resnorm_afd = norm(res_afd(:),2)./norm(a0+a1.*gsn,2);
        resnorm_afr = norm(res_afr(:),2)./norm(a0+a1.*gsn,2);
        resnorm_au1 = norm(res_au1(:),2)./norm(a0+a1.*gsn,2);
        resnorm_au2 = norm(res_au2(:),2)./norm(a0+a1.*gsn,2);
        resnorm_au3 = norm(res_au3(:),2)./norm(a0+a1.*gsn,2);
        
        ti = ti+1;
        time = time +dt;
        
    end
    fprintf(1,'   --- %d; cntrfd %1.3e; flxdiv %1.3e; fromm %1.3e; upw1 %1.3e; upw2 %1.3e; upw3 %1.3e;\n',ti,resnorm_acd,resnorm_afd,resnorm_afr,resnorm_au1,resnorm_au2,resnorm_au3);
    figure(1)
    plot(x0,z0,'*')
    hold on
    resend_acd=[resend_acd resnorm_acd];
    resend_afd=[resend_afd resnorm_afd];
    resend_afr=[resend_afr resnorm_afr];
    resend_au1=[resend_au1 resnorm_au1];
    resend_au2=[resend_au2 resnorm_au2];
    resend_au3=[resend_au3 resnorm_au3];
end

% profile report
% saveas(figure(3),[pwd '/out/', RunID, '/ time residuals ', num2str(ti), 'iterations.jpg'])
figure(2)
loglog(CFL_all,resend_acd,'-*r')
hold on
loglog(CFL_all,resend_afd,'-*b')
hold on
loglog(CFL_all,resend_afr,'-*k')
hold on
loglog(CFL_all,resend_au1,'-*c')
hold on
loglog(CFL_all,resend_au2,'-*g')
hold on
loglog(CFL_all,resend_au3,'-*y')
hold on
legend('Central difference','flux divergent','fromm','first upwind','second upwind','third upwind')
xlabel('CFL'); ylabel('Residuals');

function adv = advection(u,w,a,h,a0,output)

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

agh([1 2 end-1 end],:) = a0;
agh(:,[1 2 end-1 end]) = a0;

acc = agh(3:end-2,3:end-2);
ajp = agh(4:end-1,3:end-2);  ajpp = agh(5:end-0,3:end-2);
ajm = agh(2:end-3,3:end-2);  ajmm = agh(1:end-4,3:end-2);
aip = agh(3:end-2,4:end-1);  aipp = agh(3:end-2,5:end-0);
aim = agh(3:end-2,2:end-3);  aimm = agh(3:end-2,1:end-4);

switch output
    
    case 'centred FD'
        
        adv = vx.*(aip-aim)./h + vz.*(ajp-ajm)./h;
        
    case 'flxdiv'
        
        adv = ((ajp+acc)./2.*wp - (ajm+acc)./2.*wm)./h ...
            + ((aip+acc)./2.*up - (aim+acc)./2.*um)./h ...
            - acc.*(diff(w(:,2:end-1),1,1)./h + diff(u(2:end-1,:),1,2)./h);
        
    case 'fromm'
        
        adv   =     up .*(-(aipp-aip)./h./8 + (aip + acc)./h./2 + (acc-aim )./h./8) ...
            - abs(up).*(-(aipp-aip)./h./8 + (aip - acc)./h./4 - (acc-aim )./h./8) ...
            -     um .*(-(aip -acc)./h./8 + (acc + aim)./h./2 + (aim-aimm)./h./8) ...
            + abs(um).*(-(aip -acc)./h./8 + (acc - aim)./h./4 - (aim-aimm)./h./8) ...
            +     wp .*(-(ajpp-ajp)./h./8 + (ajp + acc)./h./2 + (acc-ajm )./h./8) ...
            - abs(wp).*(-(ajpp-ajp)./h./8 + (ajp - acc)./h./4 - (acc-ajm )./h./8) ...
            -     wm .*(-(ajp -acc)./h./8 + (acc + ajm)./h./2 + (ajm-ajmm)./h./8) ...
            + abs(wm).*(-(ajp -acc)./h./8 + (acc - ajm)./h./4 - (ajm-ajmm)./h./8) ...
            - acc.*(diff(w(:,2:end-1),1,1)./h + diff(u(2:end-1,:),1,2)./h);
        
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
