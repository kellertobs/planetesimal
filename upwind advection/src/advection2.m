function [adva,advb] = advection2(vx_out,vz_out,a,b,dx,dz,output)
w = vz_out(1:end-1,:); % vz
u = vx_out(:,1:end-1); % vx

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

% agh([1 2 end-1 end],:) = agh([4 3 end-2 end-3],:);
% agh(:,[1 2 end-1 end]) = agh(:,[4 3 end-2 end-3]);
agh([1 2 end-1 end],:) = agh([3 3 end-2 end-2],:);
agh(:,[1 2 end-1 end]) = agh(:,[3 3 end-2 end-2]);

acc = agh(3:end-2,3:end-2);
ajp = agh(4:end-1,3:end-2);  ajpp = agh(5:end-0,3:end-2);
ajm = agh(2:end-3,3:end-2);  ajmm = agh(1:end-4,3:end-2);
aip = agh(3:end-2,4:end-1);  aipp = agh(3:end-2,5:end-0);
aim = agh(3:end-2,2:end-3);  aimm = agh(3:end-2,1:end-4);

bgh                    = zeros(size(b)+2);
bgh(2:end-1,2:end-1)   = b;

% agh([1 2 end-1 end],:) = agh([4 3 end-2 end-3],:);
% agh(:,[1 2 end-1 end]) = agh(:,[4 3 end-2 end-3]);
bgh([1 2 end-1 end],:) = bgh([3 3 end-2 end-2],:);
bgh(:,[1 2 end-1 end]) = bgh(:,[3 3 end-2 end-2]);

bcc = bgh(3:end-2,3:end-2);
bjp = bgh(4:end-1,3:end-2);  bjpp = bgh(5:end-0,3:end-2);
bjm = bgh(2:end-3,3:end-2);  bjmm = bgh(1:end-4,3:end-2);
bip = bgh(3:end-2,4:end-1);  bipp = bgh(3:end-2,5:end-0);
bim = bgh(3:end-2,2:end-3);  bimm = bgh(3:end-2,1:end-4);

switch output
    
    case 'centred FD'
        % not recommended
        adva = vx.*(aip-aim)./h + vz.*(ajp-ajm)./h;
        
    case 'flxdiv'      
        adva = ((ajp+acc)./2.*wp - (ajm+acc)./2.*wm)./dz ...
            + ((aip+acc)./2.*up - (aim+acc)./2.*um)./dx ...
            - acc.*(diff(w(:,2:end-1),1,1)./dz + diff(u(2:end-1,:),1,2)./dx);
        
        advb = ((bjp+bcc)./2.*wp - (bjm+acc)./2.*wm)./dz ...
            + ((bip+bcc)./2.*up - (bim+bcc)./2.*um)./dx ...
            - bcc.*(diff(w(:,2:end-1),1,1)./dz + diff(u(2:end-1,:),1,2)./dx);
        
    case 'fromm'        
        adva   =     up .*(-(aipp-aip)./dx./8 + (aip + acc)./dx./2 + (acc-aim )./dx./8) ...
            - abs(up).*(-(aipp-aip)./dx./8 + (aip - acc)./dx./4 - (acc-aim )./dx./8) ...
            -     um .*(-(aip -acc)./dx./8 + (acc + aim)./dx./2 + (aim-aimm)./dx./8) ...
            + abs(um).*(-(aip -acc)./dx./8 + (acc - aim)./dx./4 - (aim-aimm)./dx./8) ...
            +     wp .*(-(ajpp-ajp)./dz./8 + (ajp + acc)./dz./2 + (acc-ajm )./dz./8) ...
            - abs(wp).*(-(ajpp-ajp)./dz./8 + (ajp - acc)./dz./4 - (acc-ajm )./dz./8) ...
            -     wm .*(-(ajp -acc)./dz./8 + (acc + ajm)./dz./2 + (ajm-ajmm)./dz./8) ...
            + abs(wm).*(-(ajp -acc)./dz./8 + (acc - ajm)./dz./4 - (ajm-ajmm)./dz./8) ...
            - acc.*(diff(w(:,2:end-1),1,1)./dz + diff(u(2:end-1,:),1,2)./dx);
        
        advb   =     up .*(-(aipp-aip)./dx./8 + (bip + bcc)./dx./2 + (bcc-bim )./dx./8) ...
            - abs(up).*(-(bipp-bip)./dx./8 + (bip - bcc)./dx./4 - (bcc-bim )./dx./8) ...
            -     um .*(-(bip -bcc)./dx./8 + (bcc + bim)./dx./2 + (bim-bimm)./dx./8) ...
            + abs(um).*(-(bip -bcc)./dx./8 + (bcc - bim)./dx./4 - (bim-bimm)./dx./8) ...
            +     wp .*(-(bjpp-bjp)./dz./8 + (bjp + bcc)./dz./2 + (bcc-bjm )./dz./8) ...
            - abs(wp).*(-(bjpp-bjp)./dz./8 + (bjp - bcc)./dz./4 - (bcc-bjm )./dz./8) ...
            -     wm .*(-(bjp -bcc)./dz./8 + (bcc + bjm)./dz./2 + (bjm-bjmm)./dz./8) ...
            + abs(wm).*(-(bjp -bcc)./dz./8 + (bcc - bjm)./dz./4 - (bjm-bjmm)./dz./8) ...
            - bcc.*(diff(w(:,2:end-1),1,1)./dz + diff(u(2:end-1,:),1,2)./dx);
        
    case 'first upwind'       
        axp   = (aip-acc)./dx;
        axm   = (acc-aim)./dx;
        azp   = (ajp-acc)./dz;
        azm   = (acc-ajm)./dz;
        
        daxdt = vxp.*axm + vxm.*axp;
        dazdt = vzp.*azm + vzm.*azp;
        
        adva   = daxdt + dazdt;
        
        bxp   = (bip-bcc)./dx;
        bxm   = (bcc-bim)./dx;
        bzp   = (bjp-bcc)./dz;
        bzm   = (bcc-bjm)./dz;
        
        dbxdt = vxp.*bxm + vxm.*bxp;
        dbzdt = vzp.*bzm + vzm.*bzp;
        
        advb   = dbxdt + dbzdt;
        
    case 'second upwind'        
        axp   = (-3*acc+4*aip-aipp)/2/dx;
        axm   = ( 3*acc-4*aim+aimm)/2/dx;
        azp   = (-3*acc+4*ajp-ajpp)/2/dz;
        azm   = ( 3*acc-4*ajm+ajmm)/2/dz;
        
        daxdt = vxp.*axm + vxm.*axp;
        dazdt = vzp.*azm + vzm.*azp;
        
        adva   = daxdt + dazdt;
        
        bxp   = (-3*bcc+4*bip-bipp)/2/dx;
        bxm   = ( 3*bcc-4*bim+bimm)/2/dx;
        bzp   = (-3*bcc+4*bjp-bjpp)/2/dz;
        bzm   = ( 3*bcc-4*bjm+bjmm)/2/dz;
        
        dbxdt = vxp.*bxm + vxm.*bxp;
        dbzdt = vzp.*bzm + vzm.*bzp;
        
        advb   = dbxdt + dbzdt;
        
    case 'third upwind'        
        axp   = (-2*aim-3*acc+6*aip-aipp)/6/dx;
        axm   = ( 2*aip+3*acc-6*aim+aimm)/6/dx;
        azp   = (-2*ajm-3*acc+6*ajp-ajpp)/6/dz;
        azm   = ( 2*ajp+3*acc-6*ajm+ajmm)/6/dz;
        
        daxdt = vxp.*axm + vxm.*axp;
        dazdt = vzp.*azm + vzm.*azp;
        
        adva   = daxdt + dazdt;
        
        bxp   = (-2*bim-3*bcc+6*bip-bipp)/6/dx;
        bxm   = ( 2*bip+3*bcc-6*bim+bimm)/6/dx;
        bzp   = (-2*bjm-3*bcc+6*bjp-bjpp)/6/dz;
        bzm   = ( 2*bjp+3*bcc-6*bjm+bjmm)/6/dz;
        
        dbxdt = vxp.*bxm + vxm.*bxp;
        dbzdt = vzp.*bzm + vzm.*bzp;
        
        advb   = dbxdt + dbzdt;
        
end
end