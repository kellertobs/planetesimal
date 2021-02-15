clear all
close all
L   = 1e3;
N  = 200;
h  = L/N;
movement    = 'rotation';
RunID       = 'central gaussian rotation'
% nodal grid
x   = -L/2+h/2:h:(L/2 - h/2);
z   = -L/2+h/2:h:(L/2 - h/2);
[x2d,z2d] = meshgrid(x,z);
%full nodal grid
xa  = -L/2-h/2:h:L/2+h/2;
za  = -L/2-h/2:h:L/2+h/2;
[xa2d,za2d] = meshgrid(xa,za);
% staggered grids for vx
xu  = -L/2:h:L/2;
zu  = -L/2-h/2:h:L/2+h/2;
[xu2d,zu2d] = meshgrid(xu,zu);
% staggered grids for vz
xw  = -L/2-h/2:h:L/2+h/2;
zw  = -L/2:h:L/2;
[xw2d,zw2d] = meshgrid(xw,zw);
% arrays and quantities
u   = zeros(N+2,N+1); % staggered vx
w   = zeros(N+1,N+2); % staggered vz
a   = zeros(N+2,N+2);     % advected quantity

% set initial distribution for advected quantity
a0 = 1000;
a   = a + a0;

radius = L/10;
a = a + 500.*exp(- (xa2d-L/6).^2./radius.^2 - (za2d-L/6).^2./radius.^2 ); 

figure(1)
subplot(2,2,1)
imagesc(x,z,a(2:end-1,2:end-1)); 
colorbar
title('initial')
switch movement
    case 'rotation'
        u = -zu2d;
        w = xw2d;  
end
umid = (u(:,1:end-1)+u(:,2:end))/2;
wmid = (w(1:end-1,:)+w(2:end,:))/2;


af      = a;
au1 = a;
au3     = a;
au2     = a;
nt = 2000;
time = 0;
ti = 0

while max([max(max(af))-min(min(af)) max(max(au2))-min(min(au2)) max(max(au3))-min(min(au3))]) < 600
    ti = ti+1;
% for ti = 1:1:nt
    
    dt = 0.1 * min( h/2/max(abs(u(:))) , h/2/max(abs(w(:))) );    
    
    dadtf   = advection(u,w,af,h,N,'fromm');
    dadtu1  = advection(u,w,au1,h,N,'first upwind');
    dadtu3  = advection(u,w,au3,h,N,'third upwind');
    dadtu2  = advection(u,w,au2,h,N,'second upwind');
    
    if ti==1
    af(2:end-1,2:end-1)     = af(2:end-1,2:end-1)   - dadtf*dt;
    au3(2:end-1,2:end-1)    = au3(2:end-1,2:end-1)  - dadtu3*dt;
    au2(2:end-1,2:end-1)    = au2(2:end-1,2:end-1)  - dadtu2*dt;
    au1(2:end-1,2:end-1)    = au1(2:end-1,2:end-1)  - dadtu1*dt;
    else        
    af(2:end-1,2:end-1)     = af(2:end-1,2:end-1)   - (dadtf+dadtf0)/2*dt;
    au3(2:end-1,2:end-1)    = au3(2:end-1,2:end-1)  - (dadtu3+dadtu30)/2*dt;
    au2(2:end-1,2:end-1)    = au2(2:end-1,2:end-1)  - (dadtu2+dadtu20)/2*dt;
    au1(2:end-1,2:end-1)    = au1(2:end-1,2:end-1)  - dadtu1*dt;
    end
    dadtf0 = dadtf;
    dadtu30 = dadtu3;
    dadtu20 = dadtu2;
    
    if ~mod(ti,20)    
        figure(1)
    subplot(2,2,1)
    imagesc(x,z,af(2:end-1,2:end-1)); 
    colorbar
    title('fromm scheme')
    
    subplot(2,2,2)
    imagesc(x,z,au1(2:end-1,2:end-1)); 
    colorbar
    title('1st order upwind')
    
    subplot(2,2,4)
    imagesc(x,z,au3(2:end-1,2:end-1)); 
    colorbar
    title('3rd order upwind')
    
    subplot(2,2,3)
    imagesc(x,z,au2(2:end-1,2:end-1)); 
    colorbar
    title('2nd order upwind')
    
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
    
    pause(0.1)
%     saveas(figure(1),[pwd '/out/', RunID, ' ', num2str(ti), 'iterations.jpg'])
%     saveas(figure(2),[pwd '/out/dadt ', RunID, ' ', num2str(ti), 'iterations.jpg'])
    end
    time = time +dt;
end



function adv = advection(u,w,a,h,N,output)
wp = w(2:end,2:end-1);
wm = w(1:end-1,2:end-1);
up = u(2:end-1,2:end);
um = u(2:end-1,1:end-1);

 

agh                    = zeros(size(a)+2);
agh(2:end-1,2:end-1)   = a;

 

agh([1 2 end-1 end],:) = agh([3 3 end-2 end-2],:);
agh(:,[1 2 end-1 end]) = agh(:,[3 3 end-2 end-2]);


 

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
        adv = zeros(size(a)-2);
        for i = 1:1:N
        for j = 1:1:N
            vx = (up(i,j)+um(i,j))/2;
            vz = (wp(i,j)+wm(i,j))/2;

            if vx>0
                dadx = vx*(agh(i+2,j+2)-agh(i+2,j+1))/h;
            else
                dadx = vx*(-agh(i+2,j+2)+agh(i+2,j+3))/h;
            end
            if vz>0
                dadz = vz*(agh(i+2,j+2)-agh(i+1,j+2))/h;
            else
                dadz = vz*(-agh(i+2,j+2)+agh(i+3,j+2))/h;
            end
            adv(i,j) = dadx+dadz;
        end
        end
        
    case 'second upwind'
        adv = zeros(size(a)-2);
        for i = 1:1:N
        for j = 1:1:N
            vx = (up(i,j)+um(i,j))/2;
            vz = (wp(i,j)+wm(i,j))/2;

            if vx>0
                dadx = vx*(3*agh(i+2,j+2)-4*agh(i+2,j+1)+agh(i+2,j))/2/h;
            else
                dadx = vx*(-3*agh(i+2,j+2)+4*agh(i+2,j+3)-agh(i+2,j+4))/2/h;
            end
            if vz>0
                dadz = vz*(3*agh(i+2,j+2)-4*agh(i+1,j+2)+agh(i,j+2))/2/h;
            else
                dadz = vz*(-3*agh(i+2,j+2)+4*agh(i+3,j+2)-agh(i+4,j+2))/2/h;
            end
            adv(i,j) = dadx+dadz;
        end
        end
        
    case 'third upwind'
        adv = zeros(size(a)-2);
        for i = 1:1:N
        for j = 1:1:N
            vx = (up(i,j)+um(i,j))/2;
            vz = (wp(i,j)+wm(i,j))/2;

            if vx>0
                dadx = vx*(2*agh(i+2,j+3)+3*agh(i+2,j+2)-6*agh(i+2,j+1)+agh(i+2,j))/6/h;
            else
                dadx = vx*(-2*agh(i+2,j+1)-3*agh(i+2,j+2)+6*agh(i+2,j+3)-agh(i+2,j+4))/6/h;
            end
            if vz>0
                dadz = vz*(2*agh(i+3,j+2)+3*agh(i+2,j+2)-6*agh(i+1,j+2)+agh(i,j+2))/6/h;
            else
                dadz = vz*(-2*agh(i+1,j+2)-3*agh(i+2,j+2)+6*agh(i+3,j+2)-agh(i+4,j+2))/6/h;
            end
            adv(i,j) = dadx+dadz;
        end
        end        
end
  

  
end
