N=100000; %number of photons
us=50; %scattering coefficient(1/mm)
ua = 0.2; %absorption coefficient
scatterLimit = 70; % maximum scattering times
g = 0.28;
xin = rand(N,1)*0.5;
yin = rand(N,1)*0.5;
xscattering = zeros(N,scatterLimit);
yscattering = zeros(N,scatterLimit);
zscattering = zeros(N,scatterLimit);
lscattering = zeros(N,scatterLimit);
out = zeros(N,1);
outx = zeros(N,1);
outy = zeros(N,1);
path = zeros(N,1);
for i = 1:N
    ls = -1/us*log(rand(scatterLimit,1));
    xscattering(i,1)=xin(i);
    yscattering(i,1)=yin(i);
    zscattering(i,1)=ls(1);
    fees = 2*pi*rand(scatterLimit,1);
    costhetas = 1/2./g.*(1+g.^2-((1-g.^2)./(1-g+2.*g.*rand(scatterLimit,1))).^2);
    sinthetas = sqrt(1-costhetas.^2);
    for j = 2:scatterLimit
        xlast = xscattering(i,j-1);
        ylast = yscattering(i,j-1);
        zlast = zscattering(i,j-1);
        if ((xlast<0) || (xlast>0.5) || (ylast<0) || (ylast>0.5) || (zlast> 7))
            break
        else
            if zlast<0
            out(i) = 1;
            outx(i) = xlast;
            outy(i) = ylast;
            path(i) = sum(ls(1:j));
            break   
            end
        end
        l = ls(j);
        fee = fees(j);
        costheta = costhetas(j);
        sintheta = sinthetas(j);
        dz = l*costheta;
        dx = l*sintheta*cos(fee);
        dy = l*sintheta*sin(fee);
     
        xscattering(i,j) = xscattering(i,j-1)+dx;
        yscattering(i,j) = yscattering(i,j-1)+dy;
        zscattering(i,j) = zscattering(i,j-1)+dz;
    end
end
xtotal = [xin xscattering];
ytotal = [yin yscattering];
ztotal = [zeros(N,1) zscattering];
w = exp(-ua.*path);
wtotal = sum(w);
isDynamic =( xscattering < 0.3 & xscattering > 0.2 & zscattering <0.1 & zscattering >0);
timesDynamic = sum(isDynamic,2);
figure
plot3(xtotal(5,xtotal(5,:)~=0),ytotal(5,xtotal(5,:)~=0),ztotal(5,xtotal(5,:)~=0),'g')
hold on
plot3(xtotal(6,xtotal(6,:)~=0),ytotal(6,xtotal(6,:)~=0),ztotal(6,xtotal(6,:)~=0),'b')
plot3(xtotal(7,xtotal(7,:)~=0),ytotal(7,xtotal(7,:)~=0),ztotal(7,xtotal(7,:)~=0),'r')
plot3(xtotal(8,xtotal(8,:)~=0),ytotal(8,xtotal(8,:)~=0),ztotal(8,xtotal(8,:)~=0),'y')
plot3(xtotal(9,xtotal(9,:)~=0),ytotal(9,xtotal(9,:)~=0),ztotal(9,xtotal(9,:)~=0),'m')
plot3(xtotal(10,xtotal(10,:)~=0),ytotal(10,xtotal(10,:)~=0),ztotal(10,xtotal(10,:)~=0),'c')
plot3(xtotal(11,xtotal(11,:)~=0),ytotal(11,xtotal(11,:)~=0),ztotal(11,xtotal(11,:)~=0),'k')
grid on




        
        