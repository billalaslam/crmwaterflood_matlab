%% ------- Capacitance Resistance Model Inj-Prod (CRMIP) ------- %%
% Author      : BMA (bmaslam@tm.itb.ac.id), Ogrindo ITB
% Description : CRM-IP model implementation in Matlab for Waterflood
%               Predictive Model (calibration and forecast mode)
%               Main Reference from Sayarpour (2009)
% Feature     : History matching by fmincon solver
%               Optimization mode:
%               1. Cum Oil
close all;
clear;
clc;
%% load data
PRDI=readtable('PRD_Var.xlsx');
INJI=readtable('INJ_Var.xlsx');
prd=table2array(PRDI(1:end,2:end)); %set time to start & num of prod
inj=table2array(INJI(1:end,2:end)); %set time to start & num of injt
XYII=readtable('XYI.csv');
XYPI=readtable('XYP.csv');
xyi=table2array(XYII(:,2:end));
xyp=table2array(XYPI(2:end,2:end));


ninj=size(inj,2);
npro=size(prd,2);
%% Plot history data
t=1:1:size(prd,1);
plot(t,sum(inj,2),'LineWidth',2)
hold on
scatter(t,sum(prd,2))
legend('total inj rate','total liq. prod. rate','Location','best')
ax=gca;
ax.YAxis.Exponent = 0;
hold off;

xlabel('time (month)')
ylabel('Liquid injection - production rate (BLPD)')

figure
plot(t,sum(inj,2)./sum(prd,2),'LineWidth',2)
xlabel('time (month)')
ylabel('VRR')
ylim([0 Inf])


%% Guess Initial Value of Connectivity by Well Pair Distance
close all
distance=zeros(size(xyi,1), size(xyp,1));
for k=1:size(xyp,1)
    for i =1:size(xyi,1)
        distance(i,k)=((abs(xyp(k,1)-xyi(i,1)).^2+(abs(xyp(k,2)-xyi(i,2)).^2))).^.5;
    end
end
%% ---------Creating the matrix of the connectivity-----------
fmat=zeros*distance;
for k=1:size(xyp,1)
    for i =1:size(xyi,1)
        fmat(i,k)=(1/distance(i,k))./(sum(1./distance(:,k)));
    end
end
%% ---------Creating the matrix of taos------------------------
tao=0.6*ones(size(inj,2),size(prd,2));
%qest=zeros*prd;

%% ---------CRM History Match----------------------------------

%select model
fun=@crmip_hm;

%set training interval
tstart=30;
dt=30;

%HM solver input
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
%--Design Variable
fhm=zeros(ninj,1);
taohm=zeros(ninj,1);
zhm=[fhm;taohm];

fsol=zeros(ninj,npro);
taosol=zeros(ninj,npro);
zsol=[fsol;taosol];


%--Upper & Lower Bound
lb_tao=zeros(ninj,1);
ub_tao=Inf(ninj,1); %ub_tao=Inf(ninj,1); %ub_tao=10000*ones(ninj,1);

lb_f=zeros(ninj,1);
ub_f=ones(ninj,1);

lb_z=[lb_f;lb_tao];
ub_z=[ub_f;ub_tao];

%--Constraints

nonlcon=[];
A_fmat=ones(1,ninj);
A_tao=zeros(1,ninj);
Az=[A_fmat A_tao];

b=1;

    %% optimize all simultaneously
    clc;
    for pind=1:npro
        fhm0=fmat(:,pind);
        tao0=tao(:,pind);
        z0=[fhm0; tao0];
        zhm=fmincon(@(zhm)msehm1(fun,inj,prd(:,pind),zhm,tstart,dt),z0,[],[],[],[],lb_z,ub_z,nonlcon,options);
        zsol(:,pind)=zhm;
    end
    fsol=zsol(1:ninj,1:npro);
    taosol=zsol(ninj+1:2*ninj,1:npro);

    %% optimize tao
    %%
    %    %     for pind=1:npro
    %         tao0=tao(:,pind);
    %         taohm=fmincon(@(taohm)msehm2(fun,inj,prd(:,pind),fhm,taohm,tstart,dt),tao0,[],[],[],[],lb_tao,[],nonlcon,options);
    %         taosol(:,pind)=taohm;
    %     end
    %
    %% optimize fij
    %%
    %
    %     for pind=1:npro
    %         fhm0=fmat(:,pind);
    %         taohm=taosol(:,pind);
    %         fhm=fmincon(@(fhm)msehm2(fun,inj,prd(:,pind),fhm,taohm,tstart,dt),fhm0,A_fmat,b,[],[],lb_f,ub_f,nonlcon,options);
    %         fsol(:,pind)=fhm;
    %     end
    %
%% ----------- checking and plotting-------------------------
close all
thm=tstart:size(prd,1);
qestf=zeros(size(thm,2),1);
for i =1:npro
    subplot(2,2,i)
    qest=crmip_hm(inj,prd(:,i),fsol(:,i),taosol(:,i),tstart,dt);
    qestf=qestf+qest;
    scatter(1:size(prd,1),prd(:,i),'ro')
    hold on
    plot(tstart:size(prd,1),qest(:,1),'LineWidth',2,'Color','k')
   
    ylim([0 max(prd(:))+1000])
    hold off
    xlabel('time (month)')
    ylabel('rate (bbl/d)')
    legend('Calculated', 'Simulated','Location','best')
    title(['PRD',num2str(i)])
    %lsn=ls2(fmat(:,i),tao(:,i),size(prd,2),size(inj,2),prd(:,i),inj);
end

figure
scatter(1:size(prd,1),sum(prd,2),'ro')
hold on
plot(tstart:size(prd,1),qestf,'LineWidth',2,'Color','k')

ylim([0 max(qestf)+1000])
hold off
xlabel('time (month)')
ylabel('rate (bbl/d)')
legend('Calculated', 'Simulated','Location','best')
title('Total Field Liquid Rate')
%% ---------------sectors routine----------------------------
angles=zeros(size(xyi,1),size(xyp,1));
arrlength=fsol;
axislimit=max(xyi(:));
if max(xyp(:))>axislimit
    axislimit=max(xyp(:));
end
for i = 1 : size(angles,1)%filling per row for each injector
    for j =1:size(angles,2)
        angles(i,j)=atand((xyp(j,2)-xyi(i,2))/(xyp(j,1)-xyi(i,1)));
        if (xyp(j,2)-xyi(i,2))<0 && (xyp(j,1)-xyi(i,1))<0
            angles(i,j)=angles(i,j)+180;
        end
        if (xyp(j,2)-xyi(i,2))==0 && (xyp(j,1)-xyi(i,1))<0
            angles(i,j)=angles(i,j)+180;
        end
        if (xyp(j,2)-xyi(i,2))>0 && (xyp(j,1)-xyi(i,1))<0
            angles(i,j)=angles(i,j)+180;
        end
    end
end

%% Quiver
xpro=xyp(:,1);
ypro=xyp(:,2);
xinj=xyi(:,1);
yinj=xyi(:,2);

labelinj=cell(1,ninj);
labelpro=cell(1,npro);
for i=1:ninj
    label=['I' num2str(i)];
    labelinj{i}=label;
end

u=zeros(ninj,npro);
v=zeros(ninj,npro);
for i=1:npro
    u(:,i)=fsol(:,i).*cos(deg2rad(angles(:,i)));
    v(:,i)=fsol(:,i).*sin(deg2rad(angles(:,i)));
    
    label=['P' num2str(i)];
    labelpro{i}=label;
end 


figure
hold on;
scatter(xpro,ypro,'filled')
labelpoints(xpro, ypro, labelpro, 'SE', 0.5);
scatter(xinj,yinj,'v','filled')
labelpoints(xinj, yinj, labelinj, 'SE', 0.5);
for i=1:npro
   
    q=quiver(xinj,yinj,u(:,i),v(:,i),'k');
    q.AutoScaleFactor = 0.5;
    q.LineWidth = 1;
end

axis square

hold off;
%% Function Library

function err=msehm1(fun,inj,prd,z,tstart,dt)
%calculate MSE from estimated liquid rate to observed liquid rate
ninj=size(inj,2);
npro=1;
z=reshape(z,[],2);
fmat=reshape(z(:,1),[npro,ninj]);
tao=reshape(z(:,2),[npro,ninj]);

t=tstart:size(prd,1);

qhm=fun(inj,prd,fmat,tao,tstart,dt);

tind=1:size(t,2);


qhist=prd(t(tind),1);

err=sum((qhist-qhm).^2)/size(t,2);

end


function err=msehm2(fun,inj,prd,fmat,tao,tstart,dt)
%calculate MSE from estimated liquid rate to observed liquid rate
ninj=size(inj,2);
npro=size(prd,2);

fmat=reshape(fmat,[npro,ninj]);
tao=reshape(tao,[npro,ninj]);

t=tstart:size(prd,1);

qhm=fun(inj,prd,fmat,tao,tstart,dt);

tind=1:size(t,2);

qhist=prd(t(tind),1);

err=sum((qhist-qhm).^2)/size(t,2);

end
%% 

function qhm=crmip_hm(inj,prd,fmat,tao,tstart,dt)
%calculate estimated liquid rate per producer from CRMIP equation
%constant BHP only
ninj=size(inj,2);

t=tstart:size(prd,1);
qhm=zeros*prd(t,1); %total liq rate todo: change index for multiple well 
qesti=zeros(size(t,2),ninj); %contrib rate from injector todo: change index for multiple well

%initial rates
for i=1:ninj
    qesti(1,i)=inj(tstart,i)*fmat(i);
end

qhm(1)=sum(qesti(1,:));
for i=2:size(t,2)
    
    for j=1:ninj
        qesti(i,j)=qesti(i-1,j)*exp(-dt/tao(j))+...
            ((1-exp(-dt/tao(j)))*fmat(j)*inj(t(i),j));
    end
    qhm(i)=sum(qesti(i,:));
    
end

end

