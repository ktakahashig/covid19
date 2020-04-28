% Calcular la componente anual del estimado del numero de reproduccion
% segun Wang et al usando datos de Campo de Marte
% @ktakahashig
% 27 abril 2020
%

clear

% Cargar datos horarios
dat=load('-ascii','campo_de_marte.txt');
tzero=datenum(2020,1,1); 
Rv=461; 
yy=dat(:,3)+2000; mm=dat(:,1); dd=dat(:,2); hh=dat(:,4);
t=datenum(yy,mm,dd,hh); 
T=dat(:,6);  % Temperatura  (C)
RH=dat(:,8);  % Humedad relativa (%)
e=es(T+273,15).*RH/100; % Presion de vapor (hPa)
AH=1e5*e./Rv./(T+273.15); % Humidad absoluta (g/kg)

% Calcular promedios diarios
t1=floor(t(1)); t2=floor(t(end));
Td=[]; RHd=[]; AHd=[];
for i=1:t2-t1+1
   ii=find(t>=t1+i-1&t<t1+i);
   if (length(ii)>0)
      Td=[Td; nanmean(T(ii))];
      RHd=[RHd; nanmean(RH(ii))];
      AHd=[AHd; nanmean(AH(ii))];
   else
      Td=[Td; NaN];
      RHd=[RHd; NaN];
      AHd=[AHd; NaN];
   end
end
T=Td; RH=RHd; AH=AHd;
t=[t1:t2]'; tt=t-tzero;

% Numero de reproduccion segun formulas 1 y 2 de Wang et al
R1 = 3.011 - 0.0233 * T - 0.0133 * RH;
R2 = 2.268 - 0.0704 * AH;

% Ajuste armonico
X=[cos(2*pi*tt/365.24) sin(2*pi*tt/365.24) ones(length(tt),1)];
AT=regress(T,X); ARH=regress(RH,X); AAH=regress(AH,X);
AR1=regress(R1,X); AR2=regress(R2,X);

% Reporte de parametros de beta = A*cos(2*pi*(t-t0)/365.24) 
% considerando beta=gamma*R, donde t=0 es tzero
gamma=0.111; % Mejor ajuste
tini=datenum(2020,3,4)-tzero;
A=gamma*sqrt(AR1(1)^2+AR1(2)^2); 
t0=atan2(AR1(2),AR1(1))*365.24/2/pi; if (t0<0); t0=t0+365.24; end
B=-A*cos(2*pi*(tini-t0)/365.24);
disp('Formula 1:')
disp(['A = ' num2str(A)])
disp(['t0 = ' num2str(t0)]) 
disp(['B = ' num2str(B)]) 
disp(' ')
disp('Formula 2:')
A=gamma*sqrt(AR2(1)^2+AR2(2)^2);  
t0=atan2(AR2(2),AR2(1))*365.24/2/pi; if (t0<0); t0=t0+365.24; end
B=-A*cos(2*pi*(tini-t0)/365.24);
disp(['A = ' num2str(A)])
disp(['t0 = ' num2str(t0)]) 
disp(['B = ' num2str(B)]) 


% Para graficar
t2=[t(1):datenum(2021,1,1)]'; tt2=t2-tzero;
X2=[cos(2*pi*tt2/365.24) sin(2*pi*tt2/365.24) ones(length(tt2),1)];
That=X2*AT; RHhat=X2*ARH; AHhat=X2*AAH; 
R1hat=X2*AR1; R2hat=X2*AR2; 

figure(1,'papersize',[12 10])
%
subplot(3,1,1)
plot(t,T,'b','linewidth',1,t2,That,'r','linewidth',3);
xlim([t2(1) t2(end)])
set(gca,'XTick',datenum([2017 2017 2018 2018 2019 2019 2020 2020 2021]',[1 7 1 7 1 7 1 7 1]',ones(9,1)));
datetick('mm/yy','keeplimits','keepticks')
grid on
title('a\) Temperatura diaria (C)','fontsize',15)
set(gca,'fontsize',12)
%
subplot(3,1,2)
plot(t,RH,'b','linewidth',1,t2,RHhat,'r','linewidth',3);
xlim([t2(1) t2(end)])
set(gca,'XTick',datenum([2017 2017 2018 2018 2019 2019 2020 2020 2021]',[1 7 1 7 1 7 1 7 1]',ones(9,1)));
datetick('mm/yy','keeplimits','keepticks')
grid on
title('b\) Humedad relativa diaria (%)','fontsize',15) 
set(gca,'fontsize',12)
%
subplot(3,1,3)
plot(t,AH,'b','linewidth',1,t2,AHhat,'r','linewidth',3);
xlim([t2(1) t2(end)])
set(gca,'XTick',datenum([2017 2017 2018 2018 2019 2019 2020 2020 2021]',[1 7 1 7 1 7 1 7 1]',ones(9,1)));
datetick('mm/yy','keeplimits','keepticks')
grid on
title('c\) Humedad absoluta diaria (g/m^3)','fontsize',15) 
set(gca,'fontsize',12)
%
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition", 
     [border, border, (papersize - 2*border)],...
     'paperorientation','portrait')
print('-dpdf','campo_de_marte.pdf')


figure(2,'papersize',[9 10])
subplot(2,1,1)
plot(t,R1,'b','linewidth',1,t2,R1hat,'r','linewidth',3);
title('a\) Numero efectivo de reproduccion Re de COVID-19 (Formula 1 de Wang et al, 2020)','fontsize',15)
xlim([t2(1) t2(end)])
ylim([0.9 1.6])
set(gca,'XTick',datenum([2017 2017 2018 2018 2019 2019 2020 2020 2021]',[1 7 1 7 1 7 1 7 1]',ones(9,1)));
datetick('mm/yy','keeplimits','keepticks')
grid on
set(gca,'fontsize',12)
%
subplot(2,1,2)
plot(t,R2,'b','linewidth',1,t2,R2hat,'r','linewidth',3);
title('b\) Numero efectivo de reproduccion Re de COVID-19 (Formula 2 de Wang et al, 2020)','fontsize',15)
xlim([t2(1) t2(end)])
ylim([0.9 1.6])
set(gca,'XTick',datenum([2017 2017 2018 2018 2019 2019 2020 2020 2021]',[1 7 1 7 1 7 1 7 1]',ones(9,1)));
datetick('mm/yy','keeplimits','keepticks')
grid on
set(gca,'fontsize',12)
%
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition", 
     [border, border, (papersize - 2*border)],...
     'paperorientation','landscape')
print('-dpdf','R_wang_et_al.pdf')

