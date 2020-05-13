% Genera los escenarios de R0 considerando medidas de gobierno
% y la posibles estacionalidad
% @ktakahashig
% 29 abril 2020
%
clear

tzero=datenum(2020,1,1); % t = 0
tini=datenum(2020,3,4); % Inicio de brote
t=[tini:tini+365]'; % Periodo total del brote en 2020

gamma = 0.111; % Tasa de recuperacion+fallecimiento estimado en el ajuste del modelo

Rmin1=0.8; % Minimo alcanzado (en forma sostenida) por Alemania, Francia y España
Rmin2=1.2; % Valor en Peru segun CMMID es 1.3 al 24/04/2020 (datos hasta 04/05). 
           % El Presidente Vizcarra anuncio R entre 1.1 y 1.2 el 08/05/2020


% Parametros de estacionalidad de beta
% 1: Kissler et al, 2: Neher et al, 3: Wang et al (form. 1), 4: Wang et al (form.2)
A=[0.0256 0.024427 0.00440504 0.020603 0.038796]';
t0=[196   183      248.549    234.9    166.021]';
B=[0.0168 0.011560 0.0044447  0.020254 0.0077629]';


% Numero de reproduccion (CMMID)
load Re.txt
tRe=datenum(Re(:,[1:3])); Re=Re(:,4);
ARe=regress(Re,[tRe-tzero ones(length(Re),1)]);
disp(['R0 alcanza ' num2str(Rmin1) ' en t = ' num2str(round((Rmin1-ARe(2))/ARe(1))) ])
datevec(round((Rmin1-ARe(2))/ARe(1))+tzero)
disp(['R0 alcanza ' num2str(Rmin2) ' en t = ' num2str(round((Rmin2-ARe(2))/ARe(1))) ])
datevec(round((Rmin2-ARe(2))/ARe(1))+tzero)


% Escenarios base de Re
tt=t-tzero;
R0base=[tt ones(length(t),1)]*ARe; R0base(find(tt>(Rmin1-ARe(2))/ARe(1)))=Rmin1;
R0basedum=[tt ones(length(t),1)]*ARe; R0basedum(find(tt>(Rmin2-ARe(2))/ARe(1)))=Rmin2;
R0base=[R0base R0basedum];
Nbase=size(R0base,2);

% Variaciones estacionales de beta y R0
betaest=[];
Nest=length(A);
for i=1:Nest
   betaest=[betaest A(i).*cos(2*pi*(tt-t0(i))/365.24)+B(i)];
end
R0est=betaest/gamma;

% Guardar archivos
save -ascii R0base.txt R0base
save -ascii R0est.txt R0est
save -ascii betaest.txt betaest


figure(1,'papersize',[9 4.5]),clf

labs=strvcat('A','B','C','D');
option=strvcat('r-.','g:','c--','b--','m-');
for j=1:Nbase
subplot(1,Nbase,j)
   h=plot(t,R0base(:,j),'k-');
   hold on
   for i=1:Nest
      h=[h;plot( t,R0base(:,j)+R0est(:,i),option(i,:))];
   end
plot([t(1) t(end)],[1 1],'k-','linewidth',0.8);
hold off
xlabel('Mes')
xlim([t(1) t(end)]);
ylim([0.5 2.5]);
set(h(1),'linewidth',2.8);
set(h(2:end),'linewidth',2.5);
set(h(end),'linewidth',2);
set(gca,'XTick',datenum([repmat(2020,10,1); repmat(2021,4,1)],[[3:12] [1:4]]',1))
datetick('mm','keeplimits','keepticks')
set(gca,'fontsize',12)
title(['R0: Escenarios ' labs(j)],'fontsize',14)
end
legend('Base','Kissler','Neher','Wang 1','Wang 2','Merow&Urban','location','northeast')

papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition", 
     [border, border, (papersize - 2*border)],...
     'paperorientation','landscape')
print('-dpdf','escenarios.pdf')
