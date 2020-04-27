% Ajuste de modelo SIR para Peru
% Ken Takahashi (ken.takahashi.guevara@gmail.com, ktakahashig)
% 27-04-2020
clear

% Set time
nyrs=1; % Number of years in run
dt=1; % time step in days

% Model params
N=30e6;
asympt=0.18; % Fraccion de infectados que son asintomaticos (crucero Diamond Princess, )

t=[datenum(2020,3,4):datenum(2020,4,26)]';

% Datos oficiales MINSA
% Numero de casos confirmados (acumulado)
casos=load("casos_confirmados.txt");
tdum=datenum(2020,3,6)+[0:length(casos)-1]';
casos=interp1(tdum,casos,t);
casos(1:2)=0;
dcasos=[0 ;diff(casos)]; % Casos diarios
clear tdum

% Fraccion de reporte segun CRF ajustado (CMMID)
frac=load("fraction_reported.txt");
tdum=datenum(2020,3,15)+frac(:,1)-frac(1,1);
tf=floor(tdum(end));
frac=frac(:,2); frac=interp1(tdum,frac,t);
ii=find(isna(frac)==1&t<tdum(1));
frac(ii)=frac(ii(end)+1); % Completar valores al inicio con el primer valor
clear tdum
% Corregir datos MINSA por subreporte y asintomaticos
casoscorrcfr=cumsum(dcasos*100./frac)/(1-asympt);

% Numero de reproduccion (CMMID)
load Re.txt
tdum=datenum(Re(:,[1:3])); Re=Re(:,4);
tf=min([tf tdum(end)]);
Re=interp1(tdum,Re,t);
clear tdum


% Configurar corridas para ajuste de gamma

t0s=[datenum(2020,3,9):datenum(2020,3,31)]; % Fechas iniciales a probarse
gammafit=[]; % Aqui ira el mejor gamma para cada fecha inicial
dgamma=[-0.2:0.001:0.2]'; % Rango de gammas para explorar en torno al first-guess
rmse=[]; % Aqui se guardara el error RMS para cada dia inicial 

for t0=t0s

   rmsedum=[]; % Aqui se guardara el error de simulacion para cada gamma probado
   gammas = 0.15+dgamma; % First-guess de gamma + rango a explorar

   for j=1:length(gammas)

      gamma=gammas(j);
      SS=[]; II=[]; RR=[]; 

      % Condiciones iniciales 
      I=casoscorrcfr(find(t==t0));
      S=N-I; R=0; X=[S I R];

      % Correr modelo
      for i=1:tf-t0+1
          betadum=Re(i)*gamma; 
          X=[X;rk4(X(end,:),dt,[betadum gamma N])];
      end
      S=X(:,1); I=X(:,2); R=X(:,3);

      % Calcular error RMS
      tdum=t(find(isna(casoscorrcfr)==0));
      [tdum ii jj]=intersect(t0+[1:length(S)]',t(find(isna(casoscorrcfr)==0)));
      rmsedum=[rmsedum; mean((log(R(ii)+I(ii))-log(casoscorrcfr(jj))).^2)];

   end

   rmse=[rmse rmsedum]; % Guardar el error RMS
   gammafit=[gammafit gammas(find(rmsedum==min(rmsedum)))]; % Elegir gamma con minimo RMSE
   if (length(gammafit)>1) % Usar gamma elegido o extrapolado como first guess
      gammas = gammafit(end)+(gammafit(end)-gammafit(end-1))+dgamma; 
   else
      gammas = gammafit(end)+dgamma; % Usar gamma elegido como first guess
   end
 
end

t0=datenum(2020,3,15); % Elegido al producir el error RMS localmente minimo
gamma=gammafit(find(t0s==t0)) % gamma=0.111


% Figura para asegurar que para cada fecha el RMSE minimo esta en el rango
figure(3,'papersize',[7 7]) 
[cs h]=contour(t0s,dgamma,rmse,[0.1 0.2:0.2:1.4]); clabel(cs,h); datetick
set(h,'linewidth',1);set(gca,'fontsize',14)
xlabel('Fecha inicial','fontsize',14); ylabel('Rango de gamma en torno al first-guess','fontsize',14)
title('Error RMS','fontsize',16); 
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition", 
     [border, border, (papersize - 2*border)],...
     'paperorientation','landscape')
print('-dpdf','calib_model3.pdf')




% FIGURA 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1,'papersize',[7 7]),clf
h=plot(t0s,min(rmse)','o-',t0s,gammafit,'rx-');
set(h,'linewidth',1.5)
xlabel('Fecha inicial','fontsize',14)
set(gca,'fontsize',14)
datetick
grid on
legend('RMSE(log(I+R))','gamma (1/dia)')
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition", 
     [border, border, (papersize - 2*border)],...
     'paperorientation','landscape')
print('-dpdf','calib_model1.pdf')



% FIGURA 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correr el modelo con los parametros elegidos para la figura
I=casoscorrcfr(find(t==t0)); S=N-I; R=0; X=[S I R];
for i=1:tf-t0+1
   betadum=Re(i)*gamma;
   X=[X;rk4(X(end,:),dt,[betadum gamma N])];
end
S=X(:,1); I=X(:,2); R=X(:,3);


figure(2,'papersize',[9 9]),clf
h=semilogy(t,casos,'b.-',t,casoscorrcfr,'ro-',... 
    t0+[0:length(S)-1]',I+R,'k--');
set(h,'linewidth',1.5)
set(h(end),'linewidth',3)
xlim(datenum(2020,[3 5],[1 1]))
ylim([10 5e4])
xlabel('Dia/mes 2020','fontsize',14)
set(gca,'XTick',datenum(2020,[3 3 4 4 5]',[1 15 1 15 1]'),'fontsize',14)
datetick('dd/mm','keepticks','keeplimits')
legend('Casos confirmados (MINSA)','I+R estimado',...
      ['Modelo SIR (gamma=' num2str(gamma) ' 1/dia)'],'location','southeast')
title('Infectados + Removidos','fontsize',16)
%
papersize = get (gcf, "papersize"); border = 0.25;
set (gcf, "paperposition", 
     [border, border, (papersize - 2*border)],...
     'paperorientation','landscape')
print('-dpdf','calib_model2.pdf')

