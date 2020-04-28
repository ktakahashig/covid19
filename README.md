# covid19
Codigo fuente (octave/matlab) y datos usados para el modelado SIR

deriv.m : el lado derecho de las ecuaciones del modelo SIR
rk4.m : método Runge-Kutta de cuarto orden

calib_model.m : Permite estimar el mejor valor de gamma.
          Input    casos_confirmados.txt : Datos de casos confirmados (MINSA)
                   fraction_reported.txt : % casos reportados estimado (CMMID; digitalizado aprox)
                   Re.txt : Número de reproducción en tiempo real (CMMID)
          Output   calib_model1.pdf : Figura de error RMS y mejor gamma por fecha de inicio
                   calib_model2.pdf : Figura de I+R estimado vs simulación
                   calib_model3.pdf : Figura de error para verificar estimación de gamma
                   
seas_wang.m : Estima el ciclo anual de Re según las fórmulas 1 y 2 de Wang et al (2020) y datos meteorológicos
          Input    campo_de_marte.txt : Datos horarios (2017-2020) de la estación Campo de Marte, Lima
          Output   campo_de_marte.pdf : Figura de temperatura, humedad relativa y humedad absoluta, con ciclo anual
                   R_wang_et_al.pdf   : Figura de Re y su ciclo anual con ambas formulas
                   
                   
