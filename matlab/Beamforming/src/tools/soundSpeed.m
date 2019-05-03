function [c0, trusty] = soundSpeed( t, h, p, xc_ppm )
%{
    Estimate the sound of speed in the air provided some weather
    conditions.
    
    INPUT:
    t - Temperature in Celsius degrees
    h - Relative humidity (%)
    p - Pressure (Pa)
    xc_ppm - CO2 concentration (ppm)
    
    OUTPUT:
    c0 - Sound of speed (m/s)
    trusty - Trusty model flag
                                            luis.gato@copsonic.com
%}

% Initial checks
  narginchk(0,4);
  nargoutchk(0,2);

  switch nargin
      case 0
          t = 20;         % Default temperature
          h = 0;          % Dry air
          p = 101325;     % 1 atm
          xc_ppm = 350;   % Typical value
      case 1
          h = 0;          % Dry air
          p = 101325;     % 1 atm
          xc_ppm = 350;   % Typical value
      case 2
          p = 101325;     % 1 atm
          xc_ppm = 350;   % Typical value
      case 3
          xc_ppm = 350;   % Typical value
  end

  if t <= -273.15 || p <= 0 || h < 0 || h > 100 || xc_ppm < 0
      disp(' ')
      error(' Invalid input parameters. Please check.')
  end

% Weather parameters:
  T = t + 273.15;     % Temperature in Kelvin degrees
  xc = xc_ppm/1e6;
  
% Coefficients:
  a0 = 331.5024;
  a1 = 0.603055;
  a2 = -0.000528;
  a3 = 51.471935;
  a4 = 0.1495874;
  a5 = -0.000782;
  a6 = -1.82e-7;
  a7 = 3.73e-8;
  a8 = -2.93e-10;
  a9 = -85.20931;
  a10 = -0.228525;
  a11 = 5.91e-5;
  a12 = -2.835149;
  a13 = -2.15e-13;
  a14 = 29.179762;
  a15 = 0.000486;

% Main contribution of temperature:
  ct = a0 + a1*t + a2*t^2;
  
% Main contribution of water vapor:  
  xi = 1.00062 + 3.14e-8*p + 5.6e-7*t^2;
  psv = exp(1.2811805e-5*T^2 - 1.9509874e-2*T + 34.04926034 -  6.3536311e3/T);
  xw = h*xi*psv/p/100;
  cw = (a3 + a4*t + a5*t^2)*xw;
  
% Main contribution of pressure:
  cp = (a6 + a7*t + a8*t^2)*p;
  
% Main contribution of carbon dioxide concentration:
  cc = (a9 + a10*t + a11*t^2)*xc;
  
% Crossed terms:
  cm = a12*xw^2 + a13*p^2 + a14*xc^2 + a15*xw*p*xc;
  
% Estimated sound speed:
  c0 = ct + cw + cp + cc + cm;
  
  if t < 0 || t > 30 || p < 75e3 || p > 102e3 || xc_ppm > 10e3 || xw > 0.06
      disp(' ')
      warning(' The mathematical approximation is not reliable with the provided parameters.')
      trusty = 0;
  else
      trusty = 1;
  end

end

