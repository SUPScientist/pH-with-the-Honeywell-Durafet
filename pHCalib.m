function calib = pHCalib(calEint,calEext,calpH,calT,calsal)
% Univ gas constant, Faraday constant,
R = 8.3145; F = 96487;
% Temperature dependence of standard potentials, Martz et al. 2010
dE0int = -0.00125; dE0ext = -0.001048;
% See Martz et al. 2010 for greater detail
tempK = calT+273.15; % Convert temp from C to K
S_T = (R*tempK)/F*log(10); % Nernst temp dependence
E0int = calEint-S_T*calpH; % Calc E0int from Nernst & pH @ calibration point
E0int25 = E0int+dE0int*(25-calT);
Z = 19.924.*calsal./(1000-1.005.*calsal); % Ionic strength, Dickson et al. 2007
SO4_tot = (0.14/96.062)*(calsal./1.80655);  % Total conservative sulfate
cCl = 0.99889/35.453.*calsal/1.80655; % Conservative chloride
mCl = cCl*1000/(1000-calsal.*35.165/35); % mol/kg-H2O
K_HSO4 = exp(-4276.1/tempK+141.328-23.093*log(tempK)...
          +(-13856/tempK+324.57-47.986*log(tempK))*Z^0.5...
          +(35474/tempK-771.54+114.723*log(tempK))*Z-2698/tempK*Z^1.5...
          +1776/tempK*Z^2+log(1-0.001005.*calsal)); % Bisulfate equilibrium const., Dickson et al. 2007
pHint_free = calpH+log10(1+SO4_tot/K_HSO4);
cHfree = 10^(-pHint_free); % mol/kg-sw
pHint_free = pHint_free+log10((1000-calsal.*35.165/35)/1000); % mol/kg-H2O
mHfree = 10^(-pHint_free); % mol/kg-H2O
DHconst = 0.00000343*calT^2+0.00067524*calT+0.49172143; % Debye-Huckel, Khoo et al. 1977
log10gamma_HCl = 2*(-DHconst*sqrt(Z)/(1+1.394*sqrt(Z))+(0.08885-0.000111*calT)*Z);
aHfree_aCl = mHfree*mCl*10^(log10gamma_HCl);
E0ext = calEext+S_T*log10(aHfree_aCl);
E0ext25 = E0ext+dE0ext*(25-calT);

calib = [E0int25 E0ext25];
