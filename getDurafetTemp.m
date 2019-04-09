function tempC = getDurafetTemp(VTherm,TCOffset)
% Convert Durafet thermistor voltage to temperature (C) using following polynomial
c0 = 340.9819863; c1 = -9.10257E-05; c2 = -95.08806667; c3 = 0.965370274;
RTherm = 20000./(3.3./VTherm-1); 
tempC = c0+c1*RTherm+c2*log10(RTherm)+c3*(log10(RTherm)).^3;
tempC = tempC+TCOffset;
