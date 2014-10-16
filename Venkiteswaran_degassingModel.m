function [DICconc, d13CDIC, gemass, CO2conc, d13CCO2, pH, Alkconc, volume, geiso] = degassingModel(DICconcinit, d13CDICinit, pHinit, Temperature, varargin)
%% DIC-Alk-CO2 degassing model (Venkiteswaran et al. 2014 PLoS ONE 9(7): e101756, doi: 10.1371/journal.pone.0101756)
% starting with DIC and pH and calculating from there
%
% [DICconc, d13CDIC, gemass, CO2conc, d13CCO2, pH, Alkconc] = ...
%     degassingModel(DICconcinit, d13CDICinit, pHinit, Temperature, varargin)
%     varargin = GW, k, runlength
% GW is GW added in at each timestep in L/timestep

% Handle varargin -- not great code but it works
runlength = 1000;
k = 0.02;
GW = 0;
if length(varargin) == 3
    runlength = varargin{3};
    k = varargin{2};
    GW = varargin{1};
elseif length(varargin) == 2
    k = varargin{2};
    GW = varargin{1};
elseif length(varargin) == 1
    GW = varargin{1};
end;
% disp([length(varargin) size(varargin,2)]);
% disp([ varargin runlength k GW]);

% Go!
TinC = Temperature; % Assigned by input
Rpdb = 0.0112372;
Delta = @(x) (x./Rpdb-1).*1000;
Ratio = @(x) (x./1000+1)*Rpdb;
TinK = TinC + 273.15;

% enrichment factor for CO2(d) to HCO3- Mook et al. 1974
% eb = -(9.484 +/- 0.22) .* 1E3./TinK + (23.89 +/- 0.75);
eb = -(9.484) .* 1E3./TinK + 23.89;
ab = eb./1000 + 1;

% enrichment factor for CO2(g) to CO2(d) Mook et al. 1972 refit of Vogel
% et al. 1970
% eg = -(0.373 +/- 0.07) .* 1E3./TinK + (0.19 +/- 0.23);
eg = -(0.373) .* 1E3./TinK + 0.19;
ag = eg./1000 + 1;

% enrichment factor for HCO3cl- to CO3-- (have to do by difference since
% Zhang et al. 1995 report all epsilon factors relative to CO2(g)
ec = (-0.114 .* TinC + 10.78)  - (-(0.052) .* TinC + 7.22);
%   eHCO3-g   -   eeCO3-g   notation from Zhang
ac = ec./1000 + 1;

% Ka values for CO2(d) and HCO3- Harned and Davis Jr 1943
logKa = -3404.71./TinK + 14.8435 - 0.032786.*TinK;

% Kb values for HCO3- and CO3-- Harned and Scholes Jr 1941
logKb = -2902.39./TinK + 6.4980 - 0.02379.*TinK;

% Henry's Constant for CO2 Mook et al. 1974 (mol/L/atm)
logKH = -(-2385.73./TinK + 14.0184 - 0.0152642.*TinK); % for 0 molality NaCl; extra - since paper gave -logKH

% CO2 equilibrium concentration (umol/L)
pCO2atm = 380; % ppmv
CO2satconc = 10.^logKH * pCO2atm;

% Ionization constant of H2O Dickson and Millero 1987 and Millero 1995
pKw = -log10(exp(148.9802 - 13847.26/TinK - 23.6521*log(TinK))); % assuming zero ionic strength and zero salinity

% model setup
% applied above via nargin
% runlength = 2500; % want this just long enough to get past all the field data
% k = 0.001; % smaller is more accurate, to a point

% preallocate matrices
prealloc = zeros(1,runlength);
gemass = prealloc;
geiso = prealloc;
volume = prealloc;
CO2mass = prealloc;
CO2conc = prealloc;
CO2iso = prealloc;
CO2ratio = prealloc;
HCO3mass = prealloc;
HCO3conc = prealloc;
HCO3ratio = prealloc;
HCO3iso = prealloc;
DICmass = prealloc;
DICconc = prealloc;
DICratio = prealloc;
DICiso = prealloc;
CO3mass = prealloc;
CO3conc = prealloc;
CO3ratio = prealloc;
CO3iso = prealloc;
Hmass = prealloc;
Hconc = prealloc;
pH = prealloc;
OHmass = prealloc;
OHconc = prealloc;
pOH = prealloc;
Alk = prealloc;
Alkmass = prealloc;

% Staring values (start with DIC)
% set initial volume here
volume(1) = 0.1; % 1 L
DICconc(1) = DICconcinit; % Assigned by input
DICmass(1) = DICconc(1)*volume(1);
pH(1) = pHinit; % Assigned by input
Hconc(1) = 10^-pH(1);
Hmass = Hconc(1)*volume(1);
pOH(1) = pKw - pH(1);
OHconc(1) = 10^-pOH(1);
OHmass = OHconc(1)*volume(1);
CO2conc(1) = DICconc(1) .* ((10.^-pH(1)).^2) ./ ((10.^-pH(1)).^2 + 10.^logKa .* 10.^-pH(1) + 10.^logKa.*10.^logKb);
% CO2conc(1) = DICconc(1) .* (1 + 10.^logKa./10.^-pH + 10.^logKb.*10.^logKa./(10.^-pH).^2).^-1;
CO2mass(1) = CO2conc(1)*volume(1);
HCO3conc(1) = DICconc(1) .* (10^logKa .* 10^-pH(1)) ./ ((10.^-pH(1)).^2 + 10^logKa .* 10^-pH(1) + 10^logKa .* 10^logKb);
% HCO3conc(1) = DICconc(1) .* (1 + 10.^-pH./10^logKa + 10^logKb./10^-pH).^-1;
HCO3mass(1) = HCO3conc(1)*volume(1);
CO3conc(1) = DICconc(1) .* (10^logKa .* 10^logKb) ./ ((10.^-pH(1)).^2 + 10^logKa .* 10^-pH(1) + 10^logKa .* 10^logKb);
% CO3conc(1) = DICconc(1) .* (1 + (10.^-pH(1)).^2./(10.^logKb.*10.^logKa) + 10.^-pH(1)./10.^logKb).^-1;
Alkconc(1) = HCO3conc(1) + 2*CO3conc(1) + 10^-pOH(1) - 10^-pH(1); % has to be constant in this model IF there is no GW addition
Alkmass(1) = Alkconc(1)*volume(1);
DICratio(1) = Ratio(d13CDICinit); % Assigned by input
DICiso(1) = DICmass(1) .* DICratio(1);
CO2ratio(1) = (DICconc(1)*ab*ac*DICratio(1)) / (CO3conc(1) + HCO3conc(1)*ac + CO2conc(1)*ab*ac);
CO2iso(1) = CO2mass(1) .* CO2ratio(1);
HCO3ratio(1) = CO2ratio(1)./ab;
HCO3iso(1) = HCO3mass(1) .* HCO3ratio(1);
CO3ratio(1) = HCO3ratio(1)./ac;
CO3iso(1) = CO3mass(1) .* CO3ratio(1);


% Model itself
for i = 2:runlength
    gemass(i) = k.* (CO2satconc - CO2conc(i-1)); % length/time * mass/length^3 = mass/length^2/time; assume area = 1 and thus mass/time
    geiso(i) = k .* 0.9987 .* (CO2satconc .* Ratio(-7.75) .* 0.9989   -   CO2conc(i-1) .* CO2ratio(i-1));
    % Assume:
    % 1. CO2 loss = DIC loss (both mass and isotopes)
    % 2. HCO3 and CO3 have to decline and pH has to go up to keep Alk
    % constant
    % 3. Use constant Alk with new DIC to get pH from pH_from_Alk model
    % 4. Use new pH to calculate CO2, HCO3, CO3
    volume(i) = volume(i-1)+GW;
    DICmass(i) = DICmass(i-1) + gemass(i) + DICconc(1)*GW;
    DICconc(i) = DICmass(i)/volume(i);
    DICiso(i) = DICiso(i-1) + geiso(i) + DICiso(1)/volume(1)*GW; % the GW portion has to be mol/L*L to get a mass, so use iso/vol, i.e. mass*ratio/volume
    DICratio(i) = DICiso(i)/DICmass(i);

    [pH(i) pOH(i)] = pH_from_Alk(DICconc(i), Alkconc(1), TinC); % pH accounting for DIC & Alk from GW but not the H from GW
    Hconc(i) = 10^-pH(i); % pH accounting for DIC & Alk then add the H from GW separately
    OHconc(i) = 10.^-pOH(i);
    Hmass(i) = Hconc(i)*volume(i) + Hconc(1)*GW; % Add the H from GW and then recalculate the concentration
    OHmass(i) = OHconc(i)*volume(i) + OHconc(1)*GW;
    Hconc(i) = Hmass(i)/volume(i); % recalc with the new H from GW
    OHconc(i) = OHmass(i)/volume(i);
    pH(i) = -log10(Hconc(i)); % "new" pH is thus affected by the constand Alk requirement and by new H from GW
    pOH(i) = -log10(OHconc(i));
    pOH(i) = pKw - pH(i); % assuming the pOH is not affected by new OH from GW to close the pKw:pH:pOH loop

    CO2conc(i) = DICconc(i) .* ((10.^-pH(i)).^2) ./ ((10.^-pH(i)).^2 + 10.^logKa .* 10.^-pH(i) + 10.^logKa.*10.^logKb);
    HCO3conc(i) = DICconc(i) .* (10^logKa .* 10^-pH(i)) ./ ((10.^-pH(i)).^2 + 10^logKa .* 10^-pH(i) + 10^logKa .* 10^logKb);
    CO3conc(i) = DICconc(i) - CO2conc(i) - HCO3conc(i); % forced by DIC constant
    % don't need the masses of CO2, HCO3, or CO3
    CO2ratio(i) = (DICconc(i)*ab*ac*DICratio(i)) / (CO3conc(i) + HCO3conc(i)*ac + CO2conc(i)*ab*ac);
    CO2iso(i) = CO2mass(i) .* CO2ratio(i);
    HCO3ratio(i) = CO2ratio(i)./ab;
    HCO3iso(i) = HCO3mass(i) .* HCO3ratio(i);
    CO3ratio(i) = HCO3ratio(i)./ac;
    CO3iso(i) = CO3mass(i) .* CO3ratio(i);
end

Alkconc = HCO3conc + 2*CO3conc + 10.^-pOH - 10.^-pH; % In case we need it
geratio = geiso./gemass; % to make it easier to plot, etc.
d13CDIC = Delta(DICratio);  % For output
d13CCO2 = Delta(CO2ratio);  % For output

end
%EOF
