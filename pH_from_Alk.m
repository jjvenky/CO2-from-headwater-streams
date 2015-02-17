function [modelpH, modelpOH] = pH_from_Alk(DIC, Alk, TinC)

% give DIC, Alk and Temp
% get pH and pOH using fminsearch
%
% rearrangement of [Alk] = [HCO3-] + 2[CO3--] + [OH-] - [H+]

TinK = TinC + 273.15;

% Ka values for CO2(d) and HCO3- Harned and Davis Jr 1943
logKa = -3404.71./TinK + 14.8435 - 0.032786.*TinK;

% Kb values for HCO3- and CO3-- Harned and Scholes Jr 1941
logKb = -2902.39./TinK + 6.4980 - 0.02379.*TinK;

% Ionozation constant of H2O Millero et al. 1987 and Millero 1995
pKw = -log10(exp(148.9802 - 13847.26/TinK - 23.6521*log(TinK))); % assuming zero ionic strength and zero salinity

Ka = 10^logKa;
Kb = 10^logKb;


% Call fminsearch with a random starting point.
% start_point = rand(1)*pKw; % random point within 0:pKw
% now trying 6 since that's pretty close
model = @expfun;
% options = optimset('fminsearch');   % Use FMINSEARCH defaults
% [modelpH, FVAL, EXITFLAG, OPTIONS] = fminsearch(model, start_point, options);
[modelpH, modelpOH] = fminsearch(model, 6); % no options to speed up
% disp(modelpH);
% disp(FVAL);
% disp(EXITFLAG);
% disp(OPTIONS);

function [sse, modelAlk, modelpH] = expfun(params)
    pH = params; % want to solve for this value
    pOH = pKw - pH; % need this for calculations
    H = 10^-pH; % need this for calculations
    OH = 10^-pOH; % need this for calculations
    modelAlk = DIC * (Ka*H) / (H^2 + Ka*H + Ka*Kb) + ...
        2*DIC * (Ka*Kb) / (H^2 + Ka*H + Ka*Kb) + ...
        OH - H;
    ErrorVector = modelAlk - Alk;
    sse = ErrorVector ^ 2;
    modelpH = pH;
    modelpOH = pOH;
end

end
