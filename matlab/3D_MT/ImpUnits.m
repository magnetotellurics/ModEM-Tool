function SI_factor = ImpUnits(oldUnits,newUnits)

% function SI_factor = ImpUnits(oldUnits,newUnits)
%
% Computes the value by which the data must be multiplied to convert
% from the old units to the new units. 
% The units may be any of the following.
% 1) SI units for E/B: [V/m]/[T] (used in ModEM code)
% 2) practical units for E/B: [mV/km]/[nT]
% 3) SI units for E/H: [V/m]/[A/m] = Ohm
%  (c) Anna Kelbert, 2009

% if the units are the same, nothing to convert, exit
if strcmp(oldUnits,newUnits)
    SI_factor = 1;
    return;
elseif strcmp(oldUnits,'[]') && strcmp(newUnits,'[mV/km]/[nT]')
    disp('Skipped unit conversion from apparent resistivity & phase to impedance');
    SI_factor = 1;
    return;
elseif strcmp(oldUnits,'[mV/km]/[nT]') && strcmp(newUnits,'[]')
    disp('Skipped unit conversion from impedance to apparent resistivity & phase');
    SI_factor = 1;
    return;
end

% first convert the old units to [V/m]/[T]
if ~isempty(strfind(oldUnits,'[V/m]/[T]'))
    factor1 = 1.0;
elseif ~isempty(strfind(oldUnits,'[mV/km]/[nT]'))
    factor1 = 1000.0;
elseif ~isempty(strfind(oldUnits,'[V/m]/[A/m]')) ...
        || ~isempty(strfind(oldUnits,'Ohm'))
    factor1 = 1000.0 * 10000.0/(4*pi); % approx. 796000.0
else
    warning(['Unknown old units ' oldUnits ': unit conversion failed']);
    factor1 = 1;
end

% now convert [V/m]/[T] to the new units
if ~isempty(strfind(newUnits,'[V/m]/[T]'))
    factor2 = 1.0;
elseif ~isempty(strfind(newUnits,'[mV/km]/[nT]'))
    factor2 = 1.0 / (1000.0);
elseif ~isempty(strfind(newUnits,'[V/m]/[A/m]')) ...
        || ~isempty(strfind(newUnits,'Ohm'))
    factor2 = 1.0 / (1000.0 * 10000.0/(4*pi));
else
    warning(['Unknown new units ' newUnits ': unit conversion failed']);
    factor2 = 1;
end   

SI_factor = factor1 * factor2;