function SI_factor = ImpUnits(oldUnits,newUnits)
% Computes the value by which the data must be multiplied to convert
%
%   Units may be:
%   * SI units for E/B: '[V/m]/[T]' (used in ModEM code)
%   * Practical units for E/B: '[mV/km]/[nT]'
%   * SI units for E/H: '[V/m]/[A/m]' or 'Ohm'
%   * Apparent Resistivity: '[]'
%   
% Usage:
%   factor = ImpUnits("[V/m]/[T]", "[V/m]/[A/m]")
%   factor = ImpUnits("[]", "[V/m]/[T]")
% Input Arguments:
%   oldUnits - string, valid units listed above
%   newUnits - string, valid unit listed above
%
%  (c) Anna Kelbert, 2009

APPARENT_RES = "[]";
SI_UNITS_EB = "[V/m]/[T]";
PRACTICAL_UNITS_EB = "[mV/km]/[nT]";
SI_UNITS_EH = "[V/m]/[A/m]";
OHMS = "Ohm";

VALID_UNITS = [APPARENT_RES SI_UNITS_EB PRACTICAL_UNITS_EB SI_UNITS_EH OHMS];

if ~ismember(oldUnits, VALID_UNITS)
    warning("WARNING: Unkown oldunits %s conversion failed (returning 1)", oldUnits)
    warning("WARNING: Valid units are: %s", VALID_UNITS.join(', '))
    SI_factor = 1;
    return 
end 

if ~ismember(newUnits, VALID_UNITS)
    warning("WARNING: Unkown newUnits %s conversion failed (returning 1)", newUnits)
    warning("WARNING: Valid units are: %s", VALID_UNITS.join(', '))
    SI_factor = 1;
    return
end 

 
% if the units are the same, nothing to convert, exit
if strcmp(oldUnits, newUnits)
    SI_factor = 1;
    return;
elseif strcmp(oldUnits, APPARENT_RES) && strcmp(newUnits, PRACTICAL_UNITS_EB)
    disp('Skipped unit conversion from apparent resistivity & phase to impedance');
    SI_factor = 1;
    return;
elseif strcmp(oldUnits, PRACTICAL_UNITS_EB) && strcmp(newUnits, APPARENT_RES)
    disp('Skipped unit conversion from impedance to apparent resistivity & phase');
    SI_factor = 1;
    return;
end

% first convert the old units to [V/m]/[T]
if contains(oldUnits, SI_UNITS_EB)
    factor1 = 1.0;
elseif contains(oldUnits, PRACTICAL_UNITS_EB)
    factor1 = 1000.0;
elseif contains(oldUnits,SI_UNITS_EH) || contains(oldUnits,OHMS)
    factor1 = 1000.0 * 10000.0/(4*pi); % approx. 796000.0
end

% now convert [V/m]/[T] to the new units
if contains(newUnits, SI_UNITS_EB)
    factor2 = 1.0;
elseif contains(newUnits, PRACTICAL_UNITS_EB)
    factor2 = 1.0 / (1000.0);
elseif contains(newUnits, SI_UNITS_EH) || contains(newUnits, OHMS)
    factor2 = 1.0 / (1000.0 * 10000.0/(4*pi));
end   

SI_factor = factor1 * factor2;