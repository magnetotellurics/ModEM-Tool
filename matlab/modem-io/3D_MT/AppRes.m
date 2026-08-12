%************************************************************************
%  AppRes : Calculates the Apparent Resistivity for given impedances
% USAGE:  [AppZxx, AppZxy, AppZyx, AppZyy] = AppRes(Per, Zxx, Zxy, Zyx, Zyy);
% Nx = grid dimensions in X-axis 
% Ny = grid dimensions in Y-axis 
% NPer = Number of Periods
% Per = peiods in seconds
% Zxx(Nx, Ny, NPer) = (complex) impedance values in xx for given NPer
% Zxx(Nx, Ny, NPer) = (complex) impedance values in xy for given NPer
% Zxx(Nx, Ny, NPer) = (complex) impedance values in yx for given NPer
% Zxx(Nx, Ny, NPer) = (complex) impedance values in yy for given NPer
% Nx and Ny are the values for center of the cube on the surface (Nz = 0)
% Note: Input impedances are in SI Units
% Author: Kush Tandon, OSU, 06/19/03
function [AppZxx, AppZxy, AppZyx, AppZyy] = AppRes(Per, Zxx, Zxy, Zyx, Zyy)

% Nlines keeps track of number of lines of dataset that was read
Nlines = 0;
[Nx, Ny, NPer] = size(Zxx);
for ip=1:NPer
        % read the impedances
        % rows denote the X-direction
        % colums denote the Y-direction
        for j=1:Nx
            for k = 1:Ny 
                % the SI values are to be converted to (mV/ km)/ nT
                % multiplied by a factor of 1000 for apparent resistivity
                % formula used
                Zxx(j, k, ip) = 1E3*Zxx(j, k, ip);
                Zxy(j, k, ip) = 1E3*Zxy(j, k, ip);
                Zyx(j, k, ip) = 1E3*Zyx(j, k, ip);
                Zyy(j, k, ip) = 1E3*Zyy(j, k, ip);
                AbsZxx = abs(Zxx(j,k,ip));
                AbsZxy = abs(Zxy(j,k,ip));
                AbsZyx = abs(Zyx(j,k,ip));
                AbsZyy = abs(Zyy(j,k,ip));
                AppZxx(j, k, ip) = (Per(ip)*(AbsZxx.^2))/5.;
                AppZxy(j, k, ip) = (Per(ip)*(AbsZxy.^2))/5.;
                AppZyx(j, k, ip) = (Per(ip)*(AbsZyx.^2))/5.;
                AppZyy(j, k, ip) = (Per(ip)*(AbsZyy.^2))/5.;
                Nlines = Nlines + 1;
            end %k
         end %j
    end %ip
end
