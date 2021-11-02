function [E, e, kf] = bartonBandis(stress, jrc, jcs, varargin)
% 
% SYNOPSIS:
%   [E, e, kf] = bartonBandis(stress, varargin)
%
% DESCRIPTION: Relate stress to aperture and perms empirically. Empirical relationship based on laboratory experiments, relates stress to 
%              hydraulic aperture and permeability. Note: may need to put ranges of validity as an assetion.
%    
% PARAMETERS:
%   stress - applied outer stress vector, Pa (will be converted to MPa in the code)
%   jrc    - joint roughness coefficient, dimensionless
%   jcs    - joint compressive strength, MPa
%   'mode' - 'single' calculate a single kf based on one mean Vm
%          - 'range' calculate range of kf based on range of A,B,C,D stddev  
% 
%   'sigmac'          - (optional) unconfined compression strength (rock adjacent to the wall)   
%                       if not given, then equaql to jcs
%   initial_stiffness - (optional), MPa/mm 
%   initial_aperture  - (optional), mm
%  
%
% RETURNS:
%   E  - Normal mechanical aperture 
%   e  - Hydraulic aperture
%   kf - Effective perm, based on || plates 
%
% EXAMPLE:
%
% SEE ALSO:
%
% Author: Aidan Kubeyev, initially Maier C.

%{

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


    %% Inputs 
    opt = struct('initial_stiffness', 0, ...
                 'initial_aperture' , 0, ...
                 'sigmac'           , jcs, ...      
                 'mode'             , 'single');
    
             
    opt = merge_options(opt, varargin{:});   
    
    % Converting normal stress to MPa
    sigman = stress./(mega*Pascal);
    
    % Formula is for specific range. Caution if range is outside
    if jrc < 5  ||  jrc > 15
        warning('JRC must be in range 5-15'); 
    end
       
    if jcs < 22  ||  jcs > 182
        warning('JCS ideally must be in range 22-182 Mpa'); 
    end
    
       
       
    %% Initial mechanical aperture, mm
    
       
    if(~opt.initial_aperture)
        E0 = (jrc/5)*(0.2*opt.sigmac/jcs - 0.1);
    else
        E0 = opt.initial_aperture;
    end
    
   if E0 < 0.099  ||  E0 > 0.601 % 0.1 and 0.6
        warning('Initial aperture ideally must be in range 0.1-0.6 mm'); 
   end
       
    
    %% Maximum closure, mm
    
    switch opt.mode
        
        case('range')
    
        %----------------------------------------------------------------------
        % Maximum closure

        M = [-0.1032  -0.1712  -0.0352 ;...    % A Amin Amax
             -0.0074  -0.0113  -0.0035 ;...    % B Bmin Bmax
              1.135    0.8089   1.4611 ;...    % C Cmin Cmax
             -0.251   -0.3539  -0.1481    ];   % D Dmin Dmax 


        vm = [];     
        for i = 1 : length(M(1, :))      
            for j = 1 : length(M(2, :))
                for k = 1 : length(M(3, :))
                   for l = 1 : length(M(4, :))

                       vmn = M(1, i) + M(2, j)*jrc + M(3, k)*(jcs/E0)^M(4, l); 
                       vm  = [vm; vmn]; 

                   end
                end
            end
        end

        % Delete negatives and which are vm <= E0 
        vm(vm < 0) = []; 
        vm(vm > E0) = [];
     
            
        case('single')    
        
             %      A         B             C                D   
             vm = -0.1032 - 0.0074.*jrc + 1.135.*(jcs/E0)^(-0.251);
                   
    end
    
    assert(all(vm >= 0),  'Maximum closure cannot be less than 0');
    assert(all(vm <= E0), 'Maximum closure must be less than initial aperture.');
    assert(~isempty(vm),  'Maximum closure is an empty vector'); 
    
    
    
    %% Initial Stiffness, MPa/mm
    
    if(~opt.initial_stiffness) 
        Kni = -7.15 + 1.75.*jrc + 0.02.*jcs/E0;
    else
        Kni = opt.initial_stiffness;
    end
    
    
    %% Normal mechanical aperture, mm 
    
    % For each maximum closure calculate apertures at various stress,(in mm) 
    for i = 1 : length(vm)
        
        % E n-to-m matrix, where n is results for  various vm, n - stress
        deltaVj(i, :) =  ( 1./vm(i) + Kni./sigman ).^-1;

    end
    
    
    E = E0 - deltaVj; % Both are in mm
    
    
    % Plot joint closure
    %figure(11);
    %plot(deltaVj*1000, sigman); % micrometre, MPa
    %ylabel('Stress, Mpa');
    %xlabel('deltaVj, micrometer'); 
    
    
    % Convert to micro meters. 
    Emicro = E.*milli/micro;

    
    %% Hydraulic aperture (in micro meters)
    
    emicro = (Emicro.^2)./(jrc.^2.5);
    

    %% Converting to S.I. and calculate permeability
   
    e = emicro.*micro;
    E = E.*milli;
    
    kf = (e.^2) ./ 12; % m2
    
    
    
    % Plot perm in 'log10 cm2' as in BartonBandis paper
    % kf1 = log10(kf*10000);    % cm2
    % figure(12); 
    % plot(sigman, kf1); % micrometre, MPa
    % xlabel('Stress, Mpa');
    % ylabel('kf, log cm2'); 
    
   
end