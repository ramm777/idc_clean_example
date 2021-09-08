function [perm, permf, finfo] = permStokes(G, varargin)
% 
% SYNOPSIS:
%   perm = permStokes(G, varargin); 
%
% DESCRIPTION: 
%   Calculate permeability in a random 2D section in a perpendicular to the 2D plane domain based on Stokes equation. 
%
% PARAMETERS:
%   G                - Grid desribing area in which permeability is calculated
%   matrixvol        - Volume (area here) of the RHSbody and LHSbody (i.e. matrix)
%   mode             - 'loopEachCell' or 'largePolygon'
%                    - if 'polygon' is selected, LHSbody and RHS body must be provided  
%   'getFlowInfo'    - if true, outputs additional flow info such as u(vel).
%
% RETURNS:
%   perm             - Permeability inc. matrix in G, in perpendicular direction in mD units
%   permf            - Permeability fracture only  in G, in perpendicular direction in mD units
%   finfo            - Flow information structure: velocity, corresponding grid G
%
% EXAMPLE:
%
% SEE ALSO:
% 
% AUTHOR:
%   Aman Kubeyev  

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


   
    opt = struct('matrixvol',              0, ...
                 'mode',      'loopEachCell', ...
                 'pcoords',               [], ...
                 'getFlowInfo',         false, ...
                 'nx', 500                      );


    opt = merge_options(opt, varargin{:});


    finfo.vel = []; 
    finfo.G   = [];
 
    

%% CONTROL

    disp('CALCULATING: Stokes perm...');
    mrstModule add incomp

%% Create Box Grid which superimposes all the input grid (G) geometry

   
    coords = G.nodes.coords; 

    
    a = abs((min(coords(:, 2)) - max(coords(:, 2))) / (min(coords(:, 1)) - max(coords(:, 1)))); % Aspect ratio of the Grid
    nx = opt.nx;                   
    ny = nx*ceil(a);

    
    x = linspace(min(coords(:, 1)), max(coords(:, 1)), nx+1);
    y = linspace(min(coords(:, 2)), max(coords(:, 2)), ny+1); 


    Gbox = tensorGrid(x, y);
    Gbox = computeGeometry(Gbox); 
    bccoords = Gbox.cells.centroids; 


%%  Select which cells lying outer of Gfrac (fracture)


switch opt.mode
    
    case ('loopEachCell')
    
        
        %--------------------------------------------------------------------------
        % Data preparation

        % Matrix of demcells coords for same format as in inpolygon()
        x = zeros(G.cells.num, 4); 
        y = zeros(G.cells.num, 4); 

        for z = 1:G.cells.num
            % For z cell extract coords
            coord = cellCoords(G, z);
            % Place x and y coords as 2 separate matrices 
            x(z, :) = coord(:, 1)';
            y(z, :) = coord(:, 2)'; 

        end

        % Sort for function inpolygon() to work, otherwise it is not a polygon
        sx = [x(:, 1), x(:, 2), x(:, 4), x(:, 3)]; % Sorted x
        sy = [y(:, 1), y(:, 2), y(:, 4), y(:, 3)]; % Sorted y 

        clear coord


        %--------------------------------------------------------------------------
        % Loop over all cells of ellipse to check which cells_centres are inside     

        k = 0; 
        ncells = []; 

        for i = 1 : G.cells.num 
            k = k + 1;

            
            in = inpolygon(bccoords(:, 1), bccoords(:, 2), sx(k, :), sy(k, :)); % Find cells' nodes incontact including current demcell (self-contact)
            tlcells = find(in);                                                 % To leave cells among the box cells
            ncells = [ncells; tlcells];                                         % Accumulate new coords


        end
        

    case('polygon')
        
        warning('IDC: potentially, you can check for the polygon if its self-intersecting for QC')
        
        in = inpolygon(bccoords(:, 1), bccoords(:, 2), opt.pcoords(:, 1), opt.pcoords(:, 2)); % Find cells' nodes incontact including current demcell (self-contact)
        ncells = find(in);           % To leave cells among the box cells
                   
end
        

%% Extract Subgrid


G = extractSubgrid(Gbox, ncells);
G = computeGeometry(G);


% Helper plot
% figure(1)
% plot(ncoords(:, 1), ncoords(:, 2), 'x');    
% axis equal; 
% hold on
% plot(coords(:, 1), coords(:, 2), 'o'); 


%% FLOW Parameters / gravity


perm = ones(G.cells.num, 1); 
poro = ones(G.cells.num, 1);


rock = struct('perm', reshape(perm, [], 1),     ...
              'poro', reshape(poro, [], 1), ...
              'alpha', ones(G.cells.num, 1)); 

mu = 0.001;           
fluid = initSingleFluid('mu', mu,'rho', 1000);


pressure = zeros(G.cells.num, 1);
pressure = reshape(pressure, [], 1);


state = struct('pressure', pressure, 's', ones(G.cells.num, 1), 'flux', zeros(G.faces.num, 1));


%% BC 


dpdx = 10000;       % in SI, Pa/m
c    = dpdx/(4*mu); % in SI


bfaces = boundaryFaces(G);                   % Finds boundary faces
bcfcentroids = G.faces.centroids(bfaces, :); % Get face centroids for each face


Ubc  = c.*( bcfcentroids(:, 2).^2 + bcfcentroids(:, 1).^2 ); %  Calculate Ubc
bc_p = addBC([], bfaces, 'pressure', Ubc);


%% SOLVE FLOW 

T  = computeTrans(G, rock); 


state = incompTPFA(state, G, T, fluid, 'bc', bc_p);
U = state.pressure; 
u = U - c.*(G.cells.centroids(:, 2).^2 + G.cells.centroids(:, 1).^2);


% figure(5); clf; 
% plotCellData(G, u, 'EdgeAlpha', 0.001); colorbar; axis equal
% title('Velocity, m/s'); 
% xlim([7e-3, 7.45e-3]); ylim([0, 3.5e-4]); % Image for the paper



%% CALCULATE perm

upper = mu.*sum(u.*G.cells.volumes); 

lower  = dpdx.*(sum(G.cells.volumes) + opt.matrixvol);
lowerf = dpdx.*sum(G.cells.volumes);

k  = upper/lower;
kf = upper/lowerf;

perm  = k  / (milli*darcy); % Convert to mD 
permf = kf / (milli*darcy); % Convert to mD 



if  opt.getFlowInfo    
       
    finfo.vel = u; % Flow velocity
    finfo.G   = G; 
        
end




%% PLOT bc-vectors

% Find faces on left/right of an ellipse
% l = bfaces((bcfcentroids(:, 1) < G.faces.centroids(151, 1))); 
% r = bfaces((bcfcentroids(:, 1) > G.faces.centroids(101, 1)));
% 
% % Get face centroids for each face
% lfcentroids = G.faces.centroids(l, :);
% rfcentroids = G.faces.centroids(r, :);
% 
% lUbc  = c.*( lfcentroids(:, 2).^2 + lfcentroids(:, 1).^2 );
% rUbc  = c.*( rfcentroids(:, 2).^2 + rfcentroids(:, 1).^2 );
% 
% figure(1); clf; hold on;
% plotGrid(G); 
% quiver(lfcentroids(:, 1), lfcentroids(:, 2), lUbc, zeros(length(lUbc), 1), 'Color', 'r'); 
% quiver(rfcentroids(:, 1), rfcentroids(:, 2), -rUbc, zeros(length(rUbc), 1), 'Color', 'r'); 
% axis equal 

%% APPENDIX
% Below is based on Poisson equation, not finished

%--------------------------------------------------------------------------
% Refine Grid Triangulation

% % Save original grid
% Gorig = G; 
% 
% coords = G.nodes.coords; 
% f = 20; % factor to make node equal dist on x and y
% 
% % Aspect ratio of the Grid
% a = abs((min(coords(:, 2)) - max(coords(:, 2))) / (min(coords(:, 1)) - max(coords(:, 1))));
% 
% 
% % Box coords encompassing the ellipse. f*a - factor to make node equal dist on x and y
% bx = linspace(min(coords(:, 1)), max(coords(:, 1)), f)';
% by = linspace(min(coords(:, 2)), max(coords(:, 2)), f*a)'; 
% 
% bcoords = [];
% for i = 1 : f
%     bcoords = [bcoords; [repmat(bx(i), length(by), 1), by]];
% end
% 
% % Box nodes distances from each other
% bndist = abs( (min(coords(:, 1)) - max(coords(:, 1)))/ f );
% 
% %--------------------------------------------------------------------------
% % Data preparation
% 
% % Matrix of demcells coords for same format as in inpolygon()
% x = zeros(G.cells.num, 4); % Preallocate
% y = zeros(G.cells.num, 4); % Preallocate
% 
% for z = 1:G.cells.num
%     % For z cell extract coords
%     coord = cellCoords(G, z);
%     % Place x and y coords as 2 separate matrices 
%     x(z, :) = coord(:, 1)';
%     y(z, :) = coord(:, 2)'; 
% 
% end
% 
% % Sort for function inpolygon() to work, otherwise it is not a polygon
% sx = [x(:, 1), x(:, 2), x(:, 4), x(:, 3)]; % Sorted x
% sy = [y(:, 1), y(:, 2), y(:, 4), y(:, 3)]; % Sorted y 
% 
% clear coord
% 
% %--------------------------------------------------------------------------
% % Loop over all cells of ellipse to check which nodes in coords are inside     
% 
% k = 0; 
% ncoords = []; 
% 
% for i = 1 : G.cells.num 
%     k = k + 1;
%    
%     % Find cells' nodes incontact including current demcell (self-contact)
%     in = inpolygon(bcoords(:, 1), bcoords(:, 2), sx(k, :), sy(k, :));
%     
%     % To leave coords among the box coords
%     tlcoords = bcoords(in, :);
%     
%     % Accumulate new coords
%     ncoords = [ncoords; tlcoords]; 
%   
%     
% end
% 
% %--------------------------------------------------------------------------
% % Delete nodes in ncoords that are too close to coords. Intracells are ok. 
% 
% % Calculate distances 
% for i = 1 : numel(coords(:, 1))
% 
%     distances = sqrt(sum(bsxfun(@minus, coords(i, :), ncoords).^2,2));
%     
%     if any( distances < (bndist/100) )
%         
%         % Find index of coords that is too close to ncoords
%         ind = find( distances < (bndist/100) );
%         
%         % Delete which is too close
%         ncoords(ind, :) = [];  
%         
%     end
%     
%     
% end
% 
% 
% % Finally concatinate new grid coords
% ngcoords = [coords; ncoords]; 
% 
% %--------------------------------------------------------------------------
% % Plot
% 
% % figure(1)
% % plot(ncoords(:, 1), ncoords(:, 2), 'x');    
% % axis equal; 
% % hold on
% % plot(coords(:, 1), coords(:, 2), 'o'); 

%--------------------------------------------------------------------------
% Create new grid triangles
% 
% % No need in tesselations for triangleGrid()
% %dt = delaunayTriangulation(ngcoords(:, 1), ngcoords(:, 2));  
% 
% G = triangleGrid(ngcoords); 
% G = computeGeometry(G); 