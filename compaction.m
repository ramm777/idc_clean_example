% - Model that simulates fracture compaction with permeability calculation at each simulation step. 
% - Grid can be selected from a lot of pre-defined grids using the natural Camel mudrock. Usually hard-coded grids are used.  
% - Stresses are measured as a mean resultant at the boundary for Svk but they very close to the force_bc in Pa
% - The simulation will run until specified 'final_stress'
% - Mechanics stress convention, unless specified on the figure. 
% - Author Amanzhol Kubeyev, aman85work@gmail.com


clear all 
close all
mrstModule add vemmech mrst-gui coarsegrid incomp coarsegrid % for plotCellNumbers(G)
disp(['Model: ', mfilename]);
ss = get(groot, 'Screensize'); 
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontSize', 18);


%--------------------------------------------------------------------------
% SIMULATION CONTROL (select manually)

 
contactSolver = 'LMd';                   % 1) LM - for Full Lagrange Multipliers or 2) LMd for Directional Lagrange Multipliers 
grid_type     = 'GaussGrid2';            % Carmel1a, GaussGrid2, Carmel1a, Carmel1b-JRCreduced
ubc =  9000*repmat(2.3253e-07, 1, 50);
ubc(1) = ubc(1) / 10;            
totsteps = length(ubc);
final_stress = -2.5e7;                   % At what stress simulation stops
gapr = sqrt((0.000001*milli*darcy)*12);  % Residual gap based on parallel plate and matrix perm. Is removed from mech apertures. 

save_.Ghist      = false;                    
save_.Gstokes    = false;                % This also saves fluid velocities  
save_.stresshist = true; 
save_.results    = 'none';               % 'all' (inc Ghist), 'only_roughness', 'only_stressperm', 'standard'
save_.video      = false;

calc_.Stokes     = true; 
calc_.JRC        = false; 

LASTN = maxNumCompThreads(1);            % Set the maximum number in multi-threading (parallel computation)   



%% GRID + INTRAGRID


%[G, ginfo] = createGridCarmel(grid_type, 'shiftstep', 1);
[G, ginfo] = createGridSynthetic(grid_type);
fraccells  = ginfo.fraccells; 
demcells   = ginfo.demcells;
rhsbody    = ginfo.rhsbody;
lhsbody    = ginfo.lhsbody;


prop = intracellProp(G, fraccells);   
neiintracells  = prop.neiintracells;  


fnodes = findCellNodes(G, fraccells); % fracture nodes
assert(strcmp(G.type{1}, 'tensorGrid'), 'WARNING: perm calculation is for MRST grids, see fnodes');
fnodes = sort(fnodes, 2, 'ascend');   % fnodes must always be sorted: [LHS body,  RHS body] or [LOWER, UPPER]
[B, I] = sort(fnodes(:, 1));          % fnodes also must be sorted ascending of the first comumn for lenA calculation    
fnodes = [B, fnodes(I, 2)];  


G = moveBodiesCloser(G, gapr, rhsbody, lhsbody, fnodes); 
Ginit = G;                            % Save initial Grid, with distance btwn bodies minimalised
 

% figure(2); plotGridDiscrete(G, rhsbody, lhsbody); 
    
    
clear prop grid_type B I 



%% ELASTICITY parameters

opt = struct( 'grid_type' , 'square',     ...
              'E'         , 23e9,         ...  % Young's modulus, Carmel mudrock
              'nu'        , 0.31);             % Poisson's ratio, Carmel mudrock

         
Ev  = repmat(opt.E, G.cells.num, 1);        
nuv = repmat(opt.nu, G.cells.num, 1);

Ev(fraccells) = 0;        % Put E=0 to intracells for discontinuities

C    = Enu2C(Ev, nuv, G); % Isotropic



%% GRAVITY as loading

gravity reset off
density = 0; % 2260 kg/m^3 
grav    = norm(gravity());   
load_grav    = @(x) -(grav*density)*repmat([0, 1], size(x, 1), 1);

clear density grav


%% BC Displacement (Diriclet)


bc = findBoundaryFacesAndCreateBC(G);
bc = findBoundaryNodesAndCreateBC(G, bc); 


% Set masks for BC: if 'mask' = false, then imposed bc is ignored   
bc_el_sides{1} = bc{1}; % Left
bc_el_sides{1}.el_bc.disp_bc.mask(:, :) = true;
bc_el_sides{2} = bc{2}; % Right
bc_el_sides{2}.el_bc.disp_bc.mask(:, :) = true;
bc_el_sides{3} = bc{3}; % Lower
bc_el_sides{3}.el_bc.disp_bc.mask(:, :) = false;
bc_el_sides{4} = bc{4}; % Upper
bc_el_sides{4}.el_bc.disp_bc.mask(:, :) = false;

% Specify disp value for each side BC on x and y directions [x y]
bc_el_sides{1}.el_bc.disp_bc.uu(:, :)      = [0 0].*ones(numel(bc_el_sides{1}.el_bc.disp_bc.nodes), 2);
bc_el_sides{1}.el_bc.disp_bc.uu_face(:, :) = [0 0].*ones(numel(bc_el_sides{1}.el_bc.disp_bc.faces), 2);
bc_el_sides{2}.el_bc.disp_bc.uu(:, :)      = [0 0].*ones(numel(bc_el_sides{2}.el_bc.disp_bc.nodes), 2); 
bc_el_sides{2}.el_bc.disp_bc.uu_face(:, :) = [0 0].*ones(numel(bc_el_sides{2}.el_bc.disp_bc.faces), 2); 
bc_el_sides{3}.el_bc.disp_bc.uu(:, :)      = [0 0].*ones(numel(bc_el_sides{3}.el_bc.disp_bc.nodes), 2);
bc_el_sides{3}.el_bc.disp_bc.uu_face(:, :) = [0 0].*ones(numel(bc_el_sides{3}.el_bc.disp_bc.faces), 2);
bc_el_sides{4}.el_bc.disp_bc.uu(:, :)      = [0 0].*ones(numel(bc_el_sides{4}.el_bc.disp_bc.nodes), 2);
bc_el_sides{4}.el_bc.disp_bc.uu_face(:, :) = [0 0].*ones(numel(bc_el_sides{4}.el_bc.disp_bc.faces), 2);


disp_bc = collectUBCs(bc_el_sides);
bc      = getSingleBoundaryCells(G, bc, 'Left');  

% figure(33); plotFaces(G, [bc{1}.face; bc{2}.face; bc{3}.face; bc{4}.face]); % Check manually selected faces. Color doesn't work for some reason 



%% BC Force (Neumann)

faces = [bc{1}.face; bc{2}.face; bc{3}.face; bc{4}.face];

% Specify forces to be applied on x and y directions [x y] in Pa/m^3  
force1 = [0 0].*ones(numel(bc{1}.face), 2); 
force2 = [0 0].*ones(numel(bc{2}.face), 2); 
force3 = [0 0].*ones(numel(bc{3}.face), 2);
force4 = [0 0].*ones(numel(bc{4}.face), 2);
forces = [force1; force2; force3; force4]; 
force_bc = struct('faces', faces, 'force', forces); 

el_bc = struct('disp_bc', disp_bc, 'force_bc', force_bc); % Final structure for the boundary conditions

clear force1 force2 force3 force4 forces force_bc


%% SIMULATION

disp('---============ STARTING IDC SIMULATION ============---')
s = 0; % Simulation step
tic;


% Pre-allocate 
stress = 0;
connodes = [];
gridupdated = false; 
prevnodes = [];
F = {};
jrc = nan(1, totsteps); 
z2 = nan(1, totsteps);
connodesnum  = nan(1, totsteps); 
allmechapers = nan(length(fnodes(:, 1)) - 1, totsteps); 
lensA        = nan(length(fnodes(:, 1)) - 1, totsteps);
stresshistbc = nan(length(fnodes(:, 1)) - 1, totsteps);
perm         = nan(1, totsteps); 
permf        = nan(1, totsteps);    
if save_.Ghist   == true; Ghist   = cell(1, totsteps); end 
if save_.Gstokes == true; Gstokes = cell(1, totsteps); end 
if save_.Gstokes == true; vel     = cell(1, totsteps); end   
if save_.stresshist == true; stresshist = cell(1, totsteps); end


[connodes_new, ~] = conDetect2(G, neiintracells); % Check before simulation that Grid is ok: Contact Detection Algorithm
assert(~any(connodes_new(:)), 'Grid is wrong, contact before the simulation');
clear connodes1


while s < totsteps
    
        
    % ---------------------------------------------------------------------
    % Modify BC to Dirichlet
        
    el_bc.disp_bc.mask(1: numel(bc{1}.el_bc.disp_bc.nodes), : ) = true;     % Turn-on Dirichlet BC manually in the bc{1}
    el_bc.disp_bc.uu(1: numel(bc{1}.el_bc.disp_bc.nodes), 1 ) = ubc(s + 1); % Modify BC 
    el_bc.force_bc.force = 0.*el_bc.force_bc.force;                         % Turn off Forces bc
       
    
    % ---------------------------------------------------------------------
    % SOLVE ORIGINAL
    
    [uu, extra] = VEM_linElast_f(G, C, el_bc, load_grav, 'experimental_scaling', true);
 
             
    %% CONTACT DETECTION 
    
    Gcd = G;        
    Gcd.nodes.coords = Gcd.nodes.coords + uu;           
    [connodes_new, ~] = conDetect2(Gcd, neiintracells); 
    conflag = any(connodes_new(:));                   
    
    
    connodes = [connodes; connodes_new]; 
    connodes = unique(connodes, 'rows');  
    clear connodes_new
     
    
    if conflag == false; gridupdated = false; end % If there is no contact, gridupdated must stay false to update grid at the end 
    
   
    while conflag
    %% CONTACT INTERACTION 
        
    
        [connodes, lagnum, g] = dataForContact(G, ginfo, connodes, gapr); % Prepare data to solve the Contact Problem
        %assert(all(ubc(s + 1) > g), 'Your ubc must be > gap, otherwise equations in Contact Problem may be inaccurate')
        el_bc.disp_bc.uu(1: numel(bc{1}.el_bc.disp_bc.nodes), 1 ) = ubc(s + 1);         
        checkBypassedPenetration(G, connodes, fraccells);   
        dispContactIterationStep(s, connodes); 
        connodesnum(s+1) = numel((connodes(:,1))); 
        

        % -----------------------------------------------------------------           
        % SOLVE Contact Problem 
        
        switch contactSolver
            
            case('LM') % Full Lagrange Multipliers (inc gaps)
                [uu, extra] = VEM_linElast_LM(G, C, el_bc, load_grav, 'experimental_scaling', true, 'lagnum', lagnum, ...
                                              'connodes', connodes, 'gap', g, 'ginfo', ginfo);

    
                lagmultXY = reshape(extra.disc.lagmult, G.griddim, [])'; % Get contact forces [LMx, LMy] 
                if any(lagmultXY(:, 1) < 0); warning('Negative Lagmult'); end

                
            case('LMd') % SOLVE Directional Lagrange Multipliers (inc gaps)
                [uu, extra] = VEM_linElast_LMd(G, C, el_bc, load_grav, 'experimental_scaling', true, 'lagnum', lagnum, ...
                                               'connodes', connodes, 'gap', g, 'direction', 'x', 'ginfo', ginfo);

             
                lagmultX = extra.disc.lagmult; % Get contact forces [LMx, LMy] 
                if any(lagmultX(:, 1) < 0); warning('Negative Lagmult'); end
            
            
             otherwise
                error('Select function to solve the Contact Problems: LM or LMd')
        end
        
        %------------------------------------------------------------------
        % CONTACT DETECTION (As none of the cells can be incontact after the penetration)
          
        Gcd = G; 
        Gcd.nodes.coords = Gcd.nodes.coords + uu; 
        [connodes_new, ~] = conDetect2(Gcd, neiintracells);
        conflag = any(connodes_new(:)); 
        
        
        % Continue, if conDetect keeps finding exactly same connodes as in previous
        if sum(prevnodes(:)) == sum(connodes_new(:));  conflag = false;  end
        prevnodes = connodes_new; 
           
        
        %------------------------------------------------------------------
        % If penetration after solving a Contact Problem => concatinate connodes
        
        if conflag 
            disp('WARNING: Penetration after Contact Solution=> adding additional connodes');
            connodes = [connodes_new; connodes]; 
        
        %------------------------------------------------------------------
        % If contact is solved with no penetration => success   
        
        else      
            
            
            s = s + 1; 
            
            
            % Update grid and re-compute VEM and geometry (volumes, centroids)
            G.nodes.coords = G.nodes.coords + uu;
            G = createAugmentedGrid(G); 
            G = computeGeometry(G);                                      
            gridupdated = true;
            if save_.Ghist == true; Ghist{s} = G; end 
            
            
            [allmechapers, lensA] = mechapers(G, fnodes, s, allmechapers, lensA); 
            
            if calc_.JRC == true
                Lcoords = G.nodes.coords(fnodes(:, 1), :);
                [z2(s), jrc(s)] = JRCtse(Lcoords(:, 1), Lcoords(2, 2)- Lcoords(1, 2), 1, 'mode', 'reproduce');                     
            end 
            
             
            %-------------------------------------------------------------------------
            % STOKES permeability
            
            if calc_.Stokes == true
            
                Gfrac = extractSubgrid(G, fraccells);
                matrixvol = sum(G.cells.volumes(lhsbody)) + sum(G.cells.volumes(rhsbody));
                [perm(s), permf(s), finfo] = permStokes(Gfrac, 'matrixvol', matrixvol, 'getFlowInfo', true); % in mD
                
                if save_.Gstokes == true; vel{s}     = finfo.vel; end
                if save_.Gstokes == true; Gstokes{s} = finfo.G;   end 
            
            end
                        
            %-------------------------------------------------------------------------
            % Calculate stress and plot
            
            
            [stress1, ~] = calculateStressVEM(G, uu, extra); % Note, stress-xy is different as /2
            stress = stress + stress1;                       % Cumulative stress 
            stresshistbc(:, s) = stress(bc{1}.cells, 1);     % Record stress-x at each iteration
            if save_.stresshist == true; stresshist{s} = stress; end
            clear stress1; 
            
            figure(4); clf; set(4, 'Position', [ss(1)*2310, ss(2)- 100, ss(3)/2, ss(4)/1.22]);
            plotCellData(G, stress(:, 1), demcells, 'FaceAlpha', 0.9, 'EdgeAlpha', 0.01); 
            axis tight; % zoom(10); colormap jet;  
            set(gca,'visible','off'); 
            title(['Simulation, step: ', num2str(s-1)]);
            hcb=colorbar; hcb.Title.String = "\sigma_{x}, Pa";
            caxis([-110e7, 0]);
            hold on
            
            
            if save_.video == true; F{s} = getframe(4); end
            if s == 1; saveFirstSolutionFig(G, stress, demcells, s); end
            
            
            simtime(s) = toc;     
            
            break 
            
        end
        
        
    end
    
    %--------------------------------------------------------------
    % Continue simulation unless specified cumulative stress at bc is reached 
    
    if s ~= 0
        stressm = mean(stresshistbc(:, s)); 
        disp(['Average bc stress, 1e7Pa= ', num2str(stressm/1e7)])
        if s >= totsteps && stressm > final_stress
            ubc = [ubc, ubc(end)];   
            totsteps = totsteps + 1; 
        end
    end
    
    %--------------------------------------------------------------
    % If grid is updated during solving Contact Problem, continue
    
    if gridupdated 
        continue
    else           
        G.nodes.coords = G.nodes.coords + uu;
        G = createAugmentedGrid(G);
        G = computeGeometry(G);
    end
    
    
end


total_ubc = cumsum(ubc); 
disp(['Total ubc= ', num2str(total_ubc(end))])
disp(['Simulation time (min)= ', num2str(simtime(end)/60)])



%% APERTURES / PLOT


allmechapers(allmechapers ~= 0) = allmechapers(allmechapers ~= 0)  - gapr; % Delete gapr from perm calculation where aperture is not zero
meanapers     = mean(allmechapers);                                        % Mean apertures 
mstresshistbc = -1*(mean(stresshistbc, 1));                                % Mean stress on outer bc (geomech convention)


figure(6); set(6, 'Position', [ss(1)*2150, ss(2)+20, ss(3)/3, ss(4)/2]);  
plot(mstresshistbc, perm, '--o'); hold on; title('Geomechanics convention'); 
xlabel('Stress (geomechanics convention), Pa'); ylabel('Permeability, mD' ); hold off


%% SAVE

if save_.video == true; createVideo(F, mfilename); end

clear LASTN C calc_ conflag connodes_new disp_bc Ev F faces final_stress
clear gridupdated Lcoords load_grav nuv opt simtime ss 


switch save_.results  
    case ('standard')        
        save(mfilename)     
        figure(77); clf;      % Save figure of the last timestep
        plotCellData(G, stress(:, 1), demcells, 'FaceAlpha', 0.9, 'EdgeAlpha', 0.01); colorbar;
        title(['SOLVED Contact, Step= ', num2str(s-1)]); axis tight; 
        saveas(77, [mfilename, '-final-solution']);
    case ('all')
        save('Gstokes.mat', 'Gstokes', '-v7.3'); 
        save('Ghist.mat',   'Ghist',   '-v7.3'); 
        save(mfilename, '-regexp', '^(?!(Gstokes|Ghist)$).');
    case ('only_roughness') 
        save(mfilename, 'z2', 'jrc', 'mstresshistbc');
    case ('only_stressperm') 
        save(mfilename, 'mstresshistbc', 'perm', 'permf');
    case ('none')
        fprintf('No results are saved')        
end


function saveFirstSolutionFig(G, stress, demcells, s)

        figure(7); clf;      % Save figure of the last timestep
        plotCellData(G, stress(:, 1), demcells, 'FaceAlpha', 0.9, 'EdgeAlpha', 0.01); colorbar;
        title(['SOLVED Contact, Step= ', num2str(s-1)]); axis tight; 
        saveas(7, [mfilename, '-initial-solution']);
        close(7);
        
end
