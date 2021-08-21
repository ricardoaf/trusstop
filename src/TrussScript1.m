%% ------------------------------------------------------------- INPUT DATA
Nx = 2; Ny = 1;                     % Number of cells
Lx = 2; Ly = 1;                     % Domain size
CS = [0 -1, 1 1];                   % Support data [x y, sx sy; ...]
CL = [Lx 0, 0 -100];                % Load data    [x y, Px Py; ...]
Lvl = 2;                            % Ground structure (GS) level
RestrictDomain = @(nodes,coord)[];  % Define void regions
ColTol = 0.999999;                  % Colinear tolerance for GS
E = 7e+07;                          % Young Modulus
VolMax = Lx*Ly/9000;                % Volume constraint
%% ---------------------------------------------------- CREATE 'fem' STRUCT
[Node,Elem,Supp,Load]=StructDomain(Nx,Ny,Lx,Ly,CS,CL);
[Bars,L,L2I] = GenerateGS(Node,Elem,Lvl,RestrictDomain,ColTol);
fem = struct(...
    'NNode',size(Node,1),...  % Number of nodes
    'NElem',size(Bars,1),...  % Number of elements
    'Node',Node,...           % [NNode x 2] array of nodes
    'Element',Bars,...        % [NElem x Var] cell array of elements
    'L',L,...                 % [NElem x 1] array with length of bars
    'localIntMap',L2I,...     % [NElem x 4] Local to internal coords map
    'Supp',Supp,...           % Array of supports
    'Load',Load,...           % Array of loads
    'E',E...                  % Young's modulus of material
    );
%% ---------------------------------------------------- CREATE 'opt' STRUCT
Area = VolMax/sum(fem.L);
xIni = Area*ones(fem.NElem,1);
opt = struct(...
    'xMin',Area*1e-04,...     % Lower bound for design variables
    'xMax',Area*1e+04,...     % Upper bound for design variables
    'xIni',xIni,...           % Initial design variables
    'VolMax',VolMax,...       % Specified volume fraction cosntraint
    'Tol',1e-08,...           % Convergence tolerance on design vars.
    'MaxIter',4000,...        % Max. number of optimization iterations
    'OCMove',Area*1e+04,...   % Allowable move step in OC update scheme
    'OCEta',0.5, ...          % Exponent used in OC update scheme
    'Adapt',true,...          % Adaptive Modified OC update scheme
    'UpdateScheme','OC'...    % Update scheme
    );
%% --------------------------------------------------------- RUN 'TrussTop'
opt.UpdateScheme = 'OC';
% opt.UpdateScheme = 'DirectUpdate';
opt.Adapt = true;
tic; [xHist,fHist,fem] = TrussTop(fem,opt);
Time = toc;
%% ------------------------------------------------------- CALL 'TrussView'
Filter = 0.01;  % Fraction of maximum area value (for visualization)
TrussView(xHist,fHist,Time,fem,Filter);
