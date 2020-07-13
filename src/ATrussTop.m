%% --------------------------------------------------------------- TRUSSTOP
function [x_hist,f_hist,fem] = ATrussTop(fem,opt)
Iter=0; ConsecOCIter=0; Tol=opt.Tol; Change=2*Tol; x=opt.xIni; dfdx_ = [];
x_hist = zeros(fem.NElem,opt.MaxIter+1);
f_hist = zeros(1,opt.MaxIter+1);
while (Iter<opt.MaxIter) && (Change>Tol)
    Iter = Iter + 1;
    %Compute cost functionals and analysis sensitivities
    [f,dfdx,fem] = ObjectiveFnc(fem,x);
    [g,dgdx,fem] = ConstraintFnc(fem,x,opt.VolMax);
    %Store design vars. and obj. function history
    x_hist(:,Iter) = x; f_hist(Iter) = f;
    %Update design variable and analysis parameters
    %[x,Change,x_,dfdx_] = UpdateScheme(dfdx,g,dgdx,x,dfdx_,x_,Iter,opt);
    [x,dfdx_,Change,ConsecOCIter] = HybridUpdateScheme...
        (x_hist,f_hist,dfdx,g,dgdx,dfdx_,Iter,ConsecOCIter,opt);
end
% Calculate stress and internal force
fem = Stress(fem,FEAnalysis(fem,x));
%Store last design var. and obj. function value, then resize history arrays
x_hist(:,Iter+1) = x; [f_hist(Iter+1),~,fem] = ObjectiveFnc(fem,x);
x_hist(:,Iter+2:end) = []; f_hist(Iter+2:end) = [];
%% ----------------------------------------------------- OBJECTIVE FUNCTION
function [f,dfdx,fem] = ObjectiveFnc(fem,x)
[U,fem] = FEAnalysis(fem,x);
f = dot(fem.F,U);
temp = cumsum(-U(fem.i).*fem.k.*U(fem.j));
temp = temp(cumsum(fem.ElemNDof.^2));
dfdx = [temp(1);temp(2:end)-temp(1:end-1)];
%% ---------------------------------------------------- CONSTRAINT FUNCTION
function [g,dgdx,fem] = ConstraintFnc(fem,x,VolMax)
g = dot(fem.L,x)-VolMax;
dgdx = fem.L;
%% --------------------------------------------- OPTIMALITY CRITERIA UPDATE
function [xNew,Change,x_,dfdx_] = UpdateScheme...
    (dfdx,g,dgdx,x,dfdx_,x_,Iter,opt)
xMin=opt.xMin; xMax=opt.xMax; move=opt.OCMove*(xMax-xMin); eta=opt.OCEta;
if Iter>1 && opt.Adapt
    eta=1./(1-max(min(1+log(dfdx./dfdx_)./log(x./x_),-.1),-15));
end
Bm = -dfdx./dgdx; l1=0; l2=1.2*max(Bm);
while l2-l1 > 1e-10*(1+l2)
    lmid = 0.5*(l1+l2); B = Bm/lmid;
    xCnd = xMin+(x-xMin).*B.^eta;
    xNew = max(max(min(min(xCnd,x+move),xMax),x-move),xMin);
    if (g+dgdx'*(xNew-x)>0), l1=lmid; else, l2=lmid; end
end
x_ = x; dfdx_ = dfdx; Change = max(abs(xNew-x)./(1+x));
%% --------------------------------------------------- HYBRID UPDATE SCHEME
function [xNew,dfdx,Change,ConsecOCIter] = HybridUpdateScheme...
    (xHist,fHist,dfdx,g,dgdx,dfdx_,Iter,ConsecOCIter,opt)
% init Iteration as 'OC' and current design variables values as 'x'
OCIter = true; x = xHist(:,Iter); xNew = x;
% check for asymptotic flag and minimum consecutive OC Iterations
if opt.Asymptotic && ConsecOCIter>opt.Am
    % get last 'm' obj function values and its decrements
    Fun = fHist(Iter-opt.Am+1:Iter);
    %     Decrements = Fun(1:end-1)-Fun(2:end);
    Decrements = fHist(Iter-opt.Am:Iter-1) - Fun;
    % check if all decrements are positive and sufficiently small
    %     if all(Decrements>0 && Decrements<opt.Atol*(1+Fun(2:end)))
    if all(Decrements>0 & Decrements<opt.Atol*(1+Fun))
        % get last 'n' increments and design variable values
        LastIncrements = Iter-opt.An+1:Iter;
        LastX = xHist(:,LastIncrements);
        % using Linear Least Squares (X = A + B/Inc)
        %   M * C = Y =>  (M'*M)*C = (M'*Y) => C = (M'*M)\(M'*Y)
        M = ones(opt.An,2); M(:,2) = 1./LastIncrements'; Y = LastX';
        MT = M'; MTM = MT*M;
        
        %   when Increment -> infinity => X -> A
%         C = MTM\(MT*Y); xNew = max(0, C(1,:)');
        for i = 1:length(x), Yi = Y(:,i);
            ci = MTM\(MT*Yi);
            xNew(i) = max(0, ci(1));
        end
        
        % mark iteration as Asymptotic
        OCIter = false; ConsecOCIter = 0;
        fprintf('[Asymptotic] Iter: %d\n', Iter);
    end
end
% if Asymptotic scheme cannot be applied, use OC
if OCIter, ConsecOCIter = ConsecOCIter + 1;
    xNew = OC(xHist,dfdx,dfdx_,g,dgdx,Iter,opt);
end
% calc Change on design variables
Change = max(abs(xNew-x)./(1+x));
%% --------------------------------------------- OPTIMALITY CRITERIA UPDATE
function xNew = OC(x_hist,dfdx,dfdx_,g,dgdx,Iter,opt)
xMin=opt.xMin; xMax=opt.xMax; x=x_hist(:,Iter);
move=opt.OCMove*(xMax-xMin); eta=opt.OCEta;
if Iter>1 && opt.Adapt
    x_=x_hist(:,Iter-1);
    eta=1./(1-max(min(1+log(dfdx./dfdx_)./log(x./x_),-.1),-15));
end
Bm = -dfdx./dgdx; l1=0; l2=1.2*max(Bm);
while l2-l1 > 1e-10*(1+l2)
    lmid = 0.5*(l1+l2); B = Bm/lmid;
    xCnd = xMin+(x-xMin).*B.^eta;
    xNew = max(max(min(min(xCnd,x+move),xMax),x-move),xMin);
    if (g+dgdx'*(xNew-x)>0), l1=lmid; else, l2=lmid; end
end
%% ------------------------------------------------------------ FE-ANALYSIS
function [U,fem] = FEAnalysis(fem,x)
if ~isfield(fem,'k')
    fem.ElemNDof = 4*ones(size(fem.Element,1),1);
    fem.i = zeros(sum(fem.ElemNDof.^2),1);
    fem.j=fem.i; fem.k=fem.i; fem.e=fem.i;
    idx = 0;
    for el = 1:fem.NElem
        Ke=LocalK(fem,el); NDof = fem.ElemNDof(el);
        eDof = reshape([2*fem.Element(el,:)-1;2*fem.Element(el,:)],NDof,1);
        I=repmat(eDof ,1,NDof); J=I';
        fem.i(idx+1:idx+NDof^2) = I(:); fem.j(idx+1:idx+NDof^2) = J(:);
        fem.e(idx+1:idx+NDof^2) = el; fem.k(idx+1:idx+NDof^2) = Ke(:);
        idx = idx + NDof^2;
    end
    NLoad = size(fem.Load,1);
    fem.F = zeros(2*fem.NNode,1);  %external load vector
    fem.F(2*fem.Load(1:NLoad,1)-1) = fem.Load(1:NLoad,2);  %x-crdnt
    fem.F(2*fem.Load(1:NLoad,1))   = fem.Load(1:NLoad,3);  %y-crdnt
    NSupp = size(fem.Supp,1);
    FixedDofs = [fem.Supp(1:NSupp,2).*(2*fem.Supp(1:NSupp,1)-1);
        fem.Supp(1:NSupp,3).*(2*fem.Supp(1:NSupp,1))];
    FixedDofs = FixedDofs(FixedDofs>0);
    AllDofs   = 1:2*fem.NNode;
    fem.FreeDofs = setdiff(AllDofs,FixedDofs);
end
K = sparse(fem.i,fem.j,x(fem.e).*fem.k);
U = zeros(2*fem.NNode,1);
U(fem.FreeDofs,:) = K(fem.FreeDofs,fem.FreeDofs)\fem.F(fem.FreeDofs,:);
%% ----------------------------------------------- ELEMENT STIFFNESS MATRIX
function [Ke] = LocalK(fem,el)
A = fem.localIntMap(el,:);
Ke = fem.E/fem.L(el)*(A'*A);
%% ----------------------------------------------------------------- STRESS
function fem = Stress(fem,U)
if ~isfield(fem,'stress'), fem.stress = zeros(fem.NElem,1); end
U = reshape(U,[2 fem.NNode])';
for el = 1:fem.NElem
    Ue = reshape(U(fem.Element(el,:),:)',[1 4])';
    Delta = fem.localIntMap(el,:)*Ue;
    fem.stress(el) = fem.E*Delta/fem.L(el);
end