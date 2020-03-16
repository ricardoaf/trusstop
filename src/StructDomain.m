function [NODE,ELEM,SUPP,LOAD]=StructDomain(Nx,Ny,Lx,Ly,CS,CL,DrawBC)
if nargin<7 || iesmpty(DrawBC), DrawBC=true;
elseif nargin<6, error('Not enough input arguments.'), end
%% --------------------------------- Generate structured-orthogonal domains
[X,Y] = meshgrid(linspace(0,Lx,Nx+1),linspace(0,Ly,Ny+1));
NODE = [reshape(X,numel(X),1) reshape(Y,numel(Y),1)];
k = 0; ELEM = cell(Nx*Ny,1);
for j=1:Ny
    for i=1:Nx
        k = k+1;
        n1 = (i-1)*(Ny+1)+j; n2 = i*(Ny+1)+j;
        ELEM{k} = [n1 n2 n2+1 n1+1];
    end
end
SUPP = [];
for i = 1:size(CS,1)
    node = FindNodeSet(NODE,CS(i,[1 2]));
    len = length(node);
    SUPP = [SUPP; node CS(i,3)*ones(len,1) CS(i,4)*ones(len,1)];
end
LOAD = [];
for i = 1:size(CL,1)
    node = FindNodeSet(NODE,CL(i,[1 2]));
    len = length(node);
    LOAD = [LOAD; node CL(i,3)*ones(len,1) CL(i,4)*ones(len,1)];
end
%% ----------------------------------------------- Draw boundary conditions
if DrawBC
    figure('Color',[1 1 1]);
    box on; hold on; grid off; axis equal;
    plot(NODE(:,1),NODE(:,2),'k.','MarkerSize',4,'DisplayName','Node');
    SUPP_X  = SUPP(SUPP(:,2)==1 & SUPP(:,3)==0);
    SUPP_Y  = SUPP(SUPP(:,2)==0 & SUPP(:,3)==1);
    SUPP_XY = SUPP(SUPP(:,2)==1 & SUPP(:,3)==1);
    plot(NODE(SUPP_X,1),NODE(SUPP_X,2),'cs','MarkerSize',6,...
        'MarkerFaceColor','c','DisplayName','Support 1,0');
    plot(NODE(SUPP_Y,1),NODE(SUPP_Y,2),'gs','MarkerSize',6,...
        'MarkerFaceColor','g','DisplayName','Support 0,1');
    plot(NODE(SUPP_XY,1),NODE(SUPP_XY,2),'bs','MarkerSize',6,...
        'MarkerFaceColor','b','DisplayName','Support 1,1');
    plot(NODE(LOAD(:,1),1),NODE(LOAD(:,1),2),'ro','MarkerSize',6,...
        'MarkerFaceColor','r','DisplayName','Load');
    set(legend(gca,'show'),'Orientation','horizontal',...
        'Location','southoutside');
    title('Boundary conditions'); drawnow; hold off;
end
%-------------------------------------------------------------------------%
function nset = FindNodeSet(Node,pos)
x = pos(1); y = pos(2);
if x<0 && y>=0
    dist = abs(Node(:,2)-y);
elseif x>=0 && y<0
    dist = abs(Node(:,1)-x);
elseif x<0 && y<0
    dist = zeros(size(Node,1),1);
else
    dist = sqrt((Node(:,1)-x).^2+(Node(:,2)-y).^2);
end
nset = find(dist==min(dist));