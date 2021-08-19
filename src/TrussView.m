function TrussView(xH,fH,Time,fem,Filter)
%% --------------------------------------------- Apply visualization filter
x = xH(:,end);
inactive = find(x<Filter*max(x));
active = setdiff((1:length(x))',inactive);
xEff = x; xEff(inactive) = 0;
%% -------------------------------------------------------- Summary results
fprintf('Number of bars: %d of %d\n', nnz(xEff), length(x));
fprintf('Objective function value: %g\n',fH(end));
fprintf('Number of iterations: %d\n',length(fH)-1);
fprintf('Elapsed time (TrussTop): %g secs.\n', Time);
%% --------------------------------------------- Objective function history
figure('Color',[1 1 1]); box on; grid on; hold on;
plot(fH,'r-o','MarkerEdgeColor','b','MarkerFaceColor','w','MarkerSize',4);
xlabel('Iteration'); ylabel('Objective function value');
title('Objective function history'); drawnow;
% %% --------------------------------------------------- Active areas history
% figure('Color',[1 1 1]); box on; grid on; hold on;
% for i = active'
%     plot(xH(i,:),'k-');
% end
% xlabel('Iteration'); ylabel('area value');
% title('Active areas history'); drawnow;
% %% ------------------------------------------------- Inactive areas history
% figure('Color',[1 1 1]); box on; grid on; hold on;
% for i = inactive'
%     plot(xH(i,:),'k-');
% end
% xlabel('Iteration'); ylabel('area value');
% title('Inactive areas history'); drawnow;
% %% ------------------------------------------------- Member areas histogram
% figure('Color',[1 1 1]); box on; grid on; hold on;
% bar(sort(xEff(xEff>0))./max(x),'k'); xlim([.5 nnz(xEff)+.5]);
% xlabel('Sorted active member'); ylabel('Normalized area');
% title('Active member areas'); drawnow;
%% --------------------------------------------------------- Final topology
figure('Color',[1 1 1]); box off; hold on; grid off; axis equal;
for e = find(xEff>0)'
    X = fem.Node(fem.Element(e,:),1); Y = fem.Node(fem.Element(e,:),2);
    clr = 'b'; if fem.stress(e)<0, clr='r'; end
    plot(X,Y,[clr '-'],'LineWidth',xEff(e)/max(x)*5);
end
for e = find(xEff>0)'
    X = fem.Node(fem.Element(e,:),1); Y = fem.Node(fem.Element(e,:),2);
    plot(X,Y,'ko','MarkerFaceColor','w','Markersize',5);
end
title('Final topology'); axis off; drawnow; hold off;