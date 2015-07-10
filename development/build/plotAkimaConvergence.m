clear;
%%
load('bivarAkimaInterpTest.mat');
load('bivarNativeInterpTest.mat');

figure(1);clf;
subplot(1,2,1);
loglog(bivarConv(:,1),bivarConv(:,4),'k-x',...
       bivarConv(:,1),bivarConv(:,2),'b-o',...
       bivarConv(:,1),bivarConv(:,3),'r-*','MarkerSize',12,'LineWidth',2);
hold on;
loglog(bivarConv(:,1),4000* bivarConv(:,1).^(-1),'k--',...
       bivarConv(:,1),100*bivarConv(:,1).^(-3/2),'k-');
set(gca,'FontSize',14);
xlabel('N particles');
ylabel('rel. error');
legend('scalar interp','scalar grad','scalar lap','O(N^{-1})','O(N^{-3/2})','Location','SouthWest')
xlim([10,1e6]);
ylim([1e-6,10]);
title('Akima convergence tests');   

subplot(1,2,2);
loglog(nP, interpErr, 'k-x', ...
       nP, estGradErr,'b-o', ...
       nP, estLapErr, 'r-*', 'MarkerSize',12,'LineWidth',2);
hold on;
loglog(nP, 8000*nP.^(-1),'k--',...
       nP, 1000*nP.^(-3/2),'k-');
set(gca,'FontSize',14);
xlabel('N particles');
ylabel('rel. error');   
title('LPM convergence tests');   
xlim([10,1e6]);
ylim([1e-6,10]);