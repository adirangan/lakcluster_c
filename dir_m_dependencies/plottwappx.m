% Script to plot the approximated TW distributions
%
% see paper M.Chiani, "Distribution of the largest eigenvalue for real 
% Wishart and Gaussian random matrices and a simple approximation for the 
% Tracy-Widom distribution", submitted 2012, ArXiv
%
x = -5:.01:5;
[y1,z1] = tracywidom_appx(x,1);
[y2,z2] = tracywidom_appx(x,2);
[y4,z4] = tracywidom_appx(x,4);
figure;
plot(x,y1,x,y2,x,y4);
xlabel('x')
ylabel('pdf(x)')
title('Shifted gamma approximation for {TW}_\beta')
legend('\beta=1','\beta=2','\beta=4')
set(gca,'XTick',-5:1:5)
figure;
plot(x,z1,x,z2,x,z4);
xlabel('x')
ylabel('CDF(x)')
title('Shifted gamma approximation for {TW}_{\beta}')
legend('\beta=1','\beta=2','\beta=4')
set(gca,'XTick',-5:1:5)
