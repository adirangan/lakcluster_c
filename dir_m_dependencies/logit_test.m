% 
% M. Farbood, April 30, 2012
%
% This script provides an example of logistical regression in Matlab.
% It also calculates some additional useful statistcs:
%
%
% The data used in this script is from the article, 
% "An Introduction to Logistic Regression Analysis and Reporting" in the 
% Journal of Educational Research by Peng et al., 2002.
%
% N.B. 1: The score test (Lagrange multiplier test) is not calculated here.
%
% N.B. 2: for calculating the likelihood ratio and the Wald tests of the model, I
% include the function calls (commented out) for the lratiotest 
% and waldtest functions from the Matlab Econometrics toolbox.
%
% The Hosmer-Lemeshow statistic is calculated in a separate script, HLtest.m.
%
 
function logit_test

    % This data is from the article, "An Introduction to Logistic Regression
    % Analysis and Reporting" in the Journal of Educational Research,
    % Peng et al., 2002.

    % First column is reading score; second is gender (0 or 1); last is outcome (0 or 1)
    load remedial_reading_data
    M = remedial_reading_data;

    X = M(:,1:2); % predictors: gender and reading score
    Y = M(:,3); % outcome

    readingScore = X(:,1);  
    gender = X(:,2);  
    n = length(readingScore);

    varnames = {'Constant', 'Reading', 'Gender'};
    alpha = 0.05;

    % Logistic regression
    % Note that 'link', 'logit' are not actually required since 'binomial'
    % defaults to that anyway.
    [b,dev,stats] = glmfit(X,Y, 'binomial', 'link', 'logit');
    oddsRatios = exp(stats.beta);
    betas = stats.beta; % for use later
    Wstats = stats.t.^2;

    % Print stats corresponding to upper portion of Table 3 from Peng et al., 2002
    fprintf('------------------------------------------------------------------\n');
    fprintf('%-10s %-8s %-8s %-8s %-4s %-8s %-10s\n', 'Predictor', 'Beta', 'SE', 'Wald', 'df', 'p', 'Odds ratio'); 
    fprintf('------------------------------------------------------------------\n');
    fprintf('%-10s %-8.4f %-8.4f %-8.4f %-4d %-8.4f %-10s\n', varnames{1}, stats.beta(1), ...
        stats.se(1), Wstats(1), 1, stats.p(1), 'N/A');
    fprintf('%-10s %-8.4f %-8.4f %-8.4f %-4d %-8.4f %-10.4f\n', varnames{2}, stats.beta(2), ...
        stats.se(2), Wstats(2), 1, stats.p(2), oddsRatios(2));
    fprintf('%-10s %-8.4f %-8.4f %-8.4f %-4d %-8.4f %-10.4f\n', varnames{3}, stats.beta(3), ...
        stats.se(3), Wstats(3), 1, stats.p(3), oddsRatios(3));

    % Wald test of model
    R = eye(3);
    R = R(2:end,:); %i.e., R = [0 1 0; 0 0 1];
    r = [stats.beta(2) stats.beta(3)]';
    cov = stats.covb;
    Wdf = length(r);

    % Following line calulate's Wald chi-square with the Econometrics toolbox function:
    % [h,Wp,Wstat,cValue] = waldtest(r,R,cov)
    Wstat = r'*inv(R*cov*R')*r;
    Wp = 1 - chi2cdf(Wstat,Wdf);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Some of the following code I got from this script:
    % http://springer.bme.gatech.edu/Ch17.Logistic/Logisticmat/dmdreg.m

    lin = b(1) + b(2)*readingScore + b(3)*gender;
    phat = exp(lin)./(1 + exp(lin));

    % Log likelihood
    warning off
    loglik = sum(Y .* lin - log(1 + exp(lin)));
    X0 = zeros(size(Y));
    [b0, dev0, stats0] = glmfit(X0,Y,'binomial','link','logit');
    loglik0 = sum(Y .* b0(1) - log(1 + exp(b0(1))));

    % Likelihood ratio test
    uLL = loglik;
    rLL = loglik0;
    LRdf = length(betas) - 1;

    % The following line is the call to the function lratiotest from the
    % Econometrics toolbox
    %[h,LRp,LRstat,cValue] = lratiotest(uLL,rLL,LRdf)
    LRstat = 2*(uLL - rLL);
    LRp = 1 - chi2cdf(LRstat,LRdf);

    % Score test using Econometrics toolbox (haven't implemented this yet)
    %[LMh,LMp,LMstat,LMcv] = lmtest(Rscore,REstCov1,dof)

    % Hosmer-Lemeshow goodness-of-fit test
    [HLstat HLp HLdf contingencyTable] = HLtest(M, betas',10);

    fprintf('------------------------------------------------------------------\n');
    fprintf('%-28s %-8.4f %-4d %-8.4f\n', 'Likelihood ratio test', LRstat, LRdf, LRp);
    fprintf('%-28s %-8.4f %-4d %-8.4f\n', 'Wald test', Wstat, Wdf, Wp);
    fprintf('%-28s %-8.4f %-4d %-8.4f\n', 'Hosmer & Lemeshow test', HLstat, HLdf, HLp);
    fprintf('------------------------------------------------------------------\n');

    % McFadden Pseudo R^2 
    pseudor2 = -2 * (loglik0 -  loglik)/(- 2* loglik0); % equivalent to line below
    mcfadden = 1 - loglik/loglik0; % equivalent to 1 - dev/dev0
    mcfaddenadj = 1 - (loglik - length(betas))/loglik0;  %counterpart to R2adj

    % Effron R^2
    effron = 1 -  sum((Y - phat).^2)/sum((Y - sum(Y)/n).^2);

    % Cox and Snell R^2
    coxsnell = 1 - (exp(loglik0)/exp(loglik))^(2/n);

    % Nagelkerke R^2
    nagelkerke = (1 - (exp(loglik0)/exp(loglik))^(2/n) ) / (1 - exp(loglik0)^(2/n));


    fprintf('Cox and Snell R2 =  %.4f\n', coxsnell); 
    fprintf('Nagelkerke R2 =  %.4f\n', nagelkerke); 
    fprintf('McFadden Pseudo R2 =  %.4f\n', mcfadden); 
    fprintf('McFadden Adjusted R2 =  %.4f\n', mcfaddenadj); 
    fprintf('Effron R2 =  %.4f\n', effron); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Some extra info:
    % Prints contingency talb for Hosmer-Lemeshow test
    fprintf('\nContingency Table for Hosmer-Lemeshow goodness-of-fit test:\n')
    fprintf('--------------------------------------------\n');
    fprintf('%-6s %-12s %-10s %s\n', 'Group', 'Observations', 'Expected', 'Observed');
    fprintf('--------------------------------------------\n');
    fprintf('%-6d %-12d %-10.2f %d\n', contingencyTable');
