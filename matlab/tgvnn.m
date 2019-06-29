function [L,S] = tgvnn(input)
% Low rank and sparse based reconstructions using primal-dual algorithm
% min_{L,S} = (1/2)*||A(L+S)-B||_F^2+3DTGV(S)+beta*||L||_*
%
% input:
% B: undersampled data
% A: smapling operator
% AT: inverse sampling operator
% alpha: for TV
% beta: for nuclear norm
% maxits: no. of iteration
%
% output:
% L: low rank component
% S: sparse component
%
% Dong Wang
% 06/07/2019

    % load parameters
    A = input.A;
    AT = input.AT;
    B = input.B;
    alpha0 = input.alpha0;
    alpha1 = input.alpha1;
    beta = input.beta;
    sigma = input.sigma;
    tau = input.tau;
    mu = input.mu;
    iter = input.iter;
    reduction = input.reduction;

    % shrinkage parameter
    alpha00 = alpha0;
    alpha10 = alpha1;
    alpha01 = alpha0*reduction;
    alpha11 = alpha1*reduction;

    % initialize primal variables
    [m,n,d] = size(B);

    V = zeros(m,n,d,3); % for TGV component
    L = AT(B); % for low rank component
    L = max(0, real(L));
    S = zeros(m,n,d); % for sparse component

    % initialize dual variables
    p = zeros(m,n,d,3); % for first order derivative
    q = zeros(m,n,d,6); % for second order derivative
    r = zeros(m,n,d); % for fidelity term

    % initialize intermediate variables
    Sbar = S;
    Lbar = L;
    Vbar = V;

    % run the main loop
    upd = textprogressbar(iter);
    for i = 1:iter
      % shrinkage parameter update
      alpha0 = exp(i/iter*log(alpha01) + (iter-i)/iter*log(alpha00));
      alpha1 = exp(i/iter*log(alpha11) + (iter-i)/iter*log(alpha10));
      
      % save variables
      Sold = S;
      Lold = L;
      Vold = V;

      %%% DUAL UPDATE
      % update r (fidelity term proximal)
      r = (r + sigma*(A(Lbar+Sbar)-B))/(1+sigma);
      
      % update p (first order derivative projection)
      Sx = dxp(Sbar);
      Sy = dyp(Sbar);
      St = dtp(Sbar,mu);

      p(:,:,:,1) = p(:,:,:,1) - sigma*(Sx + Vbar(:,:,:,1));
      p(:,:,:,2) = p(:,:,:,2) - sigma*(Sy + Vbar(:,:,:,2));
      p(:,:,:,3) = p(:,:,:,3) - sigma*(St + Vbar(:,:,:,3));

      absp = sqrt(abs(p(:,:,:,1)).^2 + abs(p(:,:,:,2)).^2 + ...
        abs(p(:,:,:,3)).^2);
      denom = max(1,absp/alpha1);
      p(:,:,:,1) = p(:,:,:,1)./denom;
      p(:,:,:,2) = p(:,:,:,2)./denom;
      p(:,:,:,3) = p(:,:,:,3)./denom;
      

      % update q (second order derivative projection)
      grad1 = dxm(Vbar(:,:,:,1));
      grad2 = dym(Vbar(:,:,:,2));
      grad3 = dtm(Vbar(:,:,:,3),mu);
      grad4 = (dym(Vbar(:,:,:,1)) + dxm(Vbar(:,:,:,2)))/2;
      grad5 = (dxm(Vbar(:,:,:,3)) + dtm(Vbar(:,:,:,1),mu))/2;
      grad6 = (dtm(Vbar(:,:,:,2),mu) + dym(Vbar(:,:,:,3)))/2;

      q(:,:,:,1) = q(:,:,:,1) - sigma*grad1;
      q(:,:,:,2) = q(:,:,:,2) - sigma*grad2;
      q(:,:,:,3) = q(:,:,:,3) - sigma*grad3;
      q(:,:,:,4) = q(:,:,:,4) - sigma*grad4;
      q(:,:,:,5) = q(:,:,:,5) - sigma*grad5;
      q(:,:,:,6) = q(:,:,:,6) - sigma*grad6;

      absq = sqrt(abs(q(:,:,:,1)).^2 + abs(q(:,:,:,2)).^2 + ...
         abs(q(:,:,:,3)).^2 + 2*abs(q(:,:,:,4)).^2 + ...
         2*abs(q(:,:,:,5)).^2 + 2*abs(q(:,:,:,6)).^2);
      denom = max(1,absq/alpha0);
      q(:,:,:,1) = q(:,:,:,1)./denom;
      q(:,:,:,2) = q(:,:,:,2)./denom;
      q(:,:,:,3) = q(:,:,:,3)./denom;
      q(:,:,:,4) = q(:,:,:,4)./denom;
      q(:,:,:,5) = q(:,:,:,5)./denom;
      q(:,:,:,6) = q(:,:,:,6)./denom;

      %%% PRIMAL UPDATE
      % update S
      w = AT(r);
      divp = dxm(p(:,:,:,1)) + dym(p(:,:,:,2)) + dtm(p(:,:,:,3),mu);
      S = S - tau*(w + divp);

      % update V (TGV component)
      divq1 = dxp(q(:,:,:,1)) + dyp(q(:,:,:,4)) + dtp(q(:,:,:,5),mu);
      divq2 = dxp(q(:,:,:,4)) + dyp(q(:,:,:,2)) + dtp(q(:,:,:,6),mu);
      divq3 = dxp(q(:,:,:,5)) + dyp(q(:,:,:,6)) + dtp(q(:,:,:,3),mu);

      V(:,:,:,1) = V(:,:,:,1) - tau*(divq1 - p(:,:,:,1));
      V(:,:,:,2) = V(:,:,:,2) - tau*(divq2 - p(:,:,:,2));
      V(:,:,:,3) = V(:,:,:,3) - tau*(divq3 - p(:,:,:,3));

      % update L (low rank component)
      L = L - tau*w;
      
      L = reshape(L,[m*n,d]);
      [Lu,Ls,Lv] = svd(L,'econ');
      Ls = diag(max(diag(Ls)-beta,0));
      L = Lu*Ls*Lv';
      L = reshape(L,[m,n,d]);

      %%% INTERMEDIATE UPDATE
      Sbar = 2*S - Sold;
      Vbar = 2*V - Vold;
      Lbar = 2*L - Lold;

      %%% PRINT OUT RESULTS
      upd(i);

    end
  end
