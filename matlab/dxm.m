function [dx] = dxm(u)
% Implementation of finite differences
% 
% (c) 31.8.2010
% Florian Knoll (florian.knoll@tugraz.at)
% Kristian Bredies (kristian.bredies@uni-graz.at
% Thomas Pock (pock@icg.tugraz.at)
% Rudolf Stollberger (rudolf.stollberger@tugraz.at)
%
% If you consider this code to be useful for your research, please cite:
% Knoll, F.; Bredies, K.; Pock, T.; Stollberger, R.: Second Order Total
% Generalized Variation (TGV) for MRI: Magnetic Resonance in Medicine, 
% to appear (2010)

[M N P] = size(u);
dx = cat(2,u(:,1:end-1,:),zeros(M,1,P)) - cat(2,zeros(M,1,P),u(:,1:end-1,:));
