function [dt] = dtm(u,mu)
    [M,N,~] = size(u);
    dt = cat(3,u(:,:,1:end-1),zeros(M,N,1)) - cat(3,zeros(M,N,1),u(:,:,1:end-1))/mu;
end