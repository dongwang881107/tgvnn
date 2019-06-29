function [dt] = dtp(u,mu)
    dt = cat(3,u(:,:,2:end),u(:,:,end)) - u;
    dt = dt/mu;
end
