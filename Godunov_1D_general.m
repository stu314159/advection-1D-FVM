function q_out = Godunov_1D_general(q,lamda,phi)


N = length(q);
q_out = zeros(N,1);
% compute thetas based on q and lamda
theta = get_theta(q,lamda);
% compute updated q for each element.
ind = (1:N);
qp=q(circshift(ind,[0, -1]));
qm=q(circshift(ind,[0, 1]));
thetap=theta(circshift(ind,[0 -1]));

if (lamda>0)
    for iel = 1:(N-1)
        q_out(iel)=q(iel)-lamda*(q(iel)-qm(iel))-...
            ((lamda*(1-lamda))/2)*(phi(thetap(iel))*(qp(iel)-q(iel))-...
            phi(theta(iel))*(q(iel)-qm(iel)));
    end
else
    for iel=1:(N-1)
         q_out(iel)=q(iel)-lamda*(qp(iel)-q(iel))+...
            ((lamda*(1+lamda))/2)*(phi(thetap(iel))*(qp(iel)-q(iel))-...
            phi(theta(iel))*(q(iel)-qm(iel)));
    end
end