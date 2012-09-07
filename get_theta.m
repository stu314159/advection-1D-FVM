function theta = get_theta(q,lamda)
N = length(q);

theta = zeros(N,1);

ind = (1:N);

qp=q(circshift(ind,[0 -1]));
qm=q(circshift(ind,[0 1]));
qm2=q(circshift(ind,[0 2]));

if lamda>0
    for i = 1:(N)
        if((q(i)-qm(i))==0)
            theta(i)=0;
        else
        theta(i)=(qm(i)-qm2(i))/(q(i)-qm(i));
        end
    end
else
    for i = 1:(N)
        if((q(i)-qm(i))==0)
            theta(i)=0;
        else
        theta(i)=(qp(i)-q(i))/(q(i)-qm(i));
        end
    end
end