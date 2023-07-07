function [fmin,best,cg_curve,Time]=bGEBA(N, Max_iter, dim, data,trn,vald,TFid,classifierFhd)
tic;
A=0.50;
r=0.50;
r0=0.50;
alpha=0.9;
lamda=0.9;
Qmin=0;
Qmax=2;
N_iter=0;
Q=zeros(N,1);
v=zeros(N,dim);
S=zeros(N,dim);
cg_curve=zeros(1,Max_iter);
for i=1:N
    for j=1:dim 
        if rand<=0.5
            Sol(i,j)=0;
        else
            Sol(i,j)=1;
        end
    end
end

for i=1:N
    Fitness(i)=AccSz2(Sol(i,:), data,trn,vald,classifierFhd);
end
[fmin,I]=min(Fitness);
best=Sol(I,:);
while (N_iter<Max_iter)
    N_iter=N_iter+1;
    cg_curve(N_iter)=fmin;
    for i=1:N
        for j=1:dim
            Q(i)=Qmin+(Qmin-Qmax)*rand;
            v(i,j)=v(i,j)+(Sol(i,j)-best(j))*Q(i);
            S(i,j)=trnasferFun(S(i,j),v(i,j),TFid); 
            if rand>r
                S(i,j)=best(j);
            else
                mu(1,j)=(best(1,j)+Sol(i,j))/2;
                sigma(1,j)=abs(best(1,j)-Sol(i,j));
                temp=mu(1,j)+sigma(1,j)*randn;
                S(i,j)=trnasferFun(S(i,j),temp,TFid); 
            end
        end
        Fnew=AccSz2(S(i,:),data,trn,vald,classifierFhd);
        if Fnew<=Fitness(i)
            Sol(i,:)=S(i,:);
            Fitness(i)=Fnew;
        end
        A=alpha*A;
        r=r0*(1-exp(-1*lamda*N_iter));
        if Fnew<=fmin
            best=S(i,:);
            fmin=Fnew;
        end
    end
    sum_diameter = 0;
    i_rand=ceil(rand*N/3);
    while i_rand==1
        i_rand=ceil(rand*N/3);
    end
    habitat_center(1,:) = (best(1,:)+Sol(i_rand,:))./2;
    habitat_diameter(1,:) = (best(1,:)-Sol(i_rand,:)).^2;
    sum_diameter = sum_diameter+habitat_diameter(1,:);
    sqrt_diameter = sqrt(sum_diameter);
    N_reproduced = floor(N*0.95);
    for i=N_reproduced:N
        k=rand;
        if k>=0.45
            S(i,:)=habitat_center+rand*sqrt_diameter;
        else
            S(i,:)=k.*habitat_center;
        end
        for j=1:dim
            Sol(i,j)=trnasferFun(S(i,j),Sol(i,j),TFid); 
        end
        Fnew=AccSz2(Sol(i,:),data,trn,vald,classifierFhd);
        if (Fnew<=Fitness(i))
            Sol(i,:)=S(i,:);
            Fitness(i)=Fnew;
        end
        if Fnew<=fmin
            best=S(i,:);
            fmin=Fnew;
        end
    end
end
Time = toc;

