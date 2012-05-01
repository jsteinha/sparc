%% generate the data set
N = 3;
Ps = 1.2*randn(N,2,2);
us = 5*randn(2,N);
ps = rand(N,1).*rand(N,1);
T = 5000;
v = zeros(2,N);
for t=1:T
    k = N;
    for i=1:N
        if rand < ps(i)
            k = i;
            break;
        end
    end
    v(:,t) = us(:,k)+reshape(Ps(k,:,:),[2 2])*randn(2,1);
end
clf;
plot(v(1,:),v(2,:),'.','Markersize',3)
%%
Data = v;
Data_size = size(Data,2);
Num_clusters = 3;
Dimension = 2;
init_pts = 5*randn(Dimension,Num_clusters);
zs = zeros(Data_size,1);
for i=1:Data_size
    dists = zeros(Num_clusters,1);
    for c=1:Num_clusters
        dists(c) = norm(Data(:,i)-init_pts(:,c));
    end
    [~,zs(i)] = min(dists);
end
%[~,zs] = randi(Num_clusters,Data_size,1); % cluster assignments, initialized randomly
mus = zeros(Dimension,Num_clusters); % the means for each cluster
for c=1:Num_clusters
    mus(:,c) = sum(Data(:,logical(zs==c)),2);
end
Js = zeros(Dimension,Dimension,Num_clusters); % the matrix hyperparameters for each cluster
for i=1:Data_size
    Js(:,:,zs(i)) = Js(:,:,zs(i)) + Data(:,i)*Data(:,i)';
end
counts = zeros(Num_clusters,1); % the count of the number of data assigned to each cluster
for c=1:Num_clusters
    counts(c) = sum(zs==c);
end

Num_steps = 300000;

global conc u0 s0 dof dof2;
conc = 1.0; % concentration parameter for the Dirichlet prior
u0 = zeros(Dimension,1);
s0 = eye(Dimension);
dof = 2.1; % degree-of-freedom parameter for the Wishart prior
dof2 = 0.1;

for t=1:Num_steps
    % MCMC step: pick a data point at random and propose to move it to a new
    % cluster, then update the cluster parameters
    index = randi(Data_size);
    datum = Data(:,index);
    cur_cluster = zs(index);
    prop_cluster = randi(Num_clusters+1);
%    mean = @(c) mus(:,c)/(counts(c)+dof);
%    cov = @(c) (Js(:,:,c)-mus(:,c)*mus(:,c)'/(counts(c)+1e-12) + ...
%        dof*eye(Dimension))/(counts(c)+dof);
%    logscorer = @(c) -0.5*(datum-mean(c))'*cov(c)*(datum-mean(c)) + log(counts(c)+conc) + 0.5*log(det(cov(c)/(2*pi)));
    logscorer = @(c) logpdf(datum,mus(:,c),Js(:,:,c),counts(c));
    cur_cluster_logscore = logscorer(cur_cluster);
    prop_cluster_logscore = logscorer(prop_cluster);
    if log(rand()) < prop_cluster_logscore - cur_cluster_logscore
%        disp('Accepted')
        zs(index) = prop_cluster;
        counts(cur_cluster) = counts(cur_cluster) - 1;
        mus(:,cur_cluster) = mus(:,cur_cluster) - datum;
        Js(:,:,cur_cluster) = Js(:,:,cur_cluster) - datum*datum';
        counts(prop_cluster) = counts(prop_cluster) + 1;
        mus(:,prop_cluster) = mus(:,prop_cluster) + datum;
        Js(:,:,prop_cluster) = Js(:,:,prop_cluster) + datum*datum';
    else
%        disp('Rejected')
    end
    if mod(t,3000) == 0
        getmean=@(c)getmu(mus(:,c),counts(c));
        getcov=@(c)getsig(mus(:,c),Js(:,:,c),counts(c));
        mu1=getmean(1);
        mu2=getmean(2);
        mu3=getmean(3);
        fprintf(1,'trial %d, mu1=(%f,%f), mu2=(%f,%f), mu3=(%f,%f)\n',t,mu1(1),mu1(2),...
            mu2(1),mu2(2),mu3(1),mu3(2));
        clf;
        hold on;
        plot(Data(1,logical(zs==1)),Data(2,logical(zs==1)),'r.','Markersize',1);
        drawgaussian(mu1,getcov(1),'r');
        plot(Data(1,logical(zs==2)),Data(2,logical(zs==2)),'g.','Markersize',1);
        drawgaussian(mu2,getcov(2),'g');
        plot(Data(1,logical(zs==3)),Data(2,logical(zs==3)),'b.','Markersize',1);
        drawgaussian(mu3,getcov(3),'b');
        drawnow;
%        pause;
    end
end
%%
N = 10000;
S = zeros(2,N);
S(:,1) = randn(2,1);
for i=2:N
    S(:,i) = [[0.995 -0.05];[0.05 0.995]]*S(:,i-1);
end
plot(S(1,:),S(2,:));
%%
T = S+1.0*randn(2,N);
plot(T(1,:),T(2,:))
%%
R = T;
