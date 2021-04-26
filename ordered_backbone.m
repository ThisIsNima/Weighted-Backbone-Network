%%%%%%%%%%%% upload files

dir='/' % add the data directory
cd dir
files = dir('*.mat');
for i=1:length(files)
data1= load(files(i).name); % eval(['load ' files(i).name ]);
healthy_41{i}=data1.TSData.Region_41.allVoxts;
end
  
%%%%%%%%%%%% create the dynamic network connectivity with given temporal window size

Dynamic_C_Adjacencies=ST_Preparation(w_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%START


clear Backbones

for N=1:length(files)


clear arrz
clear arrz1

for i=1:length(Dynamic_C_Adjacencies{1, 1}{1, 1}) %number of voxels of the region
for j=1:length(Dynamic_C_Adjacencies{1, 1}{1, 1})

if i~=j

for t=1:length(Dynamic_C_Adjacencies{1, 1}) %number of temporal segments
arrz(t)=Adjacencies_Healthy_41_20_5{1, N}{t, 1}(i,j);
end

normA = arrz - min(arrz(:));
arrz1 = normA ./ max(normA(:));

for t=1:20
Adjacencies_Healthy_41_normalized{1, N}{t, 1}(i,j)=arrz1(t);
end

end


end
end

data=Adjacencies_Healthy_41_normalized{1, N};
for i=1:64
for t=1:20
data{t,1}(i,i)=0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Mean
clear arr1
clear mean_vals_41
for i=1:64
for j=1:64
for t=1:20
arr1(t)=data{t, 1}(i,j);
end


mean_vals_41(i,j)=mean(arr1);
end
end
for i=1:64
mean_vals_41(i,i)=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%initial values:

clear arr1
clear arr2
clear arr3
clear init_nominator
clear init_denominator
clear sum_rows
clear total_sum

%%%%%%%nominator
for i=1:64
for t=1:20
arr1(t)=sum(data{t, 1}(i,:))-1;
end


init_nominator(i)=mean(arr1);
end

%%%%%%%denominator
for i=1:64
for t=1:20
arr2(t)=sum(data{t, 1}(i,:))-1;
end
init_nominator(i)=mean(arr2);
end
for i=1:64
for t=1:20
arr3(t)=sum(data{t, 1}(i,:))-1;
end
sum_rows(i)=sum(arr3);
end
total_sum=sum(sum_rows);
init_denominator=sqrt(total_sum);
init_vals=init_nominator/init_denominator;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MLE mean

fun = @(x) optimize_backbone(x,mean_vals_41);
    [estimated_mean,~,exitf,~] = fsolve(fun,init_vals');

%'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%std vals from estimatd means

clear ar11
clear std_vals
for i=1:64
for j=1:64

for t=1:20
arr11(t)=(data{t, 1}(i,j)-(estimated_mean(i)*estimated_mean(j)))^2 ;
end

std_vals(i,j)=sqrt(sum(arr11)/20);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MLE std
fun = @(x) optimize_backbone(x,std_vals);
    [estimated_std,~,exitf,~] = fsolve(fun,init_vals');

% '
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Statistical test

% MLE_mean_41{N, 1}=estimated_mean;
% MLE_std_41{N, 1}=estimated_std;

% clear arr1
% clear mean_vals_41
% for i=1:64
% for j=1:64
% for t=1:20
% arr1(t)=data{t, 1}(i,j);
% end
% sum_weights(i,j)=sum(arr1);
% end
% end
% for i=1:64
% sum_weights(i,i)=0;
% end
% sum_weight{N, 1}=sum_weights;


% for k=1:14
% for i=1:64
% sum_weight_nodes(i,k)= sum(sum_weight{k, 1}(i,:));
% end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




for T=1:20

mat_41=data{T,1};

for i=1:64
mat_41(i,i)=0;
end



k=0.8;

z_t_41=zeros(64,64);
for i=1:64
for j=1:64

m_41=estimated_mean(i)*estimated_mean(j);
s_41=estimated_std(i)*estimated_std(j);
x_41=norminv(k,m_41,s_41);

thresholds(i,j)=x_41;
if mat_41(i,j)>x_41
z_t_41(i,j)=1;
else
z_t_41(i,j)=0;
end

end
end


count_ones_41(T)=sum(z_t_41(:) == 1);
matrices_41{T, 1}=z_t_41;

end % for T=1:20 loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%extract the backbone network

clear count_ones

disp('extracting backbone...')



for i=1:64
for j=1:64

for T=1:20
arr4(T)= matrices_41{T, 1}(i,j);
end

count_ones(i,j)=sum(arr4(:) == 1);
end
end





Backbones{N, 1}=count_ones;



%clear healthy_41_Final_Backbone

disp('final backbone calculation...')



for i=1:64
for j=1:64

if Backbones{N,1}(i,j)>=10
healthy_41_Final_Backbone_20_95{N,1}(i,j)=1;
else
healthy_41_Final_Backbone_20_95{N,1}(i,j)=0;
end

end
end



disp('end of loop for subject')
end
