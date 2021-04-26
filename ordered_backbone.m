%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% upload files

dir='/' % enter the data directory (look at example data in the data directory)
cd dir
files = dir('*.mat');
for i=1:length(files)
data1= load(files(i).name); % eval(['load ' files(i).name ]);
healthy_41{i}=data1.TSData.Region_41.allVoxts;    %for region 41 (change the region number based on your interest)
end
  
%%%%%%%%%%%% create the dynamic network connectivity with given temporal window size

Dynamic_Adjacencies=ST_Preparation(w_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START


clear Backbones

for N=1:length(files)


clear arrz
clear arrz1
num_voxels=length(Dynamic_Adjacencies{1, 1}{1, 1});%number of voxels of the region
tau=length(Dynamic_Adjacencies{1, 1}); %number of temporal segments
for i=1:num_voxels
for j=1:num_voxels

if i~=j

for t=1:tau 
arrz(t)=Dynamic_Adjacencies{1, N}{t, 1}(i,j);
end

normA = arrz - min(arrz(:));
arrz1 = normA ./ max(normA(:));

for t=1:tau
Dynamic_Adjacencies_normalized{1, N}{t, 1}(i,j)=arrz1(t);
end

end


end
end

data=Dynamic_Adjacencies_normalized{1, N};
for i=1:length(Dynamic_C_Adjacencies{1, 1}{1, 1})
for t=1:tau
data{t,1}(i,i)=0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mean
clear arr1
clear mean_values
for i=1:num_voxels
for j=1:num_voxels
for t=1:tau
arr1(t)=data{t, 1}(i,j);
end


mean_values(i,j)=mean(arr1);
end
end
for i=1:num_voxels
mean_values(i,i)=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate Initial values:

clear arr1
clear arr2
clear arr3
clear init_nominator
clear init_denominator
clear sum_rows
clear total_sum

%%%%%%%%%%%%%% nominator
for i=1:num_voxels
for t=1:tau
arr1(t)=sum(data{t, 1}(i,:))-1;
end


init_nominator(i)=mean(arr1);
end

%%%%%%%%%%%%%% denominator
for i=1:num_voxels
for t=1:tau
arr2(t)=sum(data{t, 1}(i,:))-1;
end
init_nominator(i)=mean(arr2);
end
for i=1:num_voxels
for t=1:tau
arr3(t)=sum(data{t, 1}(i,:))-1;
end
sum_rows(i)=sum(arr3);
end
total_sum=sum(sum_rows);
init_denominator=sqrt(total_sum);
init_vals=init_nominator/init_denominator;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MLE mean optimization

fun = @(x) optimize_backbone(x,mean_values);
    [estimated_mean,~,exitf,~] = fsolve(fun,init_vals');

%'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% std vals from estimatd means

clear ar11
clear std_values
for i=1:num_voxels
for j=1:num_voxels

for t=1:tau
arr3(t)=(data{t, 1}(i,j)-(estimated_mean(i)*estimated_mean(j)))^2 ;
end

std_values(i,j)=sqrt(sum(arr3)/tau);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MLE std optimization
fun = @(x) optimize_backbone(x,std_values);
    [estimated_std,~,exitf,~] = fsolve(fun,init_vals');

% '
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistical test


for T=1:tau

mat_41=data{T,1};

for i=1:num_voxels
mat_41(i,i)=0;
end



k=0.8; %backbone threshold

z_t_41=zeros(num_voxels,num_voxels);
for i=1:num_voxels
for j=1:num_voxels

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extract the backbone network

clear count_ones

disp('extracting backbone...')



for i=1:num_voxels
for j=1:num_voxels
for T=1:tau
arr4(T)= matrices_41{T, 1}(i,j);
end
count_ones(i,j)=sum(arr4(:) == 1);
end
end
Backbones{N, 1}=count_ones;



%clear Final_Backbone

disp('final backbone calculation...')



for i=1:num_voxels
for j=1:num_voxels

if Backbones{N,1}(i,j)>=(tau/2)
Final_Backbone{N,1}(i,j)=1;
else
Final_Backbone{N,1}(i,j)=0;
end

end
end



disp('end of loop for subject')
end
