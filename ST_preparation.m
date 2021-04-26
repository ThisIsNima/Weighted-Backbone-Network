

for f=1:15

x_raw_41=healthy_41{1, f};
k=1;



for i=1:5:100
   
x_41= x_raw_41(:,i:i+19);
a_41(k).k=x_41;
k=k+1;
end


clear corr_41
clear Adj_3s
for i=1:20
z_41= a_41(i);
%previously, for was up to (527) 1000
for j=1:64
for l=1:64
temp_corr_41=corrcoef(z_41.k(j,:),z_41.k(l,:));
corr_41(j,l)=temp_corr_41(1,2);
end
end
CORR_41(i).i=corr_41;
disp('i')
disp(i)
end


for i=1:20
Adj_41{i,1}=CORR_41(i).i;
end


Adjacencies_Healthy_41_20_5{f}=Adj_41;

end
