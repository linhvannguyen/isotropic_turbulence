nc1 = netcdf('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLRfeatures/SRrec/SR_patchesLRfeatures_patchsize08_all.nc','clobber');

nc1('Nt')=0; 
nc1('Ny')=Nh;
nc1('Nz')=Nh;
nc1{'X_HR_rec'}=ncfloat('Nt','Ny','Nz');

nc2 = netcdf('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLRfeatures/SRrec/SR_patchesLRfeatures_patchsize08_all_1.nc','r');
for t=1:200
    temp1=nc2{'X_HR_rec'}(t,:,:);
    % SAVE FILE
    nc1{'X_HR_rec'}(t,:,:)=temp1;       
end
close(nc2); 


nc3 = netcdf('/data/ISOTROPIC/DICTIONARYLEARNING/final/downsample4_3/joint_HRLRfeatures/SRrec/SR_patchesLRfeatures_patchsize08_all_2.nc','r');

for t=201:592
    temp1=nc3{'X_HR_rec'}(t-200,:,:);
    % SAVE FILE
    nc1{'X_HR_rec'}(t,:,:)=temp1;
       
end
close(nc3); 

close(nc1);