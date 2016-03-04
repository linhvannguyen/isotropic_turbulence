clear all; close all; clc;
Nh=384; % Number of grid-point in HR field
k_max=192;

%% Loop all file
for file_ids=20:56
    start=tic();
    filename=strcat('FIELD-',num2str(file_ids,'%.3d'),'.nc');
    cmdString = ['scp nguyen@lmlm6-3.univ-lille1.fr:/data4/data/',filename,' /data/ISOTROPIC/data/'];
    fprintf(strcat('Loading:',filename,'\n'));
    [status, ~] = unix(cmdString);
    if status==0
        fprintf([filename,' downloaded in ',num2str(toc(start),'%.2f'),' seconds! \n']);

        nc=netcdf(['/data/ISOTROPIC/data/', filename],'r'); 
        velocity_x_resolved=nc{'velocity_x'}(:,:,:);
        velocity_y_resolved=nc{'velocity_y'}(:,:,:);
        velocity_z_resolved=nc{'velocity_z'}(:,:,:);
        close(nc)

        cmdString = ['rm ',strcat('/data/ISOTROPIC/data/',filename)];
        fprintf(strcat('Delete:',filename,'\n'));
        [status, ~] = unix(cmdString);
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%% RATIO 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        start=tic();
        % cutoff
        ratio=2;
        k_cutoff=k_max/ratio;
        Nl=2*k_cutoff; % Number of grid-point in LR field

        %% U component
        F = fftn(velocity_x_resolved);

        % Sample: remove high frequency
        ids=1:2*k_max;
        ids([k_cutoff+1:k_max k_max+1:2*k_max-k_cutoff])=[];
        F=F(ids,ids,ids);

        % Remove 3D wave number > k_cutoff
        k_1D_LR=[0:k_cutoff -k_cutoff+1:1:-1];
        
        for m=1:2*k_cutoff
            for n=1:2*k_cutoff
                for p=1:2*k_cutoff
                    if round(sqrt(k_1D_LR(m)^2+k_1D_LR(n)^2+k_1D_LR(p)^2))>k_cutoff
                        F(m,n,p)=0;
                    end
                end
            end
        end

        % ifftn back to physical space (normalize by size)
        velocity_x_downsampled=(Nl/Nh)^3*ifftn(F);
        clearvars F
        
        %% V component
        F = fftn(velocity_y_resolved);

        % remove high frequency
        ids=1:2*k_max;
        ids([k_cutoff+1:k_max k_max+1:2*k_max-k_cutoff])=[];
        F=F(ids,ids,ids);

        % Remove 3D wave number > k_cutoff
        k_1D_LR=[0:k_cutoff -k_cutoff+1:1:-1];
        for m=1:2*k_cutoff
            for n=1:2*k_cutoff
                for p=1:2*k_cutoff
                    if round(sqrt(k_1D_LR(m)^2+k_1D_LR(n)^2+k_1D_LR(p)^2))>k_cutoff
                        F(m,n,p)=0;
                    end
                end
            end
        end

        % ifftn back to physical space (normalize by size)
        velocity_y_downsampled=(Nl/Nh)^3*ifftn(F);
        clearvars F    
        
        
        %% W component
        F = fftn(velocity_z_resolved);

        % remove high frequency
        ids=1:2*k_max;
        ids([k_cutoff+1:k_max k_max+1:2*k_max-k_cutoff])=[];
        F=F(ids,ids,ids);

        % Remove 3D wave number > k_cutoff
        k_1D_LR=[0:k_cutoff -k_cutoff+1:1:-1];
        for m=1:2*k_cutoff
            for n=1:2*k_cutoff
                for p=1:2*k_cutoff
                    if round(sqrt(k_1D_LR(m)^2+k_1D_LR(n)^2+k_1D_LR(p)^2))>k_cutoff
                        F(m,n,p)=0;
                    end
                end
            end
        end

        % ifftn back to physical space (normalize by size)
        velocity_z_downsampled=(Nl/Nh)^3*ifftn(F);
        clearvars F   
        
        %% SAVE TO FILE
        
        filename=strcat('/data/ISOTROPIC/data/FIELD-',num2str(file_ids,'%.3d'),'_downsampled2.nc');
        nc = netcdf(filename,'clobber');
        nc('Nx')=Nl; 
        nc('Ny')=Nl;
        nc('Nz')=Nl;
        nc{'velocity_x'}=ncfloat('Nx','Ny','Nz');
        nc{'velocity_x'}(:,:,:)=real(velocity_x_downsampled); 
        nc{'velocity_y'}=ncfloat('Nx','Ny','Nz');
        nc{'velocity_y'}(:,:,:)=real(velocity_y_downsampled); 
        nc{'velocity_z'}=ncfloat('Nx','Ny','Nz');
        nc{'velocity_z'}(:,:,:)=real(velocity_z_downsampled); 
        close(nc)        
        
        fprintf(strcat('Downsample by 2 in:',num2str(toc(start),'%.2f'),' seconds! \n'));
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% RATIO 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        start=tic();
        % cutoff
        ratio=4;
        k_cutoff=k_max/ratio;
        Nl=2*k_cutoff; % Number of grid-point in LR field

        %% U component
        F = fftn(velocity_x_resolved);

        % Sample: remove high frequency
        ids=1:2*k_max;
        ids([k_cutoff+1:k_max k_max+1:2*k_max-k_cutoff])=[];
        F=F(ids,ids,ids);

        % Remove 3D wave number > k_cutoff
        k_1D_LR=[0:k_cutoff -k_cutoff+1:1:-1];
        
        for m=1:2*k_cutoff
            for n=1:2*k_cutoff
                for p=1:2*k_cutoff
                    if round(sqrt(k_1D_LR(m)^2+k_1D_LR(n)^2+k_1D_LR(p)^2))>k_cutoff
                        F(m,n,p)=0;
                    end
                end
            end
        end

        % ifftn back to physical space (normalize by size)
        velocity_x_downsampled=(Nl/Nh)^3*ifftn(F);
        clearvars F
        
        %% V component
        F = fftn(velocity_y_resolved);

        % remove high frequency
        ids=1:2*k_max;
        ids([k_cutoff+1:k_max k_max+1:2*k_max-k_cutoff])=[];
        F=F(ids,ids,ids);

        % Remove 3D wave number > k_cutoff
        k_1D_LR=[0:k_cutoff -k_cutoff+1:1:-1];
        for m=1:2*k_cutoff
            for n=1:2*k_cutoff
                for p=1:2*k_cutoff
                    if round(sqrt(k_1D_LR(m)^2+k_1D_LR(n)^2+k_1D_LR(p)^2))>k_cutoff
                        F(m,n,p)=0;
                    end
                end
            end
        end

        % ifftn back to physical space (normalize by size)
        velocity_y_downsampled=(Nl/Nh)^3*ifftn(F);
        clearvars F    
        
        
        %% W component
        F = fftn(velocity_z_resolved);

        % remove high frequency
        ids=1:2*k_max;
        ids([k_cutoff+1:k_max k_max+1:2*k_max-k_cutoff])=[];
        F=F(ids,ids,ids);

        % Remove 3D wave number > k_cutoff
        k_1D_LR=[0:k_cutoff -k_cutoff+1:1:-1];
        for m=1:2*k_cutoff
            for n=1:2*k_cutoff
                for p=1:2*k_cutoff
                    if round(sqrt(k_1D_LR(m)^2+k_1D_LR(n)^2+k_1D_LR(p)^2))>k_cutoff
                        F(m,n,p)=0;
                    end
                end
            end
        end

        % ifftn back to physical space (normalize by size)
        velocity_z_downsampled=(Nl/Nh)^3*ifftn(F);
        clearvars F   
        
        %% SAVE TO FILE
        
        filename=strcat('/data/ISOTROPIC/data/FIELD-',num2str(file_ids,'%.3d'),'_downsampled4.nc');
        nc = netcdf(filename,'clobber');
        nc('Nx')=Nl; 
        nc('Ny')=Nl;
        nc('Nz')=Nl;
        nc{'velocity_x'}=ncfloat('Nx','Ny','Nz');
        nc{'velocity_x'}(:,:,:)=real(velocity_x_downsampled); 
        nc{'velocity_y'}=ncfloat('Nx','Ny','Nz');
        nc{'velocity_y'}(:,:,:)=real(velocity_y_downsampled); 
        nc{'velocity_z'}=ncfloat('Nx','Ny','Nz');
        nc{'velocity_z'}(:,:,:)=real(velocity_z_downsampled); 
        close(nc)        
        
        fprintf(strcat('Downsample by 4 in:',num2str(toc(start),'%.2f'),' seconds! \n'));

        
        
        %% Clear vars
        clearvars velocity_x_resolved velocity_y_resolved velocity_z_resolved

    else
        fprintf(strcat(filename,' not exist! \n'));
    end
end

%% Plot
% fig1=figure(); imagesc(squeeze(velocity_x_resolved(:,:,1))); caxis([-4,4])
% fig2=figure(); imagesc(squeeze(real(velocity_x_filter(:,:,1)))); caxis([-4,4])









%%%%%%%%%%%%%%% COMBINE ALL DATA %%%%%%%%%%%%%%%%%%%%%%

%% Loop all file
Nh=96; % Number of grid-point in HR field
nc=netcdf('/data/ISOTROPIC/data/data_downsampled4.nc','clobber');
nc('Nt')=0; 
nc('Nx')=Nh;
nc('Ny')=Nh;
nc('Nz')=Nh;
nc{'velocity_x'}=ncfloat('Nt','Nx','Ny','Nz');
nc{'velocity_y'}=ncfloat('Nt','Nx','Ny','Nz');
nc{'velocity_z'}=ncfloat('Nt','Nx','Ny','Nz');
count=0;
for file_ids=20:56
    count=count+1;
    filename = strcat('/data/ISOTROPIC/data/FIELD-',num2str(file_ids,'%.3d'),'_downsampled4.nc');
    ncload(filename);
    nc{'velocity_x'}(count,:,:,:) = velocity_x;
    nc{'velocity_y'}(count,:,:,:) = velocity_y;
    nc{'velocity_z'}(count,:,:,:) = velocity_z;
end
close(nc);

%% Loop all file
Nh=192; % Number of grid-point in HR field

nc=netcdf('/data/ISOTROPIC/data/data_downsampled2.nc','clobber');
nc('Nt')=0; 
nc('Nx')=Nh;
nc('Ny')=Nh;
nc('Nz')=Nh;
nc{'velocity_x'}=ncfloat('Nt','Nx','Ny','Nz');
nc{'velocity_y'}=ncfloat('Nt','Nx','Ny','Nz');
nc{'velocity_z'}=ncfloat('Nt','Nx','Ny','Nz');
count=0;
for file_ids=20:56
    count=count+1;
    filename = strcat('/data/ISOTROPIC/data/FIELD-',num2str(file_ids,'%.3d'),'_downsampled2.nc');
    ncload(filename);
    nc{'velocity_x'}(count,:,:,:) = velocity_x;
    nc{'velocity_y'}(count,:,:,:) = velocity_y;
    nc{'velocity_z'}(count,:,:,:) = velocity_z;
end
close(nc);
