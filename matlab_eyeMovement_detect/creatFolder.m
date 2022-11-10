SavePath='D:\eye RMG\data\';
FolderName='EOG';

creatFolder1=[SavePath,'\',FolderName];

creatFolder2=[SavePath,'\',FolderName,'\','feat_all'];
creatFolder3=[SavePath,'\',FolderName,'\','feature'];
creatFolder4=[SavePath,'\',FolderName,'\','fig_case'];



status = mkdir(creatFolder2);
status = mkdir(creatFolder3);
status = mkdir(creatFolder4);

