%efield = h5read('CylinderModel_6per_efield.h5');

hdf5file = 'CylinderModel_6per_efield.h5';
 
H5F.open(hdf5file, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
%
    % Specify the paths to the groups
    groupPath_xreal = '/Efield/MT/Tx05/Ex/Subgrid01/xreal';
   
    
    data_xreal = h5read(hdf5File,  groupPath_xreal );
 
 
    hdf5File_2 = "ReadWrite_test_Efield.h5";
 
    H5F.open('test_efield.h5', 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
 
    % Specify the paths to the groups
    groupPath_xreal_2 = '/Efield/MT/Tx01/Ex/Subgrid01/xreal';
  
    data_xreal_2 = h5read(hdf5File,  groupPath_xreal_2 );
    %