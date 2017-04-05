function CMF_install

% add paths
disp('Setting path...')
parentPath = fileparts(which(mfilename));
cd(fullfile(parentPath))
addpath(fullfile(parentPath, 'filters'));
addpath(fullfile(parentPath, 'auxiliary'));

% compile mex files
disp('Compile mex files...')
cd(fullfile(parentPath, 'filters'))
mex CXXFLAGS='-O2 -DNDEBUG' CMF_medfiltCirc2DMex.cpp      CMF_library.cpp
mex CXXFLAGS='-O2 -DNDEBUG' CMF_medfiltCirc2DQuantMex.cpp CMF_library.cpp
mex CXXFLAGS='-O2 -DNDEBUG' CMF_medfiltGeoRN2DMex.cpp     CMF_library.cpp
cd(fullfile(parentPath))

end