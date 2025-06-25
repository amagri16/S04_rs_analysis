Some toolbox needs extra setup which is a bit complicated:
1/f spectral parameterization: on code ft extract\Feature evaluation functions\Spatial_feature_calculation.m line 60. A python lib need to be installed and path of Python need to be added. See https://github.com/fooof-tools/fooof_mat
Complexity and entropy rate: on code ft extract\Feature evaluation functions\Spatial_feature_calculation.m line 61. A file need to be compiled. See https://github.com/pmediano/EntRate
Bispectrum: on code ft extract\Feature evaluation functions\Spatial_feature_calculation.m line 62. A MATLAB toolbox is needed to be installed. Link: https://uk.mathworks.com/matlabcentral/fileexchange/3013-hosa-higher-order-spectral-analysis-toolbox
Transfer entropy/Mutual information will throw error if on windows/macOS: This issue is fixable on Windows by install the visual studio C++ AIO software. We have no idea how to fix on macOS. Linux should be fine.
