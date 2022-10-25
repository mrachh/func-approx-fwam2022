ichebfun = 0;
if(exist('chebfun','folder'))
    addpath('./chebfun')
    ichebfun = 1
end
fprintf('---------------------------------------------\n');
if(ichebfun)
    fprintf('chebfun: Installation successful\n');
else
    fprintf('chebfun: installation failure\n submodule not found\n'); 
end

