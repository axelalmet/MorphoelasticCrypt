function curveNew = InterpolateToNewMesh(meshNew, meshOld, curveOld)
% This function is designed to interpolate the curves to the meshes that
% are needed for bvp4c's multipoint functionality. For piecewise-continuous
% meshes, interp1/interparc don't quite work, so we have to split the mesh
% up and do interp1/interparc on each sub-interval. Note that this assumes
% (lazily) that we only have at most two sub-intervals, and that the divide
% occurs at x = 1.

for i = 1:(length(meshOld) - 1)
    if (meshOld(i + 1) == meshOld(i))
        contIndexOld = i;
        break
    end
end

for i = 1:(length(meshNew) - 1)
    if (meshNew(i + 1) == meshNew(i))
        contIndexNew = i;
        break
    end
end

curveInterpOne = interp1(meshOld(1:contIndexOld), curveOld(1:contIndexOld), ...
                    meshNew(1:contIndexNew));
                              
curveInterpTwo = interp1(meshOld((contIndexOld + 1):end), curveOld((contIndexOld + 1):end), ...
                    meshOld((contIndexNew + 1):end)); 
         
curveNew = [curveInterpOne, curveInterpTwo];



