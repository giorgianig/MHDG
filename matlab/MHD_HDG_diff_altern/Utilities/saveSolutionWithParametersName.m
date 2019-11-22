% Save solution

saveSolutionName = ['Saves/u_' meshName '_Diffusion_' num2str(diff_n) '_tau' num2str(tau) '_tc' num2str(testcase.n)];
if driftvel
    saveSolutionName = [saveSolutionName '_Drift'];
end
if shockcapt
    saveSolutionName = [saveSolutionName '_SC' num2str(shockcapt)];
end
if localDiffPoints
   saveSolutionName = [saveSolutionName '_LocDiffPts' num2str(localDiffPoints)]; 
end
if limitRhoMin
   saveSolutionName = [saveSolutionName '_LimRhoMin' num2str(limitRhoMin)]; 
end
if useThreshold
   saveSolutionName = [saveSolutionName '_UseTh' num2str(useThreshold)]; 
end

if converged_case
    save([saveSolutionName '_WithCurvature.mat'],'u0','u_tilde','simulation_parameters')
else
    save([saveSolutionName '_WithCurvature_NOCONV' num2str(iStep) '.mat'],'u0','u_tilde','simulation_parameters')
end


