% Save solution

saveSolutionName = ['Saves/u_' meshName '_Dpe' num2str(diff_n) '_Dpai' num2str(diff_pari) '_Dpae' num2str(diff_pare)  '_Mref' num2str(Mref) '_tau' num2str(tau) '_tc' num2str(testcase.n)];
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
    save([saveSolutionName '.mat'],'u0','u_tilde','gradU','simulation_parameters')
else
    save([saveSolutionName '_NOCONV' num2str(iStep) '.mat'],'u0','u_tilde','gradU','simulation_parameters')
end


