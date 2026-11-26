% Bounded verification only: no cohort driver or preprocessing is invoked.
root='D:\projects\TMS_GWT';
run=strtrim(fileread(fullfile(root,'output','targeted_information_flow','LATEST.txt')));
folder=fullfile(root,'output','targeted_information_flow',run);
addpath(fullfile(root,'pipeline','docs','luppi_supps','software','MMI_PhiID'));
d=load(fullfile(folder,'matlab_validation_input.mat'));
actual=zeros(size(d.expected));
for k=1:numel(d.vectors)
    x=d.vectors{k};
    a=PhiIDFull(x(1,:),x(2,:),x(3,:),x(4,:));
    for j=1:numel(d.atom_names)
        actual(k,j)=a.(d.atom_names{j});
    end
end
error_max=max(abs(actual-d.expected),[],'all');
assert(error_max<2e-8,'Python/MATLAB atom mismatch');
save(fullfile(folder,'matlab_validation_output.mat'),'actual','error_max');
fid=fopen(fullfile(folder,'matlab_verification.json'),'w');
fprintf(fid,'{"status":"passed","cases":%d,"atoms_per_case":16,"max_absolute_error":%.17g}',size(actual,1),error_max);
fclose(fid);
fprintf('MATLAB verification passed: %d cases, maximum error %.4g nats\n',size(actual,1),error_max);
