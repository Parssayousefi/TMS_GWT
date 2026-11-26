param(
    [string]$Subject = '06',
    [ValidateRange(1,22)][int]$NProcs = 10,
    [ValidateRange(1,10)][int]$OmpThreads = 4,
    [ValidateRange(4096,131072)][int]$MemoryMB = 65536
)
$ErrorActionPreference = 'Stop'
if ($Subject -notin @('01','03','04','05','06','07','08','09','10','11','12','13','14','15','16')) { throw 'Only included Aachen subjects are authorized' }
$taskRoot = [IO.Path]::GetFullPath('D:\projects\TMS_GWT')
$taskOutput = [IO.Path]::GetFullPath((Join-Path $taskRoot 'data\aachen_itbs\derivatives\fmriprep_skullstrip_skip_20261125'))
$taskWork = [IO.Path]::GetFullPath((Join-Path $taskRoot ('tmp\fmriprep_aachen_skullstrip_skip_20261125\sub-' + $Subject)))
foreach ($taskPath in @($taskOutput,$taskWork)) {
    if (-not $taskPath.StartsWith($taskRoot + '\',[StringComparison]::OrdinalIgnoreCase)) { throw 'Path outside project' }
    [void](New-Item -ItemType Directory -Path $taskPath -Force)
}
$taskContainer = 'aachen_skip_' + $Subject + '_20261125'
$taskBids = (Join-Path $taskRoot 'data\aachen_itbs') -replace '\\','/'
$taskOut = $taskOutput -replace '\\','/'
$taskTmp = $taskWork -replace '\\','/'
$taskLicense = (Join-Path $taskRoot 'license.txt') -replace '\\','/'
$taskTemplate = (Join-Path $taskRoot 'tmp\templateflow') -replace '\\','/'
if (-not (Test-Path -LiteralPath $taskTemplate)) { throw 'Expected prepopulated TemplateFlow cache missing' }
if (Test-Path -LiteralPath (Join-Path $taskOutput ('sub-' + $Subject + '.html'))) { throw 'Pilot already has a report; inspect it rather than overwriting' }
$taskExisting = docker ps -a --format '{{.Names}}'
if ($taskContainer -in $taskExisting) { throw 'Pilot container already exists; inspect it instead of replacing' }
docker run --rm -v "${taskOut}:/out" --entrypoint /bin/bash nipreps/fmriprep:23.2.1 -c 'mkdir -p /out/sourcedata/freesurfer && cp -rn /opt/freesurfer/subjects/fsaverage /out/sourcedata/freesurfer/'
if ($LASTEXITCODE -ne 0) { throw 'fsaverage preparation failed' }
$taskDockerArgs = @('run','--detach','--name',$taskContainer,
    '-v',"${taskBids}:/data:ro",'-v',"${taskOut}:/out",'-v',"${taskTmp}:/work",
    '-v',"${taskLicense}:/opt/freesurfer/license.txt:ro",'-v',"${taskTemplate}:/templateflow",
    '-e','TEMPLATEFLOW_HOME=/templateflow','nipreps/fmriprep:23.2.1','/data','/out','participant',
    '--participant-label',$Subject,'-w','/work','--fs-license-file','/opt/freesurfer/license.txt',
    '--output-spaces','MNI152NLin2009cAsym:res-2','--nprocs',$NProcs.ToString(),'--omp-nthreads',$OmpThreads.ToString(),'--mem',$MemoryMB.ToString(),
    '--skip_bids_validation','--fs-no-reconall','--skull-strip-t1w','skip','--stop-on-first-crash','--notrack')
$taskId = docker @taskDockerArgs
if ($LASTEXITCODE -ne 0) { throw 'Failed to launch pilot' }
$taskRecord = [ordered]@{subject=$Subject;container=$taskContainer;container_id=$taskId;created=(Get-Date).ToUniversalTime().ToString('o');output=$taskOutput;work=$taskWork;docker_arguments=$taskDockerArgs;status='started';note='Clean isolated rerun. Original derivatives and downstream results retained. Validate anatomy and target alignment before promotion.'}
$taskRecord | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $taskRoot ('PROVENANCE\aachen_skullstrip_pilot_' + $Subject + '.json')) -Encoding utf8
Write-Output ('Started ' + $taskContainer + ' -> ' + $taskOutput)
