param(
    [ValidateSet('s001','s022')][string]$Subject = 's001',
    [ValidateRange(1,20)][int]$NProcs = 12,
    [ValidateRange(1,8)][int]$OmpThreads = 4,
    [ValidateRange(8192,98304)][int]$MemoryMB = 49152,
    [switch]$PrepareOnly
)
$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$taskRaw = Join-Path $taskRoot 'data\ds001927'
$taskOutput = Join-Path $taskRaw 'derivatives\fmriprep_23_2_1'
$taskWork = Join-Path $taskRoot ('tmp\fmriprep_ds001927_20261126\sub-' + $Subject)
$taskProvenance = Join-Path $taskRoot 'PROVENANCE\ds001927_preprocessing_20261126'
$taskSettings = Join-Path $taskRoot 'reruns\ds001927_20261126\settings'
$taskReceipt = Join-Path $taskProvenance ('pilot_' + $Subject + '.json')
if (Test-Path -LiteralPath $taskReceipt) {
    $taskPrevious = Get-Content -LiteralPath $taskReceipt -Raw | ConvertFrom-Json
    if ($PrepareOnly -or $taskPrevious.status -ne 'prepared') {
        throw ('Existing receipt preserved (' + $taskPrevious.status + '); inspect it rather than repeat preparation or launch')
    }
}
foreach ($taskPath in @($taskOutput,$taskWork,$taskProvenance,$taskSettings)) {
    $taskResolved = [IO.Path]::GetFullPath($taskPath)
    if (-not $taskResolved.StartsWith($taskRoot + '\',[StringComparison]::OrdinalIgnoreCase)) { throw 'Path outside project' }
    [void](New-Item -ItemType Directory -Path $taskResolved -Force)
}
$taskValidation = Get-Content -LiteralPath (Join-Path $taskRoot 'output\ds001927_input_audit_20261126\bids_validation.json') -Raw | ConvertFrom-Json
if (@($taskValidation.issues.errors).Count -gt 0) { throw 'Resolve BIDS errors before launching' }
$taskAudit = Get-Content -LiteralPath (Join-Path $taskRoot 'output\ds001927_input_audit_20261126\summary.json') -Raw | ConvertFrom-Json
if ($taskAudit.structural_images -ne 69 -or $taskAudit.bold_runs -ne 138 -or -not $taskAudit.all_inspected_voxels_finite) { throw 'Input audit incomplete' }
$taskAnatSessions = @('ctr','fro','occ')
if ($Subject -eq 's022') { $taskAnatSessions = @('fro','occ') }
$taskFilter = [ordered]@{
    t1w = [ordered]@{datatype='anat';suffix='T1w';session=$taskAnatSessions}
    bold = [ordered]@{datatype='func';suffix='bold';session=@('ctr','fro','occ');task=@('restpre','restpost')}
}
$taskFilterFile = Join-Path $taskSettings ('sub-' + $Subject + '_bids_filter.json')
$taskFilterJson = $taskFilter | ConvertTo-Json -Depth 5
if (Test-Path -LiteralPath $taskFilterFile) {
    $taskOldFilter = Get-Content -LiteralPath $taskFilterFile -Raw
    if ($taskOldFilter.Trim() -ne $taskFilterJson.Trim()) { throw 'Existing input filter differs; preserve and review it' }
} else {
    $taskFilterJson | Set-Content -LiteralPath $taskFilterFile -Encoding utf8
}
$taskContainer = 'ds001927_pilot_' + $Subject + '_20261126'
$taskPosixRaw = $taskRaw -replace '\\','/'
$taskPosixOut = $taskOutput -replace '\\','/'
$taskPosixWork = $taskWork -replace '\\','/'
$taskPosixSettings = $taskSettings -replace '\\','/'
$taskLicense = (Join-Path $taskRoot 'license.txt') -replace '\\','/'
$taskTemplate = (Join-Path $taskRoot 'tmp\templateflow') -replace '\\','/'
if (-not (Test-Path -LiteralPath $taskTemplate)) { throw 'Prepopulated TemplateFlow cache missing' }
if (Test-Path -LiteralPath (Join-Path $taskOutput ('sub-' + $Subject + '.html'))) { throw 'An existing report requires review; do not overwrite or treat it as complete' }
$taskArgs = @('run','--detach','--name',$taskContainer,
    '-v',"${taskPosixRaw}:/data:ro",'-v',"${taskPosixOut}:/out",'-v',"${taskPosixWork}:/work",
    '-v',"${taskPosixSettings}:/settings:ro",'-v',"${taskLicense}:/opt/freesurfer/license.txt:ro",
    '-v',"${taskTemplate}:/templateflow",'-e','TEMPLATEFLOW_HOME=/templateflow',
    'nipreps/fmriprep:23.2.1','/data','/out','participant','--participant-label',$Subject,
    '--bids-filter-file',('/settings/sub-' + $Subject + '_bids_filter.json'),
    '-w','/work','--fs-license-file','/opt/freesurfer/license.txt',
    '--output-spaces','MNI152NLin2009cAsym:res-2','--nprocs',$NProcs.ToString(),
    '--omp-nthreads',$OmpThreads.ToString(),'--mem',$MemoryMB.ToString(),
    '--fs-no-reconall','--skull-strip-t1w','force','--longitudinal',
    '--ignore','fieldmaps','slicetiming','--random-seed','20261126',
    '--stop-on-first-crash','--notrack')
$taskRecord = [ordered]@{
    subject=$Subject; container=$taskContainer; created_utc=(Get-Date).ToUniversalTime().ToString('o')
    output=$taskOutput; work=$taskWork; docker_arguments=$taskArgs
    anatomical_sessions=$taskAnatSessions; expected_bold_runs=6; expected_volumes_per_run=600; tr_seconds=1.25
    input_decisions='Defaced, skull-on T1: extract once. Unbiased shared anatomical reference from usable T1s. Exclude only sub-s022 ses-ctr damaged T1. All six BOLD runs retained for QC. No supplied SliceTiming, PE direction or fieldmaps; no invented timing or distortion fields.'
    status='prepared'; no_downstream_analysis_started=$true
}
if ($PrepareOnly) {
    $taskRecord | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $taskReceipt -Encoding utf8
    Write-Output ('Prepared, not launched: ' + $taskContainer)
    exit 0
}
$taskDrive = [IO.DriveInfo]::new([IO.Path]::GetPathRoot($taskRoot))
if ($taskDrive.AvailableFreeSpace -lt 70GB) { throw 'Pilot requires at least 70 GiB free; do not remove other project data automatically' }
$taskExisting = docker ps -a --format '{{.Names}}'
if ($LASTEXITCODE -ne 0) { throw 'Cannot inspect Docker containers' }
if ($taskContainer -in $taskExisting) { throw 'Pilot container exists; inspect its status rather than replacing it' }
$taskActive = @(docker ps --format '{{.Names}}' | Where-Object { $_ -like 'ds001927_pilot_*' })
if ($LASTEXITCODE -ne 0) { throw 'Cannot inspect active Docker containers' }
if ($taskActive.Count -gt 0) { throw 'Another ds001927 pilot is active; wait for its disk and quality review' }
# fMRIPrep 23.2.1 still expects fsaverage when surface reconstruction is disabled.
docker run --rm -v "${taskPosixOut}:/out" --entrypoint /bin/bash nipreps/fmriprep:23.2.1 -c 'mkdir -p /out/sourcedata/freesurfer && cp -rn /opt/freesurfer/subjects/fsaverage /out/sourcedata/freesurfer/'
if ($LASTEXITCODE -ne 0) { throw 'fsaverage preparation failed' }
$taskImageId = docker image inspect nipreps/fmriprep:23.2.1 --format '{{.Id}}'
if ($LASTEXITCODE -ne 0) { throw 'Cannot identify installed image' }
$taskId = docker @taskArgs
if ($LASTEXITCODE -ne 0) { throw 'Pilot launch failed' }
$taskRecord['container_id'] = $taskId
$taskRecord['image_id'] = $taskImageId
$taskRecord['status'] = 'running; six-run completion and anatomical release pending'
$taskRecord | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $taskReceipt -Encoding utf8
Write-Output ('Started ' + $taskContainer + ' -> ' + $taskOutput)
