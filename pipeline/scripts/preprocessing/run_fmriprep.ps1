<#
.SYNOPSIS
    Parallel fMRIPrep preprocessing via Docker containers.

.DESCRIPTION
    Dispatcher-worker script that runs fMRIPrep on multiple subjects in parallel.
    Each worker runs a separate Docker container processing one subject at a time.
    Automatically skips already-processed subjects (resume support).

.PARAMETER Dataset
    Dataset name from dataset_configs.yaml (e.g., hcp, berlin_ketamine).

.PARAMETER Parallel
    Number of simultaneous Docker containers. Default: 4.

.PARAMETER WslMemoryGB
    Total WSL memory in GB. Parsed from .wslconfig if not specified.

.PARAMETER DryRun
    Show plan without executing.

.PARAMETER NoResume
    Don't skip already-processed subjects.

.PARAMETER MonitorIntervalSec
    Progress check interval in seconds. Default: 120.

.EXAMPLE
    .\run_fmriprep.ps1 -Dataset hcp -Parallel 4 -DryRun
    .\run_fmriprep.ps1 -Dataset hcp -Parallel 4
    .\run_fmriprep.ps1 -Dataset berlin_ketamine -Parallel 2
#>

param(
    [Parameter(Mandatory=$true)]
    [string]$Dataset,

    [int]$Parallel = 4,

    [int]$WslMemoryGB = 0,

    [switch]$DryRun,

    [switch]$NoResume,

    [int]$MonitorIntervalSec = 120
)

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$Timestamp = Get-Date -Format "yyyyMMdd_HHmmss"

# ============================================================================
# Section A: Validation
# ============================================================================

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host " fMRIPrep Parallel Dispatcher" -ForegroundColor Cyan
Write-Host "========================================`n" -ForegroundColor Cyan

# Check Docker is running
Write-Host "[A] Validating environment..." -ForegroundColor Yellow
try {
    $dockerInfo = docker info 2>&1
    if ($LASTEXITCODE -ne 0) {
        Write-Host "ERROR: Docker is not running. Start Docker Desktop first." -ForegroundColor Red
        exit 1
    }
    Write-Host "  Docker: OK" -ForegroundColor Green
} catch {
    Write-Host "ERROR: Docker command not found. Is Docker Desktop installed?" -ForegroundColor Red
    exit 1
}

# Check fMRIPrep image exists
$fmriprepImage = "nipreps/fmriprep:23.2.1"
$imageCheck = docker images --format "{{.Repository}}:{{.Tag}}" | Where-Object { $_ -eq $fmriprepImage }
if (-not $imageCheck) {
    Write-Host "ERROR: fMRIPrep image '$fmriprepImage' not found. Pull it with:" -ForegroundColor Red
    Write-Host "  docker pull $fmriprepImage" -ForegroundColor Yellow
    exit 1
}
Write-Host "  fMRIPrep image: OK ($fmriprepImage)" -ForegroundColor Green

# Check disk space (need at least 50 GB free on D:)
$drive = Get-PSDrive D -ErrorAction SilentlyContinue
if ($drive) {
    $freeGB = [math]::Round($drive.Free / 1GB, 1)
    if ($freeGB -lt 50) {
        Write-Host "WARNING: Only ${freeGB} GB free on D:. Recommend at least 50 GB." -ForegroundColor Yellow
    } else {
        Write-Host "  Disk space: ${freeGB} GB free on D:" -ForegroundColor Green
    }
}

# ============================================================================
# Section B: Config loading
# ============================================================================

Write-Host "`n[B] Loading configuration for dataset '$Dataset'..." -ForegroundColor Yellow

$configLoaderPath = Join-Path $ScriptDir "fmriprep_config_loader.py"
if (-not (Test-Path $configLoaderPath)) {
    Write-Host "ERROR: Config loader not found at $configLoaderPath" -ForegroundColor Red
    exit 1
}

$configArgs = @("$configLoaderPath", "--dataset", $Dataset, "--output", "json")
if (-not $NoResume) {
    $configArgs += "--check-processed"
}

$configJson = python @configArgs 2>&1
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Config loader failed:" -ForegroundColor Red
    Write-Host $configJson -ForegroundColor Red
    exit 1
}

$config = $configJson | ConvertFrom-Json

# Extract paths
$BIDS_ROOT = $config.paths.bids_root
$OUTPUT_DIR = $config.paths.fmriprep_dir
$LICENSE_FILE = $config.paths.license_file
$WORK_DIR = $config.paths.work_dir
$LOG_DIR = $config.paths.log_dir
$TEMPLATEFLOW_DIR = $config.paths.templateflow_dir
$FMRIPREP_VERSION = $config.fmriprep.version
$OUTPUT_SPACES = $config.fmriprep.output_spaces
$EXTRA_ARGS = if ($config.fmriprep.extra_args) { $config.fmriprep.extra_args } else { "" }

# Determine subject list
if ($NoResume) {
    $subjects = @($config.subjects)
    $processedCount = 0
} else {
    $subjects = @($config.pending)
    $processedCount = $config.processed_count
}

$totalSubjects = @($config.subjects).Count

Write-Host "  Dataset: $($config.dataset_name) ($Dataset)" -ForegroundColor Green
Write-Host "  BIDS root: $BIDS_ROOT" -ForegroundColor Green
Write-Host "  Output dir: $OUTPUT_DIR" -ForegroundColor Green
Write-Host "  TemplateFlow: $TEMPLATEFLOW_DIR" -ForegroundColor Green
if ($EXTRA_ARGS) {
    Write-Host "  Extra args: $EXTRA_ARGS" -ForegroundColor Green
}

if (-not $NoResume -and $processedCount -gt 0) {
    Write-Host "  Already processed: $processedCount / $totalSubjects" -ForegroundColor Green
}
Write-Host "  Subjects to process: $($subjects.Count)" -ForegroundColor Green

# Validate paths exist
$pathErrors = @()
if (-not (Test-Path $BIDS_ROOT)) { $pathErrors += "BIDS root not found: $BIDS_ROOT" }
if (-not (Test-Path $LICENSE_FILE)) { $pathErrors += "FreeSurfer license not found: $LICENSE_FILE" }
if ($TEMPLATEFLOW_DIR -and -not (Test-Path $TEMPLATEFLOW_DIR)) {
    Write-Host "  WARNING: TemplateFlow dir not found: $TEMPLATEFLOW_DIR (will download templates)" -ForegroundColor Yellow
}

if ($pathErrors.Count -gt 0) {
    foreach ($err in $pathErrors) {
        Write-Host "ERROR: $err" -ForegroundColor Red
    }
    exit 1
}

if ($subjects.Count -eq 0) {
    Write-Host "`nAll subjects already processed! Nothing to do." -ForegroundColor Green
    exit 0
}

# ============================================================================
# Section C: Resource calculation
# ============================================================================

Write-Host "`n[C] Calculating resources..." -ForegroundColor Yellow

# Parse .wslconfig for memory if not specified
if ($WslMemoryGB -eq 0) {
    $wslConfigPath = "$env:USERPROFILE\.wslconfig"
    if (Test-Path $wslConfigPath) {
        $wslContent = Get-Content $wslConfigPath -Raw
        if ($wslContent -match 'memory\s*=\s*(\d+)GB') {
            $WslMemoryGB = [int]$Matches[1]
            Write-Host "  Parsed from .wslconfig: ${WslMemoryGB} GB" -ForegroundColor Green
        }
    }
    if ($WslMemoryGB -eq 0) {
        # Default: assume half of system RAM
        $totalPhysicalMem = [math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB)
        $WslMemoryGB = [math]::Floor($totalPhysicalMem / 2)
        Write-Host "  No .wslconfig memory setting found. Assuming ${WslMemoryGB} GB." -ForegroundColor Yellow
    }
}

# Parse CPU count from .wslconfig or use system count
$totalCPUs = 0
$wslConfigPath = "$env:USERPROFILE\.wslconfig"
if (Test-Path $wslConfigPath) {
    $wslContent = Get-Content $wslConfigPath -Raw
    if ($wslContent -match 'processors\s*=\s*(\d+)') {
        $totalCPUs = [int]$Matches[1]
    }
}
if ($totalCPUs -eq 0) {
    $totalCPUs = (Get-CimInstance Win32_ComputerSystem).NumberOfLogicalProcessors
}

# Reserve some resources for system overhead
$usableMemoryGB = [math]::Floor($WslMemoryGB * 0.9)  # 10% overhead for OS/Docker
$usableCPUs = $totalCPUs - 2  # Leave 2 cores for system

$cpusPerWorker = [math]::Max(1, [math]::Floor($usableCPUs / $Parallel))
$memGBPerWorker = [math]::Max(8, [math]::Floor($usableMemoryGB / $Parallel))

Write-Host "  WSL memory: ${WslMemoryGB} GB (usable: ${usableMemoryGB} GB)" -ForegroundColor Green
Write-Host "  Total CPUs: $totalCPUs (usable: $usableCPUs)" -ForegroundColor Green
Write-Host "  Workers: $Parallel" -ForegroundColor Green
Write-Host "  Per worker: $cpusPerWorker CPUs, $memGBPerWorker GB RAM" -ForegroundColor Green

# Sanity checks
if ($memGBPerWorker -lt 8) {
    Write-Host "WARNING: Less than 8 GB per worker. fMRIPrep may fail or be very slow." -ForegroundColor Yellow
    Write-Host "  Consider reducing -Parallel or increasing -WslMemoryGB." -ForegroundColor Yellow
}

# ============================================================================
# Section D: Queue split (round-robin)
# ============================================================================

Write-Host "`n[D] Splitting $($subjects.Count) subjects across $Parallel workers..." -ForegroundColor Yellow

$queues = @{}
for ($i = 0; $i -lt $Parallel; $i++) {
    $queues[$i] = @()
}

for ($i = 0; $i -lt $subjects.Count; $i++) {
    $workerIdx = $i % $Parallel
    $queues[$workerIdx] += $subjects[$i]
}

for ($i = 0; $i -lt $Parallel; $i++) {
    $queueList = $queues[$i] -join ", "
    Write-Host "  Worker $i ($($queues[$i].Count) subjects): $queueList" -ForegroundColor Green
}

# ============================================================================
# Section E: Work dirs and log dirs
# ============================================================================

Write-Host "`n[E] Preparing directories..." -ForegroundColor Yellow

# Create per-worker work directories
for ($i = 0; $i -lt $Parallel; $i++) {
    $workerWorkDir = Join-Path $WORK_DIR "worker_$i"
    if (-not (Test-Path $workerWorkDir)) {
        New-Item -ItemType Directory -Path $workerWorkDir -Force | Out-Null
    }
}
Write-Host "  Work dirs: $WORK_DIR\worker_0 .. worker_$($Parallel - 1)" -ForegroundColor Green

# Create output dir if needed
if (-not (Test-Path $OUTPUT_DIR)) {
    New-Item -ItemType Directory -Path $OUTPUT_DIR -Force | Out-Null
    Write-Host "  Created output dir: $OUTPUT_DIR" -ForegroundColor Green
}

# Pre-seed fsaverage to avoid race condition when multiple containers
# try to copy it simultaneously into the shared output directory
$fsAvgDir = Join-Path $OUTPUT_DIR "sourcedata\freesurfer\fsaverage\surf"
if (-not (Test-Path $fsAvgDir)) {
    Write-Host "  Pre-seeding fsaverage directory (one-time setup)..." -ForegroundColor Yellow
    $dockerOutputForSeed = $OUTPUT_DIR -replace '\\','/'
    docker run --rm `
        -v "${dockerOutputForSeed}:/out" `
        --entrypoint /bin/bash `
        $fmriprepImage `
        -c "mkdir -p /out/sourcedata/freesurfer && cp -rn /opt/freesurfer/subjects/fsaverage /out/sourcedata/freesurfer/"
    Write-Host "  fsaverage seeded" -ForegroundColor Green
} else {
    Write-Host "  fsaverage: already exists" -ForegroundColor Green
}

# Create log dir
$runLogDir = Join-Path $LOG_DIR "fmriprep_${Dataset}_${Timestamp}"
New-Item -ItemType Directory -Path $runLogDir -Force | Out-Null
Write-Host "  Log dir: $runLogDir" -ForegroundColor Green

# ============================================================================
# DRY RUN: Stop here if -DryRun
# ============================================================================

if ($DryRun) {
    Write-Host "`n========================================" -ForegroundColor Cyan
    Write-Host " DRY RUN SUMMARY" -ForegroundColor Cyan
    Write-Host "========================================" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "Would process $($subjects.Count) subjects with $Parallel parallel workers."
    Write-Host "Per worker: $cpusPerWorker CPUs, $memGBPerWorker GB RAM."
    Write-Host ""

    # Show sample Docker command for first subject
    $sampleSubject = $subjects[0]
    $sampleWorkDir = Join-Path $WORK_DIR "worker_0"

    # Convert Windows paths to Docker-compatible format
    $dockerBidsRoot = $BIDS_ROOT -replace '\\','/'
    $dockerOutputDir = $OUTPUT_DIR -replace '\\','/'
    $dockerLicenseFile = $LICENSE_FILE -replace '\\','/'
    $dockerWorkDir = $sampleWorkDir -replace '\\','/'
    $dockerTemplateFlow = $TEMPLATEFLOW_DIR -replace '\\','/'
    Write-Host "Sample Docker command (worker 0, sub-$sampleSubject):" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "docker run --rm ``" -ForegroundColor White
    Write-Host "    --name `"fmriprep_w0_$sampleSubject`" ``" -ForegroundColor White
    Write-Host "    -v `"${dockerBidsRoot}:/data:ro`" ``" -ForegroundColor White
    Write-Host "    -v `"${dockerOutputDir}:/out`" ``" -ForegroundColor White
    Write-Host "    -v `"${dockerLicenseFile}:/opt/freesurfer/license.txt:ro`" ``" -ForegroundColor White
    Write-Host "    -v `"${dockerWorkDir}:/work`" ``" -ForegroundColor White
    if ($TEMPLATEFLOW_DIR -and (Test-Path $TEMPLATEFLOW_DIR)) {
        Write-Host "    -v `"${dockerTemplateFlow}:/templateflow`" ``" -ForegroundColor White
        Write-Host "    -e `"TEMPLATEFLOW_HOME=/templateflow`" ``" -ForegroundColor White
    }
    Write-Host "    $fmriprepImage ``" -ForegroundColor White
    Write-Host "    /data /out participant ``" -ForegroundColor White
    Write-Host "    --participant-label $sampleSubject ``" -ForegroundColor White
    Write-Host "    -w /work ``" -ForegroundColor White
    Write-Host "    --fs-license-file /opt/freesurfer/license.txt ``" -ForegroundColor White
    Write-Host "    --output-spaces $OUTPUT_SPACES ``" -ForegroundColor White
    Write-Host "    --nprocs $cpusPerWorker ``" -ForegroundColor White
    Write-Host "    --mem $(($memGBPerWorker * 1024)) ``" -ForegroundColor White
    Write-Host "    --skip_bids_validation ``" -ForegroundColor White
    Write-Host "    --fs-no-reconall ``" -ForegroundColor White
    if ($EXTRA_ARGS) {
        Write-Host "    $EXTRA_ARGS ``" -ForegroundColor White
    }
    Write-Host "    --notrack" -ForegroundColor White
    Write-Host ""

    # Estimate wall time
    $estPerSubjectHours = 2.5  # conservative estimate
    $subjectsPerWorker = [math]::Ceiling($subjects.Count / $Parallel)
    $estTotalHours = [math]::Round($subjectsPerWorker * $estPerSubjectHours, 1)
    Write-Host "Estimated wall time: ~$estTotalHours hours ($subjectsPerWorker subjects/worker x ~${estPerSubjectHours}h each)" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "To run for real, remove the -DryRun flag." -ForegroundColor Cyan
    exit 0
}

# ============================================================================
# Section F: Generate worker scripts
# ============================================================================

Write-Host "`n[F] Generating worker scripts..." -ForegroundColor Yellow

$workerScripts = @()
$workerLogFiles = @()

for ($w = 0; $w -lt $Parallel; $w++) {
    $workerQueue = $queues[$w]
    if ($workerQueue.Count -eq 0) { continue }

    $workerScript = Join-Path $runLogDir "worker_${w}.ps1"
    $workerLog = Join-Path $runLogDir "worker_${w}.log"
    $workerWorkDir = Join-Path $WORK_DIR "worker_$w"

    # Convert Windows paths to Docker-compatible format
    $dockerBidsRoot = $BIDS_ROOT -replace '\\','/'
    $dockerOutputDir = $OUTPUT_DIR -replace '\\','/'
    $dockerLicenseFile = $LICENSE_FILE -replace '\\','/'
    $dockerWorkDir = $workerWorkDir -replace '\\','/'
    $dockerTemplateFlow = $TEMPLATEFLOW_DIR -replace '\\','/'
    # Build the TemplateFlow volume mount lines
    $templateFlowLine1 = ""
    $templateFlowLine2 = ""
    if ($TEMPLATEFLOW_DIR -and (Test-Path $TEMPLATEFLOW_DIR)) {
        $templateFlowLine1 = "    -v `"${dockerTemplateFlow}:/templateflow`" ``"
        $templateFlowLine2 = "    -e `"TEMPLATEFLOW_HOME=/templateflow`" ``"
    }

    # Build per-subject processing blocks
    $subjectBlocks = ""
    $subjectIdx = 0
    foreach ($subId in $workerQueue) {
        $subjectIdx++
        $containerName = "fmriprep_w${w}_${subId}"

        # Build docker run command as an array of lines
        $dockerLines = @(
            "docker run --rm ``",
            "    --name `"$containerName`" ``",
            "    -v `"${dockerBidsRoot}:/data:ro`" ``",
            "    -v `"${dockerOutputDir}:/out`" ``",
            "    -v `"${dockerLicenseFile}:/opt/freesurfer/license.txt:ro`" ``",
            "    -v `"${dockerWorkDir}:/work`" ``"
        )
        if ($templateFlowLine1) {
            $dockerLines += $templateFlowLine1
            $dockerLines += $templateFlowLine2
        }
        $dockerLines += @(
            "    $fmriprepImage ``",
            "    /data /out participant ``",
            "    --participant-label $subId ``",
            "    -w /work ``",
            "    --fs-license-file /opt/freesurfer/license.txt ``",
            "    --output-spaces $OUTPUT_SPACES ``",
            "    --nprocs $cpusPerWorker ``",
            "    --mem $(($memGBPerWorker * 1024)) ``",
            "    --skip_bids_validation ``",
            "    --fs-no-reconall ``"
        )
        # Add dataset-specific extra args (e.g., --skull-strip-t1w skip)
        if ($EXTRA_ARGS) {
            $dockerLines += "    $EXTRA_ARGS ``"
        }
        $dockerLines += @(
            "    --notrack"
        )
        $dockerCmd = $dockerLines -join "`n"

        $subjectBlocks += @"

# --- Subject $subjectIdx / $($workerQueue.Count): sub-$subId ---
Write-Host ""
Write-Host "================================================================" -ForegroundColor Cyan
Write-Host " Worker $w - Starting sub-$subId ($subjectIdx / $($workerQueue.Count)) at `$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" -ForegroundColor Cyan
Write-Host "================================================================" -ForegroundColor Cyan

# Check if already processed (in case another worker finished it)
if (Test-Path "$($OUTPUT_DIR -replace '\\','/')/sub-${subId}.html") {
    Write-Host "  sub-$subId already has HTML report - SKIPPING" -ForegroundColor Yellow
} else {

# Remove any leftover container with same name
docker rm -f "$containerName" 2>`$null | Out-Null

`$startTime = Get-Date

$dockerCmd

`$exitCode = `$LASTEXITCODE
`$elapsed = (Get-Date) - `$startTime
`$elapsedStr = "{0:hh\:mm\:ss}" -f `$elapsed

if (`$exitCode -eq 0) {
    Write-Host "  sub-$subId COMPLETED in `$elapsedStr" -ForegroundColor Green
    `$completed++
} else {
    Write-Host "  sub-$subId FAILED (exit code `$exitCode) after `$elapsedStr" -ForegroundColor Red
    `$failed++
}

}  # end else (not already processed)

Write-Host "  Worker $w progress: `$(`$completed + `$failed) / `$totalInQueue (completed: `$completed, failed: `$failed)" -ForegroundColor Cyan

"@
    }

    $workerContent = @"
# Auto-generated worker script for fMRIPrep parallel processing
# Worker: $w | Dataset: $Dataset | Generated: $Timestamp
# Subjects: $($workerQueue -join ', ')

`$ErrorActionPreference = "Continue"

Write-Host "Worker $w started at `$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" -ForegroundColor Cyan
Write-Host "Queue: $($workerQueue -join ', ')" -ForegroundColor Cyan
Write-Host "Resources: $cpusPerWorker CPUs, $memGBPerWorker GB RAM" -ForegroundColor Cyan

`$workerStart = Get-Date
`$totalInQueue = $($workerQueue.Count)
`$completed = 0
`$failed = 0
$subjectBlocks

`$workerElapsed = (Get-Date) - `$workerStart
`$workerElapsedStr = "{0:dd\.hh\:mm\:ss}" -f `$workerElapsed

Write-Host ""
Write-Host "================================================================" -ForegroundColor Cyan
Write-Host " Worker $w FINISHED at `$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" -ForegroundColor Cyan
Write-Host " Total time: `$workerElapsedStr" -ForegroundColor Cyan
Write-Host " Completed: `$completed / `$totalInQueue | Failed: `$failed" -ForegroundColor Cyan
Write-Host "================================================================" -ForegroundColor Cyan
"@

    Set-Content -Path $workerScript -Value $workerContent -Encoding UTF8
    $workerScripts += $workerScript
    $workerLogFiles += $workerLog
    Write-Host "  Worker $w script: $workerScript" -ForegroundColor Green
}

# ============================================================================
# Section G: Launch workers
# ============================================================================

Write-Host "`n[G] Launching $($workerScripts.Count) workers..." -ForegroundColor Yellow

$workerProcesses = @()
for ($w = 0; $w -lt $workerScripts.Count; $w++) {
    # Wrap worker script to tee output to log file while still showing in window
    $wrapperArgs = "-NoProfile -ExecutionPolicy Bypass -Command `"& { . '$($workerScripts[$w])' } 2>&1 | Tee-Object -FilePath '$($workerLogFiles[$w])'`""

    $proc = Start-Process -FilePath "powershell.exe" `
        -ArgumentList $wrapperArgs `
        -PassThru `
        -WindowStyle Normal

    $workerProcesses += $proc
    Write-Host "  Worker $w launched (PID: $($proc.Id), log: $($workerLogFiles[$w]))" -ForegroundColor Green

    # Small delay between launches to avoid Docker race conditions
    Start-Sleep -Seconds 2
}

$launchTime = Get-Date
Write-Host "`n  All workers launched at $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" -ForegroundColor Cyan

# ============================================================================
# Section H: Monitor progress
# ============================================================================

Write-Host "`n[H] Monitoring progress (every ${MonitorIntervalSec}s)..." -ForegroundColor Yellow
Write-Host "  Press Ctrl+C to stop monitoring (workers continue in background)`n" -ForegroundColor Gray

$totalToProcess = $subjects.Count

# Monitor loop
while ($true) {
    Start-Sleep -Seconds $MonitorIntervalSec

    # Count HTML reports
    $htmlReports = @(Get-ChildItem -Path $OUTPUT_DIR -Filter "sub-*.html" -ErrorAction SilentlyContinue)
    $currentProcessed = $htmlReports.Count

    # Check running Docker containers
    $runningContainers = @(docker ps --filter "name=fmriprep_w" --format "{{.Names}}" 2>$null)

    # Check worker processes
    $aliveWorkers = @($workerProcesses | Where-Object { -not $_.HasExited })

    $elapsed = (Get-Date) - $launchTime
    $elapsedStr = "{0:dd\.hh\:mm\:ss}" -f $elapsed

    # Calculate newly processed since launch
    $initialProcessed = $processedCount
    $newlyProcessed = $currentProcessed - $initialProcessed
    $remaining = $totalToProcess - $newlyProcessed

    Write-Host "[$(Get-Date -Format 'HH:mm:ss')] Elapsed: $elapsedStr | Processed: $newlyProcessed / $totalToProcess | Running containers: $($runningContainers.Count) | Active workers: $($aliveWorkers.Count)" -ForegroundColor Cyan

    if ($runningContainers.Count -gt 0) {
        Write-Host "  Containers: $($runningContainers -join ', ')" -ForegroundColor Gray
    }

    # ETA estimate
    if ($newlyProcessed -gt 0) {
        $avgSecsPerSubject = $elapsed.TotalSeconds / $newlyProcessed
        $etaSecs = $remaining * $avgSecsPerSubject / [math]::Max(1, $runningContainers.Count)
        $eta = [TimeSpan]::FromSeconds($etaSecs)
        $etaStr = "{0:dd\.hh\:mm\:ss}" -f $eta
        Write-Host "  ETA: ~$etaStr" -ForegroundColor Gray
    }

    # Check if all workers are done
    if ($aliveWorkers.Count -eq 0) {
        Write-Host "`n  All workers have finished." -ForegroundColor Green
        break
    }
}

# ============================================================================
# Section I: Summary
# ============================================================================

$totalElapsed = (Get-Date) - $launchTime
$totalElapsedStr = "{0:dd\.hh\:mm\:ss}" -f $totalElapsed

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host " PROCESSING COMPLETE" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

# Final count of HTML reports
$htmlReports = @(Get-ChildItem -Path $OUTPUT_DIR -Filter "sub-*.html" -ErrorAction SilentlyContinue)
$finalProcessed = $htmlReports.Count
$thisRunProcessed = $finalProcessed - $processedCount

Write-Host ""
Write-Host "  Dataset: $($config.dataset_name) ($Dataset)" -ForegroundColor White
Write-Host "  Total time: $totalElapsedStr" -ForegroundColor White
Write-Host "  Workers used: $Parallel" -ForegroundColor White
Write-Host "  Subjects attempted: $totalToProcess" -ForegroundColor White
Write-Host "  Subjects completed this run: $thisRunProcessed" -ForegroundColor White
Write-Host "  Total subjects with reports: $finalProcessed / $totalSubjects" -ForegroundColor White

# Identify failed subjects
$processedSubjectIds = $htmlReports | ForEach-Object {
    if ($_.Name -match '^sub-(.+)\.html$') { $Matches[1] }
}

$failedSubjects = @($subjects | Where-Object { $_ -notin $processedSubjectIds })

if ($failedSubjects.Count -gt 0) {
    Write-Host ""
    Write-Host "  FAILED subjects ($($failedSubjects.Count)):" -ForegroundColor Red
    foreach ($fs in $failedSubjects) {
        Write-Host "    - sub-$fs" -ForegroundColor Red
    }
    Write-Host ""
    Write-Host "  Check logs at: $runLogDir" -ForegroundColor Yellow
    Write-Host "  To retry failed subjects, run this script again (auto-resume)." -ForegroundColor Yellow
} else {
    Write-Host ""
    Write-Host "  All subjects processed successfully!" -ForegroundColor Green
}

# Check worker exit codes
Write-Host ""
Write-Host "  Worker exit codes:" -ForegroundColor White
for ($w = 0; $w -lt $workerProcesses.Count; $w++) {
    $exitCode = $workerProcesses[$w].ExitCode
    $status = if ($exitCode -eq 0) { "OK" } else { "Exit code: $exitCode" }
    $color = if ($exitCode -eq 0) { "Green" } else { "Red" }
    Write-Host "    Worker $w : $status" -ForegroundColor $color
}

Write-Host ""
Write-Host "  Logs: $runLogDir" -ForegroundColor White
Write-Host ""
