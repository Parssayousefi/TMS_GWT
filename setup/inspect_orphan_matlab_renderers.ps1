[CmdletBinding()]
param(
    [switch]$Apply,
    [datetime]$CreatedAfter = [datetime]'2026-11-25T09:20:00Z'
)

# Only obsolete CEF renderer processes are eligible. Never stop a MATLAB
# numerical process, a browser parent, or a renderer with a surviving parent.
$ErrorActionPreference = 'Stop'
$taskRoot = Split-Path -Parent $PSScriptRoot
$taskExpectedExe = 'C:\Program Files\MATLAB\R2025a\bin\win64\MATLABWindow.exe'
$taskProcesses = @(Get-CimInstance Win32_Process)
$taskProcessIds = @($taskProcesses.ProcessId)
$taskCandidates = @($taskProcesses | Where-Object {
    $_.Name -eq 'MATLABWindow.exe' -and
    $_.ExecutablePath -eq $taskExpectedExe -and
    $_.CommandLine -match '--type=renderer' -and
    $_.CreationDate -ge $CreatedAfter -and
    $_.ParentProcessId -notin $taskProcessIds
})
$taskRecords = @()
foreach ($taskCandidate in $taskCandidates) {
    $taskCurrent = Get-CimInstance Win32_Process -Filter "ProcessId=$($taskCandidate.ProcessId)"
    if (-not $taskCurrent) { continue }
    # Recheck identity to reject PID reuse and refresh the missing-parent check.
    if ($taskCurrent.CreationDate -ne $taskCandidate.CreationDate -or
        $taskCurrent.ExecutablePath -ne $taskExpectedExe -or
        $taskCurrent.CommandLine -notmatch '--type=renderer' -or
        $taskCurrent.ParentProcessId -ne $taskCandidate.ParentProcessId) { continue }
    $taskParent = Get-CimInstance Win32_Process -Filter "ProcessId=$($taskCurrent.ParentProcessId)"
    if ($taskParent) { continue }
    if ($Apply) { Stop-Process -Id $taskCurrent.ProcessId }
    $taskRecords += [pscustomobject]@{
        pid = $taskCurrent.ProcessId
        created = $taskCurrent.CreationDate
        missing_parent_pid = $taskCurrent.ParentProcessId
        action = $(if ($Apply) { 'stopped_orphan_renderer' } else { 'eligible_orphan_renderer' })
        checked_utc = [datetime]::UtcNow.ToString('o')
    }
}
if ($Apply -and $taskRecords.Count) {
    $taskStamp = [datetime]::UtcNow.ToString('yyyyMMddTHHmmss_fffffffZ')
    $taskLog = Join-Path $taskRoot "logs/tms_pipeline_20261125/orphan_renderers_$taskStamp.json"
    ConvertTo-Json -InputObject $taskRecords -Depth 3 | Set-Content -LiteralPath $taskLog
}
ConvertTo-Json -InputObject $taskRecords -Compress
