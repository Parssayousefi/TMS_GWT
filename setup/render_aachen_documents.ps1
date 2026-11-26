param([ValidateSet('all', 'main', 'supplement')][string]$Kind = 'all', [string]$Suffix = '')
$ErrorActionPreference = 'Stop'
if ($Suffix -notmatch '^[A-Za-z0-9_-]*$') { throw 'Invalid suffix' }
$taskRoot = [IO.Path]::GetFullPath('D:\projects\TMS_GWT')
$taskFolder = Join-Path $taskRoot 'manuscript\aachen_brief_report'
$taskWord = New-Object -ComObject Word.Application
$taskWord.Visible = $false
$taskWord.DisplayAlerts = 0
$taskWord.AutomationSecurity = 3
try {
    foreach ($taskItem in @(@(('Aachen_Brief_Research_Report' + $Suffix + '.docx'), 'main'), @(('Supplementary_Material' + $Suffix + '.docx'), 'supplement'))) {
        if ($Kind -ne 'all' -and $taskItem[1] -ne $Kind) { continue }
        $taskInput = [IO.Path]::GetFullPath((Join-Path $taskFolder $taskItem[0]))
        $taskOutputFolder = [IO.Path]::GetFullPath((Join-Path $taskFolder ('qa\' + $taskItem[1])))
        if (-not $taskInput.StartsWith($taskRoot + '\', [StringComparison]::OrdinalIgnoreCase) -or -not $taskOutputFolder.StartsWith($taskRoot + '\', [StringComparison]::OrdinalIgnoreCase)) { throw 'Output is outside the project' }
        [void](New-Item -ItemType Directory -Path $taskOutputFolder -Force)
        $taskPdf = Join-Path $taskOutputFolder 'document.pdf'
        $taskDocument = $taskWord.Documents.Open($taskInput, $false, $true, $false)
        try {
            $taskDocument.Repaginate()
            $taskDocument.ExportAsFixedFormat($taskPdf, 17)
            Write-Output ('Exported ' + $taskPdf)
        } finally {
            $taskDocument.Close(0)
            [void][Runtime.InteropServices.Marshal]::ReleaseComObject($taskDocument)
        }
    }
} finally {
    $taskWord.Quit(0)
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($taskWord)
}
