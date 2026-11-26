$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath('D:\projects\TMS_GWT\manuscript\aachen_brief_report\corrected_20261126')
$taskWord = New-Object -ComObject Word.Application
$taskWord.Visible = $false
$taskWord.DisplayAlerts = 0
$taskWord.AutomationSecurity = 3
try {
    foreach ($taskItem in @(@('Aachen_Brief_Research_Report_Corrected.docx','main'), @('Supplementary_Material_Corrected.docx','supplement'))) {
        $taskInput = Join-Path $taskRoot $taskItem[0]
        $taskFolder = Join-Path $taskRoot ('qa\' + $taskItem[1])
        [void](New-Item -ItemType Directory -Path $taskFolder -Force)
        $taskDocument = $taskWord.Documents.Open($taskInput,$false,$true,$false)
        try {
            $taskDocument.Repaginate()
            $taskDocument.ExportAsFixedFormat((Join-Path $taskFolder 'document.pdf'),17)
            Write-Output ('Rendered ' + $taskItem[1])
        } finally {
            $taskDocument.Close(0)
            [void][Runtime.InteropServices.Marshal]::ReleaseComObject($taskDocument)
        }
    }
} finally {
    $taskWord.Quit(0)
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($taskWord)
}
