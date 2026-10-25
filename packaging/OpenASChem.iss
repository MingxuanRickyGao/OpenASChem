[Setup]
AppId={code:ApplicationId}
UsePreviousLanguage=no
AppName=OpenASChem
AppVersion=1.3.0
AppPublisher=Mingxuan Gao
DefaultDirName={localappdata}\Programs\OpenASChem-1.3.0
DefaultGroupName=OpenASChem 1.3.0
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0
WizardStyle=modern
SetupIconFile=..\src\openaschem\assets\openaschem.ico
UninstallDisplayIcon={app}\openaschem.ico
Compression=lzma2/fast
SolidCompression=yes
OutputDir=..\..\release
OutputBaseFilename=OpenASChem-1.3.0-Windows-x64-Setup
LicenseFile=..\LICENSE
InfoBeforeFile=INSTALLATION.txt

[Files]
Source: "{#RuntimeDir}\*"; DestDir: "{app}\environment"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "..\LICENSE"; DestDir: "{app}"
Source: "..\NOTICE"; DestDir: "{app}"
Source: "..\src\openaschem\assets\openaschem.ico"; DestDir: "{app}"

[Tasks]
Name: desktopicon; Description: "Create a desktop shortcut"

[Icons]
Name: "{group}\OpenASChem"; Filename: "{app}\environment\pythonw.exe"; Parameters: "-E -c ""from openaschem.gui import main; main()"""; WorkingDir: "{userdocs}"; IconFilename: "{app}\openaschem.ico"
Name: "{group}\OpenASChem file filling"; Filename: "{app}\environment\pythonw.exe"; Parameters: "-E -c ""from openaschem.gui import filling_main; filling_main()"""; WorkingDir: "{userdocs}"; IconFilename: "{app}\openaschem.ico"
Name: "{group}\Uninstall OpenASChem"; Filename: "{uninstallexe}"
Name: "{userdesktop}\OpenASChem 1.3.0"; Filename: "{app}\environment\pythonw.exe"; Parameters: "-E -c ""from openaschem.gui import main; main()"""; WorkingDir: "{userdocs}"; IconFilename: "{app}\openaschem.ico"; Tasks: desktopicon

[UninstallDelete]
Type: filesandordirs; Name: "{app}\environment"
Type: filesandordirs; Name: "{app}\bin"
Type: files; Name: "{app}\openaschem-install.json"
Type: files; Name: "{app}\installation.log"

[Code]
var InstallFailed: Boolean;
function ApplicationId(Param: String): String;
begin
  if ExpandConstant('{param:QAMODE|0}') = '1' then
    Result := '{A481E0E2-8148-4325-9BCB-7D2B7D0DB207}'
  else Result := '{C1E5521A-14B3-4E12-AAD5-3DB875E05BB8}';
end;
function PrepareToInstall(var NeedsRestart: Boolean): String;
var Marker: AnsiString;
begin
  Result := '';
  if DirExists(ExpandConstant('{app}\environment')) then begin
    if not LoadStringFromFile(ExpandConstant('{app}\openaschem-install.json'), Marker) then
      Result := 'Choose a new empty installation directory.'
    else if (Pos('bundled', String(Marker)) = 0) and (Pos('fully-offline', String(Marker)) = 0) then
      Result := 'Choose a new directory for this installation.';
  end;
end;
function GetCustomSetupExitCode: Integer;
begin
  if InstallFailed then Result := 1 else Result := 0;
end;
procedure CurStepChanged(CurStep: TSetupStep);
var ExitCode: Integer; Args: String;
begin
  if CurStep = ssPostInstall then begin
    WizardForm.StatusLabel.Caption := 'Preparing bundled Python and checking OpenASChem. No downloads are required.';
    Args := '-E "' + ExpandConstant('{app}\environment\openaschem-finalize.py') + '" --prefix "' + ExpandConstant('{app}') + '"';
    if not Exec(ExpandConstant('{app}\environment\python.exe'), Args, '', SW_HIDE, ewWaitUntilTerminated, ExitCode) then begin
      InstallFailed := True;
      RaiseException('Could not start bundled Python.');
    end;
    if ExitCode <> 0 then begin
      InstallFailed := True;
      RaiseException('Local environment check failed. See installation.log.');
    end;
  end;
end;
