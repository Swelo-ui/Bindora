; Inno Setup Script for Bindora Dock Standalone Desktop Installer
; Version: 2.0
; Prepared for: Himanshu Sharma (NexPharmaTech)

#define MyAppName "Bindora Dock"
#define MyAppVersion "2.0"
#define MyAppPublisher "NexPharmaTech"
#define MyAppURL "https://github.com/Swelo-ui/Bindora"
#define MyAppExeName "bindora_launcher.exe"

[Setup]
AppId={{9F7B16E4-749F-4CB3-9564-96E81A6F824A}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}/releases
DefaultDirName={autopf}\Bindora
DefaultGroupName=Bindora
DisableProgramGroupPage=yes
LicenseFile=LICENSE
OutputDir=installer\Output
OutputBaseFilename=BindoraDock-Setup
SetupIconFile=frontend\assets\branding\favicon.ico
WizardImageFile=installer\wizard-image.bmp
WizardSmallImageFile=installer\wizard-small.bmp
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
ArchitecturesInstallIn64BitMode=x64compatible
PrivilegesRequired=admin
UninstallDisplayIcon={app}\{#MyAppExeName}

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[Files]
Source: "dist\bindora_launcher\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"
Name: "{group}\{cm:UninstallProgram,{#MyAppName}}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "{cm:LaunchProgram,{#StringChange(MyAppName, '&', '&&')}}"; Flags: nowait postinstall skipifsilent
