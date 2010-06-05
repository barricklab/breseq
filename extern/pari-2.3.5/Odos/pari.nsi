!include "MUI.nsh"
Name "PARI 2.3.3"
!define dll "libpari-gmp.dll"
!define objdir "..\Ocygwin-i686-gmp"

;--No need to modify things below --
AutoCloseWindow false

OutFile "Pari.exe"
InstallDir "$PROGRAMFILES\PARI"
InstallDirRegKey HKLM "Software\PARI" ""

!define MUI_ABORTWARNING

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "..\COPYING"
!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_WELCOME
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_UNPAGE_FINISH

!insertmacro MUI_LANGUAGE "English"
;--------------------------------
;Installer Sections

!define uninst "Software\Microsoft\Windows\CurrentVersion\Uninstall\PARI"

Section "pari (required)" SecCopy
  SetOutPath "$INSTDIR"
  File /oname=gp.exe "${objdir}\gp-dyn.exe"
  File "makegprc"
  File "..\misc\tex2mail"
  File "${objdir}\${dll}"
  File "\cygwin\bin\cygcrypt-0.dll"
  File "\cygwin\bin\cygfltknox-1.1.dll"
  File "\cygwin\bin\cygiconv-2.dll"
  File "\cygwin\bin\cygintl-8.dll"
  File "\cygwin\bin\cygncurses-8.dll"
  File "\cygwin\bin\cygreadline6.dll"
  File "\cygwin\bin\cygperl5_8.dll"
  File "\cygwin\bin\cygwin1.dll"
  File "\cygwin\bin\perl.exe"
  File "\cygwin\bin\sh.exe"

  WriteRegStr HKCU "Software\PARI" "" $INSTDIR
  WriteRegStr HKLM ${uninst} "DisplayName" "PARI (remove only)"
  WriteRegStr HKLM ${uninst} "UninstallString" '"$INSTDIR\uninstall.exe"'
  
  WriteUninstaller "$INSTDIR\Uninstall.exe"

  SetOutPath "$INSTDIR"
  ExecWait "perl makegprc $INSTDIR"
  Delete "$INSTDIR\makegprc"
SectionEnd

Section "Galois files" SecGAL
  SetOutPath "$INSTDIR\galdata"
  File "..\data\galdata\*"
SectionEnd

Section "documentation" SecDOC
  SetOutPath "$INSTDIR"
  File "..\doc\gphelp"
  SetOutPath $INSTDIR\doc
  File "acro.exe"
  File "..\doc\translations"
  File "..\doc\*.tex"
  File "..\doc\*.pdf"
SectionEnd

Section "examples" SecEX
  SetOutPath "$INSTDIR"
  File "..\doc\gphelp"
  SetOutPath $INSTDIR\examples
  File "..\examples\EXPLAIN"
  File "..\examples\Inputrc"
  File "..\examples\*.gp"
  File "..\examples\*.c"
  File "..\examples\Makefile.cygwin-i686"
SectionEnd

Function .onInstSuccess
  MessageBox MB_OK "Thank you for using PARI/GP! Double-click on 'gp' to start the calculator.$\r$\nTweak $INSTDIR\.gprc to customize GP: colors, script search path, etc."
  ExecShell "open" "$INSTDIR"
FunctionEnd

!define short "$SMPROGRAMS\PARI"
  
Section "shortcuts" SecSM
  CreateDirectory "${short}"
  CreateShortCut "${short}\gp.lnk" "$INSTDIR\gp.exe" "" "$INSTDIR\gp.exe" 0
  CreateShortCut "${short}\users.lnk" "$INSTDIR\doc\users.pdf" "" "$INSTDIR\doc\users.pdf" 0
  CreateShortCut "${short}\libpari.lnk" "$INSTDIR\doc\libpari.pdf" "" "$INSTDIR\doc\libpari.pdf" 0
  CreateShortCut "${short}\tutorial.lnk" "$INSTDIR\doc\tutorial.pdf" "" "$INSTDIR\doc\tutorial.pdf" 0
  CreateShortCut "${short}\refcard.lnk" "$INSTDIR\doc\refcard.pdf" "" "$INSTDIR\doc\refcard.pdf" 0
  WriteINIStr "${short}\PARI pages.url" "InternetShortcut" "URL" "http://pari.math.u-bordeaux.fr"
  CreateShortCut "${short}\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0
  CreateShortCut "$DESKTOP\PARI.lnk" "$INSTDIR\gp.exe"
SectionEnd

;--------------------------------
;Descriptions

LangString DESC_SecCopy ${LANG_ENGLISH} "Copy pari files to application folder."
LangString DESC_DOC ${LANG_ENGLISH} "Install documentation and online help."
LangString DESC_EX ${LANG_ENGLISH} "Install sample GP scripts."
LangString DESC_GAL ${LANG_ENGLISH} "Install Galois data files (degree > 7)."
LangString DESC_SM ${LANG_ENGLISH} "Add PARI shortcuts to Start Menu and desktop."

!insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
  !insertmacro MUI_DESCRIPTION_TEXT ${SecCopy} $(DESC_SecCopy)
  !insertmacro MUI_DESCRIPTION_TEXT ${SecGAL} $(DESC_GAL)
  !insertmacro MUI_DESCRIPTION_TEXT ${SecSM} $(DESC_SM)
  !insertmacro MUI_DESCRIPTION_TEXT ${SecDOC} $(DESC_DOC)
  !insertmacro MUI_DESCRIPTION_TEXT ${SecEX} $(DESC_EX)
!insertmacro MUI_FUNCTION_DESCRIPTION_END
 
;--------------------------------
Section "Uninstall"
  Delete "$INSTDIR\gp.exe"
  Delete "$INSTDIR\.gprc"
  Delete "$INSTDIR\gphelp"
  Delete "$INSTDIR\tex2mail"
  Delete "$INSTDIR\${dll}"
  Delete "$INSTDIR\cygcrypt-0.dll"
  Delete "$INSTDIR\cygfltknox-1.1.dll"
  Delete "$INSTDIR\cygiconv-2.dll"
  Delete "$INSTDIR\cygintl-8.dll"
  Delete "$INSTDIR\cygncurses-8.dll"
  Delete "$INSTDIR\cygreadline6.dll"
  Delete "$INSTDIR\cygperl5_8.dll"
  Delete "$INSTDIR\cygwin1.dll"
  Delete "$INSTDIR\perl.exe"
  Delete "$INSTDIR\sh.exe"

  Delete "$INSTDIR\Uninstall.exe"
  RMDir /r "$INSTDIR\doc"
  RMDir /r "$INSTDIR\examples"
  RMDir /r "$INSTDIR\galdata"

  DeleteRegKey HKLM ${uninst}
  DeleteRegKey /ifempty HKLM "Software\PARI"

  RMDir /r "$SMPROGRAMS\PARI"
  Delete "$DESKTOP\PARI.lnk"
  RMDir "$INSTDIR"
SectionEnd
