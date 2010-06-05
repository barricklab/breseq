;; $Id: pari-messages.el 4286 2003-03-01 16:18:39Z karim $
;; pari-messages.el --  part of pari.el GP/PARI editing support package.

;; documentation functions
;; See pariemacs.txt  for more details.
 
(provide 'pari-messages)
;; The only function of this script which is used in pari.el is gp-messager.

(defcustom gp-language 'francais
"Best to re-tart your emacs session if you change this value."
:type '(choice (const francais) (const english) (const deutsch))
:group 'gp)

(defconst gp-messages-list
  '((francais  .
      ("Nous utilisons le choix par defaut pour la completion"  ;; no 1
       "Elimination de %s"                            ;; no 2
       " Sauvegarde des couleurs ? "                  ;; no 3
       "M-o ou ESC-o pour oter la fenetre d'aide"     ;; no 4
       "APPUYER SUR UNE TOUCHE POUR CONTINUER..."     ;; no 5
       "termine."                                     ;; no 6
       "En attente de la reponse de gp ..."           ;; no 7
       "Impossible de lancer gp."                     ;; no 8
       "Expression incomplete."                       ;; no 9
       "Ce nouveau prompt peut conduire a une erreur. Mieux vaut le changer via M-\\p"  ;; no 10
        "Version numero %s."                          ;; no 11
        ""                         ;; no 12
        ""                      ;; no 13
        "Echange les fonctions des touches C-p/M-p and C-n/M-n."
        "gp essaie de completer ..."                  ;; no 15
        "C-u M-o pour sortir de l'edition."           ;; no 16
        "Lancement de "                               ;; no 17
        "SPC=suivant DEL=precedent RET=selectionner s=survey-menu q=quitter"
        "Il n'y a rien a selectionner ici"            ;; no 19
        "Fonction inconnue : %s"                      ;; no 20
        "Fonction"                                    ;; no 21
        "Variable utilisateur : %s"                   ;; no 22
        "Aucune occurence de \"%s\" n'a ete trouvee." ;; no 23
        "Chargement de pari-colors.el a partir de "   ;; no 24
        "### Variables globales : (une par ligne)"    ;; no 25
        "### Titres de chapitre :"                    ;; no 26
        "### Mots-cles interessants :"                ;; no 27
        "D'humeur aventureuse ? Determiner tout d'abord le type d'objets que vous souhaitez `colorer'. Si ce type est signale comme ayant une valeur par defaut, vous devez alors choisir soit de changer la valeur par defaut (ce qui la changera partout ; c'est la strategie conseillee), soit de changer son equivalent local (par exemple gp-string localement et font-lock-string-face globalement). Si vous optez pour une modification globale, il faudra alors redemarrer emacs. Les modifications locales prennent le pas sur les modifications globales. Pour eliminer une specification locale, il vous faut eliminer la ligne correspondante de votre fichier `.emacs'.\n\nChoisissez ensuite une couleur parmi celles ci-dessous (vous pouvez aussi vous referer au fichier rgb.txt si il existe...).De plus, vous pouvez ajouter l'option gras, italique et/ou souligne : cliquez sur la case concernee et sur Toggle de facon a ce que le symbole devienne t (au lieu de nil).\n Liste des couleurs :\n" ;; no 28
        "Emacs utilise un fichier general contenant tous les noms des fonctions de PARI. En surplus, gp utilises un fichier appele nom-de-fichier.cpl des que nom-de-fichier est edite. Pour creer ce fichier, vous pouvez utilise l'item [GP Completion-File Edit-File...] de la barre de menu qui creera un fichier au format adequat contenant les noms des fonctions et des variables globales de votre programme. L'edition de ce fichier se fait via l'item [GP Completion-File Edit File...] de la barre de menu et vous pouvez aussi decider d'utiliser un autre fichier a l'aide de [GP Completion-File Use Also File...]."  ;; no 29
        "Rend les noms de fonctions et ceux des variables globales du programme %s (tel qu'il existe en ce moment) utilisables pour la completion. Ils seront stockes dans `%s.cpl' des que ce fichier sera edite. Le fichier `%s.cpl' est au format d'un fichier de completion (i.e. format du fichier gp-menu) et est automatiquement utilise pour la completion lorsque %s est edite."  ;; no 30
        "Fonctions/Sections dans la description desquelles \"%s\" apparait :"  ;;  no 31
        "Sujet"
        "Nom du fichier de completion : "                ;; no 33
        "Erreur introuvable."                                 ;;  no 34
        "Probable typo."                                 ;;  no 35
        (concat "Aucune erreur a localiser ou buffer (" gp-reads-this-buffer ") absent") ;; no 36
        "Le guide est en construction..." ;; no 37
        "Mouse-2 ou Return pour selectionner un item." ;; no 38
        "Le fichier pariemacs.txt n'est pas dans votre load-path. Vous devriez decouvrir ou il se situe, disons dans le directory /usr/local/lib/pari/emacs/ et ajouter la ligne\n (setq load-path (concat load-path \"/usr/local/lib/pari/emacs/\"))\ndans votre fichier .emacs (creez-le si besoin est)." ;; no 39
        "Nouvelle couleur :" ;; no 40
        "Appele avec : %s\n\n"       ;; no 41 
        "" ;; no 42
        "Couleurs" ;; no 43
        "Mise a jour"  ;; no 44
        "Automatique" ;; no 45
        "Commutateur" ;; no 46
        "Tout recolorier" ;; no 47
        "MetaTouches" ;; no 48 "Metakeys"
        "Lire du Fichier..." ;; no 49 "Read from File..."
        "Ecrire sur le Fichier..." ;; no 50 "Write to File..."
        "Imprimer en..." ;; no 51 "Print in..."
        "Nouveau Prompt" ;; no 52 "New Prompt"
        "Simple"         ;; no 53 "Simple"
        "Avec l'heure"   ;; no 54 "With Time"
        "Avec la date"             ;; no 55
        "Avec un separateur"        ;; no 56
        "Lancer GP"      ;; no 57
        ""              ;; no 58
        ""  ;; no 59
        "GP evalue le fichier..."     ;; no 60
        "GP evalue la region"      ;; no 61
        "Quitter GP"       ;; no 62
        "Aide sur ce mode"     ;; no 63
        "Fichier de completion"       ;; no 64
        "Utiliser aussi le fichier..."      ;; no 65
        "Editer le fichier..."          ;; no 66
        "Creer/Mettre a jour"           ;; no 67
        "Ajouter un commentaire"           ;; no 68
        "Nom du programme:"                ;; no 69
        ""            ;; no 70
        "Disposition precedente"           ;; no 71
        "Completer"                   ;; no 72
        "Aller a l'erreur"              ;; no 73
        "Entree/Sortie..."                  ;; no 74
        "Recopier la derniere entree"            ;; no 75
        "Oter le dernier resultat"         ;; no 76
        "Oter la derniere action"         ;; no 77
        "Basculer les touches"              ;; no 78
        "Modifier..."                       ;; no 79
        "Je passe en mode confiance. Restorez le mode usuel\n en entrant /* Trust=Off */." ;; no 80
        "Tout !" ;; no 81
        "Lock mode" ;; no 82
        "Trust mode" ;; no 83
        "Lock semi-mode %s" ;; no 84
        "Trust semi-mode %s" ;; no 85
        ))
     (english .
       ("Using default choice for completion"            ;; no 1
        "Removing %s"
        " Save Colors ? "
        "M-o or ESC-o will remove the help window"
        "PRESS ANY KEY TO CONTINUE..."                   ;; no 5
        "done."
        "Waiting for gp output ..."
        "Could not start gp."
        "Incomplete expression : Not sent to gp."        ;; no 9
        "New prompt may lead to an error. Better to set it interactively via M-\\p"
        "Version number %s."
        ""                     ;; no 12
        ""                  ;; no 13
        "Exchange the bindings of the keys C-p/M-p and C-n/M-n."
        "Waiting for gp completion ..."               ;; no 15
        "C-u M-o to exit editing."                    ;; no 16
        "Starting "                                   ;; no 17
        "SPC=next DEL=previous RET=select s=survey-menu q=quit"
        "Nothing to be selected here"                 ;; no 19
        "Unknown function: %s"                        ;; no 20
        "Function"                                    ;; no 21
        "User Variable: %s"                           ;; no 22
        "No occurence of \"%s\" found."               ;; no 23
        "Loading pari-colors.el from "                ;; no 24
        "### Global Variables : (one per line)"       ;; no 25
        "### Chapter Headings:"                       ;; no 26
        "### Interesting Keywords:"                   ;; no 27
        "Feeling adventurous ? First select the kind of objects you want to `color'. If this type is said to have a default value, you then have to decide either to change this default value (a change that will concern all modes; the recommended strategy) either to change only the local value (for instance gp-string locally and font-lock-string-face globally). If you choose a global modification then you should leave emacs and re-enter it for the new value to take effect. Local modifications take precedence over global ones, and thus to go back to global (or default) determination once you've changed it locally, you should erase the corresponding line in your `.emacs' file.\n\nChoose a color among the list below (you can also have a look at the file rgb.txt if it exists). You can require your characters to be in bold face, underlined or in italic. Simply select the proper square and use Toggle to make the symbol t appear (instead of nil).\n Colors List:\n" ;; no 28
        "A general completion file  containing the name of all the PARI functions is always used. In addition to this file, gp uses a file named your-file-name.cpl when you edit your-file-name. To create this file, you can use the menu-bar item [GP Completion-File Edit-File...] which will create the proper completion-file and introduce the names of the functions and of the global variables of your program. You edit the file by using the item \"Edit File...\" and you can decide to use another completion-file as well through the item \"Use Also File...\"."  ;; no 29
        " Makes the names of functions and global variables of %s available for completion. They will be stored in `%s.cpl' as soon as this file is required for editing. The file `%s.cpl' has the format of a completion file (i.e. a gp-menu file) and is automatically used as a completion file when %s is edited."  ;; no 30
        "Functions or Sections in whose description \"%s\" appears:"  ;;  no 31
        "Subject"                                                     ;;  no 32
        "Name of the completion file: "                               ;;  no 33
        "Could not locate the error."                                 ;;  no 34
        "Probable mistake."                                           ;;  no 35
        (concat "No error to be found or missing buffer (" gp-reads-this-buffer ")") ;; no 36
        "The browser is being build..." ;; no 37
        "Mouse-2 or Return to select an item." ;; no 38
        "The file pariemacs.txt is not in your load-path. You should discover where it is, say in the directory /usr/local/lib/pari/emacs/ and add the line\n (setq load-path (concat load-path \"/usr/local/lib/pari/emacs/\"))\nto your .emacs file (create it if it doesn't already exist)." ;; no 39
        "Use face:"     ;; no 40
        "Called with: %s\n\n"       ;; no 41 
        ""      ;; no 42
        "Colors"        ;; no 43                  http://www.local.attac.org/petition/
        "Update"        ;; no 44 
        "Automatic"          ;; no 45 
        "Hilit Switch"  ;; no 46 
        "Refontify All"   ;; no 47 
        "Metakeys"  ;; no 48
        "Read from File..." ;; no 49 
        "Write to File..." ;; no 50 
        "Print in..." ;; no 51 
        "New Prompt" ;; no 52 
        "Simple"         ;; no 53 
        "With Time"   ;; no 54 
        "With Date"             ;; no 55
        "With Separator"        ;; no 56
        "Start GP session"      ;; no 57
        ""              ;; no 58
        ""  ;; no 59
        "Run GP on file..."     ;; no 60
        "Run GP in region"      ;; no 61
        "Quit GP session"       ;; no 62
        "Help on this mode"     ;; no 63
        "Completion File"       ;; no 64
        "Use Also File..."      ;; no 65
        "Edit File..."          ;; no 66
        "Make/Update"           ;; no 67
        "Add Comment"           ;; no 68
        "Name of the GP programm: "  ;; no 69
        ""            ;; no 70
        "Previous Setting"           ;; no 71
        "Complete"                   ;; no 72
        "Skip to error"              ;; no 73
        "In/Out..."                  ;; no 74
        "Copy Last Input"            ;; no 75
        "Remove Last Output"         ;; no 76
        "Remove Last Action"         ;; no 77
        "Exchange keys"              ;; no 78
        "Customize..."                  ;; no 79
        "Entering trust mode. Recover usual behaviour\n by entering /* Trust=Off */." ;; no 80
        "All" ;; no 81
        "Lock mode" ;; no 82
        "Trust mode" ;; no 83
        "Lock semi-mode %s" ;; no 84
        "Trust semi-mode %s" ;; no 85
        ))
     (deutsch  .
       ("Using default choice for completion"         ;; no 1
        "Removing %s"
        " Save Colors ? "
        "M-o or ESC-o will remove the help window"
        "PRESS ANY KEY TO CONTINUE..."                ;; no 5
        "done."
        "Waiting for gp output ..."
        "Could not start gp."
        "Incomplete expression : Not sent to gp."     ;; no 9
        "New prompt may lead to an error. Better to set it interactively via M-\\p"
        "Version number %s."                          ;; no 11
        ""                     ;; no 12
        ""                  ;; no 13
        "Exchange the bindings of the keys C-p/M-p and C-n/M-n."
        "Waiting for gp completion ..."               ;; no 15
        "C-u M-o to exit editing."                    ;; no 16
        "Starting "                                   ;; no 17
        "SPC=next DEL=previous RET=select s=survey-menu q=quit"
        "Nothing to be selected here"                 ;; no 19
        "Unknown function: %s"                        ;; no 20
        "Function"                                    ;; no 21
        "User Variable: %s"                           ;; no 22
        "No occurence of \"%s\" found."               ;; no 23
        "Loading pari-colors.el from "                ;; no 24
        "### Global Variables : (one per line)"       ;; no 25
        "### Chapter Headings:"                       ;; no 26
        "### Interesting Keywords:"                   ;; no 27
        "Feeling adventurous ? First select the kind of objects you want to `color'. If this type is said to have a default value, you then have to decide either to change this default value (a change that will concern all modes; the recommended strategy) either to change only the local value (for instance gp-string locally and font-lock-string-face globally). If you choose a global modification then you should leave emacs and re-enter it for the new value to take effect. Local modifications take precedence over global ones, and thus to go back to global (or default) determination once you've changed it locally, you should erase the corresponding line in your `.emacs' file.\n\nChoose a color among the list below (you can also have a look at the file rgb.txt if it exists). You can require your characters to be in bold face, underlined or in italic. Simply select the proper square and use Toggle to make the symbol t appear (instead of nil).\n Colors List:\n" ;; no 28
        "A general completion file  containing the name of all the PARI functions is always used. In addition to this file, gp uses a file named  your-file-name.cpl  when  you edit your-file-name. To create this file, you can use the menu-bar item [GP Completion-File Edit-File...] which will create the proper completion-file and introduce the names of the functions and of the global variables of your program. You edit the file by using the item \"Edit File...\" and you can decide to use another completion-file as well through the item \"Use Also File...\"."
        " Makes the names of functions and global variables of %s available for completion. They will be stored in `%s.cpl' as soon as this file is required for editing. The file `%s.cpl' has the format of a completion file (i.e. a gp-menu file) and is automatically used as a completion file when %s is edited."
        "Functions or Sections in whose description \"%s\" appears:"  ;;  no 31
        "Subject"
        "Name of the completion file: "
        "Could not locate the error."                                 ;;  no 34
        "Probable mistake."                                 ;;  no 35
        (concat "No error to be found or missing buffer (" gp-reads-this-buffer ")") ;; no 36
        "The browser is being build..." ;; no 37
        "Mouse-2 or Return to select an item." ;; no 38
        "The file pariemacs.txt is not in your load-path. You should discover where it is, say in the directory /usr/local/lib/pari/emacs/ and add the line\n (setq load-path (concat load-path \"/usr/local/lib/pari/emacs/\"))\nto your .emacs file (create it if it doesn't already exist)." ;; no 39
        "Use face:"     ;; no 40
        "Called with: %s\n\n"       ;; no 41 
        ""              ;; no 42
        "Colors"                ;; no 43
        "Update"                ;; no 44 
        "Automatic"                  ;; no 45 
        "Hilit Switch"          ;; no 46 
        "Refontify All"           ;; no 47 
        "Metakeys"              ;; no 48
        "Read from File..."     ;; no 49 
        "Write to File..."      ;; no 50 
        "Print in..."           ;; no 51 
        "New Prompt"            ;; no 52 
        "Simple"                ;; no 53 
        "With Time"             ;; no 54 
        "With Date"             ;; no 55
        "With Separator"        ;; no 56
        "Start GP session"      ;; no 57
        ""              ;; no 58
        ""  ;; no 59
        "Run GP on file..."     ;; no 60
        "Run GP in region"      ;; no 61
        "Quit GP session"       ;; no 62
        "Help on this mode"     ;; no 63
        "Completion File"       ;; no 64
        "Use Also File..."      ;; no 65
        "Edit File..."          ;; no 66
        "Make/Update"           ;; no 67
        "Add Comment"           ;; no 68
        "Name of the GP programm: "  ;; no 69
        ""            ;; no 70
        "Previous Setting"           ;; no 71
        "Complete"                   ;; no 72
        "Skip to error"              ;; no 73
        "In/Out..."                  ;; no 74
        "Copy Last Input"            ;; no 75
        "Remove Last Output"         ;; no 76
        "Remove Last Action"         ;; no 77
        "Exchange keys"              ;; no 78
        "Customize..."                ;; no 79
        "Entering trust mode. Recover usual behaviour\n by entering /* Trust=Off */." ;; no 80
        "All" ;; no 81
        "Lock mode" ;; no 82
        "Trust mode" ;; no 83
        "Lock semi-mode %s" ;; no 84
        "Trust semi-mode %s" ;; no 85
))))

;; The only function used outside :

(defun gp-messager (msg-number)
  (eval (nth msg-number (assq gp-language gp-messages-list))))

;; pari-messages.el ends here.
