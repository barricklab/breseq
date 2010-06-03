;; $Id: pari.el 4608 2003-04-18 08:51:46Z ramare $
;; pari.el -- GP/PARI editing support package.

;; Major mode for editing GP scripts. It provides functions for editing
;; the code and evaluating it . See the documentation of gp-script-mode
;; and read the file pariemacs.txt.

;; The original pari.el was written by Annette Hoffman.
;; Modified by David Carlisle (JANET: carlisle@uk.ac.man.cs).
;; Modified by Karim Belabas (belabas@math.u-bordeaux.fr) for gp 2.xxx.
;; Modified by Olivier Ramare (ramare@gat.univ-lille1.fr).

;; Maintainer (01-March-2003): Olivier Ramare (ramare@agat.univ-lille1.fr).

;; KNOWN DEFICIENCIES:
;;  -- The fontify part may have troubles with `}'. A `}' followed by
;;     a newline indicates the end of a function-definition starting with
;;     `{'. Spaces, or tab are *not* allowed. So if you use `}' as a string
;;     DON'T have it followed by a newline.

;; This file is split in six parts :
;;   PART I : MAIN CONSTANTS (contains a macro).
;;            Some of them may have to be modified by the user.
;; PART  II : KEYMAPS AND OTHER VARIABLES
;;            including 'gp-define-locked-keys.
;; PART III : gp-mode AND gp-script-mode
;;            Also the gp-locked*
;; PART  IV : GENERAL FUNCTIONS
;;            Contains: HANDLING THE WINDOWS ...
;;                      THE GP PROCESS
;;                      META-COMMANDS
;;  PART VI : MENU-BAR
;;            Contains: MENU BUILDERS (contains 3 constants)
;;                      MENU-BAR ITEM USED IN GP-SCRIPT-MODE
;;                      MENU-BAR ITEM USED IN GP-MODE

;; Note: emacs version should be higher than 20.3

(provide 'pari)

(defgroup gp-indentation nil
"GP customization subgroup concerning indentation
and furthering of constructs"
:group 'gp :prefix "gp-")

(defgroup gp-shell nil
"GP customization subgroup specific to gp-shell-mode"
:group 'gp :prefix "gp-")

(defgroup gp-font-lock-and-completion nil
"GP customization subgroup concerning colors and completion"
:group 'gp :prefix "gp-")

(defgroup gp-miscellana nil
"GP customization subgroup dedicated to less important switches"
:group 'gp :prefix "gp-")

(eval-and-compile
;; The next variables are here to pacify the compiler !
;; Do *not* assign any value to them or they may override ....
(defvar gp-tutorial-requiredp)
(defvar gp-auto-indent);
(defvar block-comment-end);
(defvar block-comment-start);

(defvar gp-prompt-pattern
  "^\\([?>]\\) [\n\t ]*"
  "Regexp used to match gp prompts.
Can be set with gp-set-prompt (bound to M-\\ p)")
(require 'pari-conf)
(require 'pari-messages)   ;; Provides: functions: gp-messager.
(require 'pari-completion) ;; Provides: functions: gp-quit-cpl-edit.
;; The following file uses variable gp-c-array which is defined and
;; created by 'pari-completion:
(require 'pari-help)       ;; Provides: functions: gp-menu-quit.
(require 'pari-fontification)
  ;; Provides: variable:  gp-fontification-keywords
  ;;           functions: gp-update-fontification, gp-find-global-var
(require 'sli-tools))

(unless (fboundp 'gp-update-fontification)
  (defun gp-update-fontification nil nil))
(unless (boundp 'gp-fontification-keywords)
  (defvar gp-fontification-keywords nil nil))
;; The use of gp-find-global-var if protected by a fboundp.

(eval-and-compile
  (require 'imenu)
  (require 'backquote)) ; This file is used in macros.

;;--------------------------
;; PART I : MAIN CONSTANTS
;;--------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Individual users may want to re-set some of the variables in this section
;; in a gp-mode-hook in their .emacs file (see pariemacs.txt for examples).

(defcustom gp-stack-size 10000000
"Default stack size passed to gp."
:type 'integer   :group 'gp-shell)

(defcustom gp-prime-limit 500000
"Default prime limit passed to gp."
:type 'integer   :group 'gp-shell)

(defcustom gp-prompt-for-args nil
  "*A non-nil value makes M-x gp act like C-u M-x gp, 
ie prompt for the command line arguments."
:type 'boolean   :group 'gp-shell)

(defcustom gp-keep-PARI-buffer-when-quitting t
"T means what it says..."
:type 'boolean   :group 'gp-shell)

(defcustom gp-locked-modep t
  "t means you cannot write above the last prompt.
If you try to modify an earlier input, emacs will automatically copy
it at the bottom of your file."
:type 'boolean
:initialize 'custom-initialize-default ;if you use :set, you should specify :initialize!
:set (lambda (sym val) (setq gp-locked-modep val) (gp-define-locked-keys))
:group 'gp-miscellana)

(defcustom gp-tutorial-requiredp t
"T if comments should be given for some functions."
:type 'boolean    :group 'gp-miscellana)
;; The functions concerned are : 'gp-make-cpl-file

(defcustom gp-menu-barp t
"A nil value means that we do not want any menu-bar"
:type 'boolean    :group 'gp-miscellana)

(defcustom gp-separate-window-for-mistakes nil
"T means errors under the gp calculator will be
displayed on a separate window."
:type 'boolean    :group 'gp-miscellana)

(defcustom gp-worryp t
"In gp-mode, finding \"input\" sets trust mode automatically,
except if this value is nil."
:type 'boolean :group 'gp-miscellana)

(defconst gp-temp-directory "/tmp/"
  "*Directory in which to create temporary files.")

(defvar gp-temp-file
  (expand-file-name (concat gp-temp-directory (make-temp-name "gp_#")))
  "Temporary file name used for text being sent as input to GP.")

(defvar gp-el-temp-file
  (expand-file-name (concat gp-temp-directory (make-temp-name "gp_#.el")))
  "Temporary file name used for text being sent as input to emacs.")

(defconst gp-max-saved-wind-conf 30
  "Maximal number of saved window configurations")

;;----------------------------------------
;; PART  II : KEYMAPS AND OTHER VARIABLES
;;----------------------------------------

(defvar gp-input-filter-hook nil
  "Hook run in `gp-input-filter'.")

(defvar gp-process nil "t if a GP process is running.")

(defvar gp-input-start nil
"Beginning of the expression to be send to GP. See `gp-copy-input'.")
(defvar gp-input-end nil
"End of the expression to be send to GP. See `gp-copy-input'.")
(defvar gp-complete-expression nil
"t if expression to be send to GP is complete. See `gp-copy-input'.")
(defvar gp-input-start-bracketp nil
"t if expression to be send to GP starts with a {.")
(defvar gp-reads-this-buffer nil
"name of the buffer gp is interpreting.")
(defvar gp-latest-error nil
"Regexp matching latest execution error. It contains a grouping
whose closing parenthesis corresponds to the point where gp
has detected a mistake.")
(defvar gp-registers-list nil
"List of registers from 0 to (1- gp-max-saved-wind-conf)
where window-configurations are stored.
See `gp-store-wind-conf' and `gp-restore-wind-conf'.")
(defvar gp-should-wait-for-outputp t
"t if gp should wait for output and fontify it
in `gp-send-input'. Automatically reset to t after each
input. See also `gp-input-filter'.")
(defvar gp-trust-mode nil
"nil is usual value.
If set to t, then the user is on its own, which means:
anything between (process-mark) and \\n is send to gp.
See also `gp-worryp'.")

(defconst gp-separator (list "----------") "")

(defconst gp-letters-list
  (string-to-list "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_=+-*/|^:!#()[]{}~%$,;.&?'`<> \"\\")
"See `gp-define-locked-keys'.")

(defvar gp-syntax-table nil
  "Syntax table in use in gp-mode and gp-script-mode buffers.")

(when (null gp-syntax-table)
  (setq gp-syntax-table (make-syntax-table))
  (mapcar (lambda (acons) (modify-syntax-entry (car acons) (cdr acons) gp-syntax-table))
          '((?( . "()") (?) .  ")(") (?[ . "(]") (?] . ")[") (?{ . "(}") (?} . "){") ; parenthesis
            (?# . ".") (?~ . "_") (?! . "_") (?% . "_")                          ; symbol constituent
            (?\\ . ". 12b") (?/ . ". 14") (?* . ". 23")            ; comments
            (?> . "." ) (?| . "." ) (?+ . ".") (?- . ".") (?= . ".") (?< . "." ) ; ponctuation
            (?. . "w") (?' . "w") (?$ . "w") (?_ . "w")))                        ; word constituent

  (if (string-match "XEmacs" emacs-version)
      (progn
	(modify-syntax-entry ?\n ">b"  gp-syntax-table)
	;; Give CR the same syntax as newline, for selective-display
	(modify-syntax-entry ?\^m ">b" gp-syntax-table))
    (modify-syntax-entry ?\n "> b"  gp-syntax-table)
    ;; Give CR the same syntax as newline, for selective-display
    (modify-syntax-entry ?\^m "> b" gp-syntax-table)))

(defvar gp-map nil
  "Local keymap used in buffer *PARI*.")

(defun gp-define-locked-keys nil
  (mapcar
    (lambda (achar)
      (define-key gp-map (vector achar)
                  (if gp-locked-modep
                      'gp-locked-self-insert-command
                      'self-insert-command)))
    gp-letters-list)
  (if gp-locked-modep
      (progn
        (define-key gp-map [mouse-2] 'gp-locked-mouse-2)
        (define-key gp-map "\C-?"    'gp-locked-backward-delete-char-untabify)
        (define-key gp-map "\C-d"    'gp-locked-delete-char)
        (define-key gp-map "\C-k"    'gp-locked-kill-line)
        (define-key gp-map "\C-y"    'gp-locked-yank))
      (define-key gp-map [mouse-2] 'mouse-yank-at-click)
      (define-key gp-map "\C-?"    'backward-delete-char-untabify)
      (define-key gp-map "\C-d"    'delete-char)
      (define-key gp-map "\C-k"    'kill-line)
      (define-key gp-map "\C-y"    'yank)))

(when (null gp-map)
(let ((map (make-sparse-keymap)))
(define-key map "\C-m"    (function gp-send-local-input))
(define-key map "\M-c"    (function gp-copy-input))
(define-key map "\M-\C-m" (function gp-C-j))
(define-key map "\C-j"    (function gp-C-j))
(define-key map "\C-c"    (function gp-interrupt))
(define-key map "\M-\\\\" (function gp-break-long-line))
(define-key map "\M-\\a"  (function gp-meta-a))
(define-key map "\M-\\b"  (function gp-meta-b))
(define-key map "\M-\\d"  (function gp-meta-d))
(define-key map "\M-\\m"  (function gp-meta-m))
(define-key map "\M-\\p"  (function gp-set-prompt))
(define-key map "\M-\\q"  (function gp-meta-q))
(define-key map "\M-\\r"  (function gp-meta-r))
(define-key map "\M-\\s"  (function gp-meta-s))
(define-key map "\M-\\t"  (function gp-meta-t))
(define-key map "\M-\\v"  (function gp-meta-v))
(define-key map "\M-\\w"  (function gp-meta-w))
(define-key map "\M-\\x"  (function gp-meta-x))
(define-key map "\C-a"    (function gp-beginning-of-line))
(define-key map [kp-home] (function gp-beginning-of-line))
(define-key map [home]    (function gp-beginning-of-line))
(define-key map "\C-p"    (function previous-line))
(define-key map "\C-n"    (function next-line))
(define-key map "\M-p"    (function gp-previous-cmd))
(define-key map "\M-n"    (function gp-next-cmd))
(define-key map "\M-s"    (function gp-skip-to-error))
(define-key map [C-kp-subtract] (function gp-remove-last-output))
(define-key map [M-kp-subtract] (function gp-remove-last-action))
(setq gp-map map)
(gp-define-locked-keys)))

(defvar gp-script-map nil
  "Local keymap used in gp-script-mode.")

(when (null gp-script-map)
(let ((map (make-sparse-keymap)))
(define-key map "\C-c\C-e"  (function sli-maid))
(define-key map "\C-c\C-f"  (function sli-tutor))
(define-key map "\M-\\\\"   (function gp-break-long-line))
(define-key map "\M-\\d"    (function gp-meta-d))
(define-key map "\M-\\t"    (function gp-meta-t))
(define-key map "\M-\\v"    (function gp-meta-v))
(define-key map "\M-\\z"    (function gp-run-in-region))
(define-key map "\M-s"      (function gp-skip-to-error))
(define-key map "\C-c\C-c"  (function gp-run-gp))
(setq gp-script-map map)))

;; Global keys. They *should* be global.

(define-key esc-map "o" (function gp-restore-wind-conf))

(define-key completion-list-mode-map [mouse-2] (function gp-mouse-2))

;; Maps used for the menu-bar.

(defvar GP-menu-map nil
"Keymap used for the menu-bar item GP in gp-mode")

(defvar GP-script-menu-map nil
"Keymap used for the menu-bar item GP in gp-script-mode")

;;---------------------------------------
;; PART ??? : sli-tools
;;---------------------------------------

(defcustom gp-tab-always-indent t
"Non-nil means TAB in MuPAD-mode should always reindent the current line,
regardless of where in the line point is when the TAB command is used."
:type 'boolean :group 'gp-indentation)

(defcustom gp-indent-level 3
"Indentation used after \"{\"."
:type 'integer :group 'gp-indentation)

(defun gp-set-and-recompute-indentation (sym val)
  (set sym val)
  (save-current-buffer
   (mapcar 
    (lambda (bf)
      (set-buffer bf)
      (when (eq major-mode 'gp-script-mode)
        (gp-learns-indentation)))
    (buffer-list))))

(defcustom gp-structures
  '((["for(" head 3] [")" end])
    (["return(" head 3] [")" end])
    (["(" head 1] [")" end])
    (["[" head 1] ["]" end])
    (["{" head gp-indent-level] ["}" end])
    ;(["{" head 0] ["local" strong 0] ["}" end])
    (["=" math-relation 1]) ;that's the last item of any relation, like in '=='
    (["<" math-relation 1])
    ([">" math-relation 1])
   )
"See `sli-structures'."
:type '(repeat (repeat (restricted-sexp :match-alternatives (vectorp listp))))
:initialize 'custom-initialize-default
:set 'gp-set-and-recompute-indentation
:group 'gp-indentation)

(defcustom gp-shift-alist '()
"See `sli-shift-alist'."
:type '(repeat (cons (vector string string) sexp))
:initialize 'custom-initialize-default
:set 'gp-set-and-recompute-indentation
:group 'gp-indentation)

(defcustom gp-no-heredity-list '()
"See `sli-no-heredity-alist'."
:type '(repeat (cons (vector string string) sexp))
:initialize 'custom-initialize-default
:set 'gp-set-and-recompute-indentation
:group 'gp-indentation)

(defvar gp-separators '(";" ",")
"See `sli-separators'.")

(defcustom gp-fixed-keys-alist '(("{" . 0))
;'(("local" . gp-indent-level) ("}" . 0))
"See `sli-fixed-keys-alist'."
:type '(repeat (cons string sexp))
:initialize 'custom-initialize-default
:set 'gp-set-and-recompute-indentation
:group 'gp-indentation)

(defcustom gp-keys-with-newline '(";")
"See `sli-keys-with-newline'."
:type '(repeat string)
:initialize 'custom-initialize-default
:set 'gp-set-and-recompute-indentation
:group 'gp-indentation)

(defcustom gp-add-to-key-alist '()
"See `sli-add-to-key-alist'."
:type '(repeat (cons string string))
:initialize 'custom-initialize-default
:set 'gp-set-and-recompute-indentation
:group 'gp-indentation)

(defcustom gp-more-maidp t
"Set it to nil if do not want `sli-maid'
to use `gp-add-to-key-alist'. Thus
if so 'end_proc' will not be followed by
a ':' and so on. See `sli-more-maidp'."
:type 'boolean
:initialize 'custom-initialize-default
:set 'gp-set-and-recompute-indentation
:group 'gp-indentation)

;;---------------------------------------
;; PART III : gp-mode AND gp-script-mode
;;---------------------------------------

(defun file-really-exists-p (file)
  (and (not (string= file "")) (file-exists-p file)))

(defun gp-kill-buffer-safely (abuffer)
  (let ((b (get-buffer abuffer)))
       (if b (kill-buffer b))))

(defun gp-learn-sexp nil
  "To teach emacs some elements of gp-syntax."
  ;; Treat comments as white spaces in sexp:
  (make-local-variable 'parse-sexp-ignore-comments)
  (setq parse-sexp-ignore-comments t)
  ;; Care about capital or not (always local):
  (setq case-fold-search nil)
  ;; Comments in sexp (We handle only one kind of comments):
  (make-local-variable 'comment-start)
  (setq comment-start "\\\\")  ;; A *string*, NOT a regexp.
  (make-local-variable 'comment-end)
  (setq comment-end "")
  (make-local-variable 'comment-start-skip)
  (setq comment-start-skip "\\\\\\\\.*$\\|/\\*\\([^\\*]\\|\\*[^/]\\)*\\*/"))

(defun pari-mode ()
  "Common part of 'gp-mode and 'gp-script-mode"
  (gp-learn-sexp)
  (set-syntax-table gp-syntax-table))

(defun gp-learns-indentation nil
  (require 'sli-tools)
  (sli-tools gp-structures gp-shift-alist gp-separators
             'sli-is-a-separatorp-default
             gp-fixed-keys-alist
             "^/\\*--+--\\*/\\|^}[ \t]*"
             gp-keys-with-newline nil gp-add-to-key-alist
             '("//" "\\\\") gp-no-heredity-list)
  (setq sli-more-maidp gp-more-maidp
        sli-tab-always-indent gp-tab-always-indent))

;; The line ";;;###autoload" is useless.
;; It will be useful when pari.el will be part
;; of the usual distribution of emacs.
;;;###autoload
(defun gp-script-mode nil
  "Major mode for editing GP input files.

The following bindings are available:
\\{gp-script-map}"

  (interactive)
  (kill-all-local-variables) ; exit from previous mode
  (setq major-mode 'gp-script-mode mode-name "GP script")
  (run-hooks 'pari-mode-hook)
  (run-hooks 'gp-script-mode-hook) ; Set up user preferences.
  (pari-mode)
  ; buffer-local:
  (setq imenu-generic-expression
        '((nil "^[{\t ]*\\([a-zA-Z]\\w*\\)(\\([^)]*\\)) *=[^=]" 1))
        imenu-case-fold-search nil)
  (setq block-comment-start "/*" block-comment-end "*/")
  (set (make-local-variable 'comment-indent-function) 'gp-indent-comment)
  (gp-learns-indentation)

  (gp-update-fontification)
  (use-local-map gp-script-map) ; Make gp-script-map the local map in this mode.
  (gp-add-imenu-index)
  (gp-init-script-menu-bar)         ; Start menu-bar.
  )

;; The line ";;;###autoload" is useless.
;; It will be useful when pari.el will be part
;; of the usual distribution of emacs.
;;;###autoload
(defun gp-mode ()
  "Major mode for running a gp-process.

The following bindings are available:
\\{gp-map}"

  (interactive)
  (kill-all-local-variables) ; exit from previous mode
  (setq major-mode 'gp-mode mode-name "GP")
  (run-hooks 'pari-mode-hook)
  (run-hooks 'gp-mode-hook) ; Set up user preferences.
  (pari-mode)
  ; buffer-local:
  (setq imenu-generic-expression
        '((nil "\\<\\([a-zA-Z]\\w*\\)(\\([^)]*\\)) *=[^=]" 1))
        imenu-case-fold-search nil)

  (gp-update-fontification)
  (use-local-map gp-map)    ; Make gp-map the local map of buffer *PARI*.
  (gp-add-imenu-index)
  (gp-init-menu-bar)        ; Start menu-bar.
  )

(defun gp-add-imenu-index nil
   (if (and gp-menu-barp
            (progn (require 'easymenu) (featurep 'easymenu)))
       (imenu-add-to-menubar "GP-functions")))

(defun gp-clear-temp-files nil
  "Remove temporary files that may have been created" 
   (if (file-exists-p gp-temp-file)
       (progn (delete-file gp-temp-file) 
              (message (gp-messager 2) gp-temp-file)))
   (if (file-exists-p gp-el-temp-file) 
       (progn (delete-file gp-el-temp-file) 
              (message (gp-messager 2) gp-el-temp-file))))

(defun gp-save-setting-kill-emacs nil
  "Remove temporary files."
  (gp-clear-temp-files))

(add-hook 'kill-emacs-hook (function gp-save-setting-kill-emacs))

(defun gp-displace-input nil
  (if (and (save-excursion (re-search-forward gp-prompt-pattern nil t))
           (save-excursion (re-search-backward gp-prompt-pattern nil t)))
      (let ((where (point)))
        ;(re-search-backward gp-prompt-pattern nil t)
        (goto-char (match-end 0))
        (setq where (- where (point)))
        (gp-copy-input)
        (re-search-backward gp-prompt-pattern nil t)
        (goto-char (+ where (match-end 0))))))

(defun gp-beginning-of-line nil
  (interactive)
  (beginning-of-line)
  (when (looking-at gp-prompt-pattern)
    (goto-char (match-end 0))))

(defun gp-locked-self-insert-command nil
  (interactive)
  (gp-displace-input)
  (insert-char last-command-char 1))

(defun gp-locked-mouse-2 (anevent)
  (interactive "e")
  (mouse-set-point anevent)
  (gp-displace-input)
  (yank))

(defun gp-locked-yank nil
  (interactive)
  (gp-displace-input)
  (yank))

(defun gp-locked-backward-delete-char-untabify nil
  (interactive)
  (gp-displace-input)
  (backward-delete-char-untabify 1))

(defun gp-locked-kill-line nil
  (interactive)
  (gp-displace-input)
  (kill-line))

(defun gp-locked-delete-char nil
  (interactive)
  (gp-displace-input)
  (delete-char 1))

;;-----------------------------
;; PART IV : GENERAL FUNCTIONS
;;-----------------------------

;;--------------------------
;; HANDLING THE WINDOWS ...
;;--------------------------
;; At the beginning, the user has asked for one window, but s/he may well
;; have introduced another window in-between (or even several ones).
;; We should then use only one other fixed window for everything else.
;; But since the list of the buffers displayed in a window does not exist,
;; and since the user may well change of window by ITself, we can't do much.

(defun gp-depile-wind-conf nil (setq gp-registers-list (cdr gp-registers-list)))

(defun gp-backward-wind-conf nil
  "Restore previously stored window configuration."
 (if (not (equal gp-registers-list nil))
     (progn
       (jump-to-register (car gp-registers-list))
       (setq gp-registers-list (cdr gp-registers-list)))))

(defun gp-store-wind-conf nil
  "Add a the current window configuration to the pile. If the pile
has more than 'gp-max-saved-wind-conf items
(0,1,...,(1- gp-max-saved-wind-conf)) then the first item is lost."
  (if (= (length gp-registers-list) gp-max-saved-wind-conf)
      (setq gp-registers-list (nreverse (cdr (nreverse gp-registers-list)))))
  (let ((next (if (equal gp-registers-list nil) 0
                  (if (= (car gp-registers-list) (1- gp-max-saved-wind-conf)) 0
                      (1+ (car gp-registers-list))))))
       (window-configuration-to-register next)
       (setq gp-registers-list (cons next gp-registers-list))))

(defun gp-restore-wind-conf (&optional arg)
  "Restore the previous window-configuration, killing the *gp-help* buffer
if it was and is no more displayed. When called with prefix C-u, end the
edition of the completion-file (if any were edited)."
  (interactive "P")
  (if (and arg (= (car arg) 4)) ;; Meaning that the call has been C-u M-o
      (gp-quit-cpl-edit)
      (let ((had-help-windowp (and (get-buffer "*gp-help*")
                                   (get-buffer-window "*gp-help*")))
            (had-message-windowp (and (get-buffer "*gp-messages*")
                                      (get-buffer-window "*gp-messages*"))))
           (gp-backward-wind-conf)
           ;; Kill the buffer *gp-help* if it is not displayed anymore:
           (if had-help-windowp
             (if (not (get-buffer-window "*gp-help*"))
                 (kill-buffer "*gp-help*")))
           (if had-message-windowp
             (if (not (get-buffer-window "*gp-messages*"))
                 (kill-buffer "*gp-messages*"))))
      ;; When called from menu-bar, write nothing in the minibuffer:
      (message "")))

(defun gp-info-wind-conf nil (message (gp-messager 4)))

(defun buffer-visiblep (abuffer-name)
  (if (get-buffer-window abuffer-name) t nil))

(defun gp-pgrmp (abuffer)
  "Set buffer ABUFFER and return t if ABUFFER is in gp-script-mode."
  (set-buffer abuffer) (eq major-mode 'gp-script-mode))

(defun gp-possible-file-name nil
  "Try to guess the name of a likely gp-program."
  ;; First tries the existing windows, then the existing buffers.
  (let ((pgrm nil))
       (walk-windows
         (lambda (wind)
          (if (gp-pgrmp (window-buffer wind))
              (setq pgrm
                 (cons (buffer-name (window-buffer wind)) pgrm)))))
       (if pgrm (car pgrm) ;Return value if a window is displaying
                           ;a candidate gp-program.
           (mapcar
             (lambda (abuffer)
              (if (gp-pgrmp abuffer)
                  (setq pgrm (cons (buffer-name abuffer) pgrm))))
             (buffer-list))
           (if pgrm (car pgrm) ;Return value if a buffer is a candidate gp-program.
                    nil        ;Return value if fail.
   ))))

(defun gp-window-manager (my-buffer-name option)
"Takes care of the windows in gp-mode and gp-script-mode.
Displays the buffer MY-BUFFER-NAME in a proper window.
The variable OPTION is
  -- gp-beginning when we handle the beginning of a procedure. If a buffer
                  already exists with this name, only store the wind-conf.
  -- gp-beginning-temp when we handle the beginning of a procedure. If a
                       buffer already exists with this name, store it.
  -- gp-remove-help-now to remove help-window,
  -- gp-remove-help-old-config to wait and remove help-window without
                               touching to the other windows.
  -- gp-remove-help-now-old-config to remove help-window without
                               touching to the other windows.
  -- gp-show-help which is similar to gp-beginning for the help buffer
                  except that we do not erase the content of this buffer.
  -- nil when it is the end of a call.
The variable MY-BUFFER-NAME is one of 
\"*PARI*\"  \"*gp-help*\" \"*gp-menu*\". "

  (cond ((and (string= my-buffer-name "*PARI*")
              (eq option 'gp-beginning)
              (get-buffer-window "*PARI*"))
         ;; We go to *PARI* and a window already exists with this buffer.
         (gp-store-wind-conf)
         (select-window (get-buffer-window "*PARI*")))
        
        ((and (string= my-buffer-name "*PARI*")
              (eq option 'gp-beginning)
              (not (get-buffer-window "*PARI*")))
         ;; We go to *PARI* and a window doesn't exist with this buffer.
         (if (= (count-windows) 1)
             ;; If there is only 1 window containing anything but *scratch*
             ;; split the window in 2, else use this window:
             (progn (if (not (string= (buffer-name) "*scratch*"))
                        (select-window (split-window-vertically)))
                    (switch-to-buffer "*PARI*"))
             ;; At least two windows exist. Do not create another one
             ;; and first try to use the help window, else the
             ;; starting window.
             (gp-store-wind-conf)
             (cond ((get-buffer-window "*gp-help*")
                    (select-window (get-buffer-window "*gp-help*"))
                    (switch-to-buffer "*PARI*"))
                   (t (switch-to-buffer-other-window "*PARI*")))))

        ((and (string= my-buffer-name "*PARI*")
              (not option)
              (get-buffer "*PARI*"))
         ;; We want to exit from *PARI*.
         (if (> (count-windows) 1)
             (delete-windows-on "*PARI*")
             ;; Else only one window.
             (if (string= (buffer-name (window-buffer)) "*PARI*")
                 ;; This only window displays "*PARI*"
                 (let ((next-buffer (gp-possible-file-name)))
                      (if next-buffer (switch-to-buffer next-buffer)
                          ;; Else, don't know what to do !
                          (gp-restore-wind-conf)
                          ))))
         (unless gp-keep-PARI-buffer-when-quitting
           (with-current-buffer (get-buffer "*PARI*")
             (let ((inhibit-read-only t))
               (remove-text-properties (point-min) (point-max) '(read-only nil))))
           (kill-buffer "*PARI*")))

        ((and (get-buffer my-buffer-name)
              (member my-buffer-name '("*gp-help*" "*gp-menu*"))
              (eq option 'gp-remove-help-now)
              (get-buffer-window my-buffer-name))
         ;; A buffer displaying "*gp-help*" or "*gp-menu*" exists.
         ;; We want to remove the message.
         (if (or (string= my-buffer-name "*gp-help*")
                 (not (get-buffer "*gp-help*")))
             ;; Exit from help or the gp-menu is alone:
             (gp-restore-wind-conf)
             (if (string= my-buffer-name "*gp-menu*")
             ;; The previous condition should always be verified!
             ;; We should remove the window displaying gp-menu:
                 (progn
                   (if (and (= (count-windows) 2)
                            (get-buffer "*gp-help*"))
                       (progn
                         (gp-depile-wind-conf)
                         (switch-to-buffer "*gp-help*")
                         (other-window 1))
                       (gp-restore-wind-conf)))))
         ;; We have to kill the buffer (in any case) and select
         ;; a proper buffer for this window in case this killing
         ;; made something weird appear:
         (gp-kill-buffer-safely my-buffer-name)
         ;; since it may have been destroyed by 'gp-restore-wind-conf.
         (let ((buffer-to-select ""))
              (save-excursion
               (let ((abufferlist (buffer-list)))
                    (while (and (string= buffer-to-select "")
                                abufferlist)
                      (set-buffer (car abufferlist))
                      (if (memq major-mode '(gp-script-mode gp-mode))
                          (setq buffer-to-select (buffer-name)))
                      (setq abufferlist (cdr abufferlist)))))
              ;; Last weird case to handle: the buffer we have selected
              ;; is already being shown on another window.
              ;; Then kill our window.
              (if nil ;(buffer-visiblep buffer-to-select)
                  (delete-window)
                  (or (string= buffer-to-select "") ;; Let it be !
                      (switch-to-buffer buffer-to-select)))))

        ((and (get-buffer my-buffer-name)
              (member my-buffer-name '("*gp-help*" "*gp-menu*"))
              (memq option '(gp-remove-help-old-config
                             gp-remove-help-now-old-config)))
         ;; A buffer displaying "*gp-help*" or gp-menu exists.
         ;; We want to remove the message without touching
         ;; to the window-configuration.
         (cond ((eq option 'gp-remove-help-old-config)
                (message (gp-messager 5))
                (read-event)))
         (kill-buffer my-buffer-name))

        ((and (string= my-buffer-name "*gp-help*")
              (memq option '(gp-beginning gp-show-help))
              (get-buffer-window "*gp-help*"))
         ;; We go to *gp-help* and a window already exists with this buffer.
         (select-window (get-buffer-window "*gp-help*"))
         (or (eq option 'gp-show-help) (erase-buffer)))

        ((and (string= my-buffer-name "*gp-help*")
              (eq option 'gp-beginning-temp)
              (get-buffer-window "*gp-help*"))
         ;; We go temporarily to *gp-help* and a window already exists with
         ;; this buffer.
         (gp-store-wind-conf)
         (select-window (get-buffer-window "*gp-help*"))
         (erase-buffer))

        ((and (get-buffer my-buffer-name)
              (member my-buffer-name '("*gp-help*" "*gp-menu*"))
              (eq option 'gp-remove-help-now))
         ;; Since it got here, my-buffer-name is not displayed.
         (gp-kill-buffer-safely my-buffer-name))

        ((and (string= my-buffer-name "*gp-help*")
              (memq option '(gp-beginning gp-beginning-temp gp-show-help))
              (not (get-buffer-window "*gp-help*")))
         ;; We go to *gp-help* and a window doesn't exist with this buffer.
         (gp-store-wind-conf)
         (if (= (count-windows) 1)
             (progn (select-window (split-window-vertically))
                    (switch-to-buffer "*gp-help*"))
             (cond ((and (get-buffer-window "*PARI*")
                     (not (eq (get-buffer-window "*PARI*") (selected-window))))
                    (select-window (get-buffer-window "*PARI*"))
                    (switch-to-buffer "*gp-help*"))
                   (t (switch-to-buffer-other-window "*gp-help*"))))
         (or (eq option 'gp-show-help) (erase-buffer)))

        ((and (string= my-buffer-name "*gp-menu*")
              (eq option 'gp-beginning))
         ;; We go to gp-menu.
              (if (get-buffer "*gp-menu*")
                  ;; A gp-menu already exists. Kill it first:
                  (save-excursion
                    (set-buffer "*gp-menu*")
                    (gp-menu-quit)))
              (gp-store-wind-conf)
              (if (get-buffer-window "*gp-help*")
                  (progn
                    (select-window (get-buffer-window "*gp-help*"))
                    (switch-to-buffer
                       (get-buffer-create "*gp-menu*"))
                    (kill-buffer "*gp-help*"))
                  (if (= (count-windows) 1)
                      (split-window-vertically))
                  (switch-to-buffer-other-window
                    (get-buffer-create "*gp-menu*"))))
        ))  ; end of 'gp-window-manager

;;----------------
;; THE GP PROCESS
;;----------------

(defun gp-make-gp-prompt-pattern (a-pattern)
  "Add regexp a-pattern at beginning of line followed by any
amount of space/tab/newline to gp-prompt-pattern."
;; gp-prompt-pattern matches:
;; (New prompt plus any following white space) OR (Old pattern).
  (let ((aux (concat "^\\(" a-pattern "\\)[\n\t ]*")))
    (setq gp-prompt-pattern (concat aux "\\|" gp-prompt-pattern)
          gp-fontification-keywords
           (append
             (list aux (list 1 'gp-prompt 'prepend t))
             gp-fontification-keywords
             ))))

(defun gp-beginning-of-last-line nil
  (goto-char (point-max))
  (re-search-backward gp-prompt-pattern)
  (goto-char (match-end 0)))

(defun gp-stiffen-prompt nil
  (save-excursion
    (re-search-backward gp-prompt-pattern nil t) ; should be beginning of line ...
    (let ((inhibit-read-only t)) ; in case of "quit" command.
      (put-text-property (1- (match-end 0)) (match-end 0) 'rear-nonsticky t)
      (put-text-property (1- (match-beginning 0)) (match-end 0) 'read-only t))))

(defun gp-wait-for-output (point-init &optional nomessage process nostiff)
"Hang around until the prompt appears.
PROCESS defaults to gp-process."
  (let ((notdone t))
  (or process (setq process gp-process))
  (setq nostiff (or nostiff (not (eq process gp-process))))
  (while notdone
    ;; Wait till something comes out:
    (while (and (not (accept-process-output process 0 300))
                (not (= point-init (point)))
                ;; Following line is required for the \q command:
                (eq 'run (process-status process))))
    (let ((p (point)))
      (if (or
            ;; Following lines are required for the \q command:
	    (not (and (processp process) 
		      (eq 'run (process-status process))))
            (save-excursion
              (if (re-search-backward gp-prompt-pattern nil t)
                  (= (match-end 0) (point-max))
                  nil)))
  ;; If gp is not running, or the prompt has appeared, stop.
	(progn (or nomessage (message (gp-messager 6)))
               (setq notdone nil))
  ;; Else flush the buffer and wait a bit longer.
	(progn (or nomessage (message (gp-messager 7)))))
      (goto-char p))))
    (sit-for 0)
    (goto-char (point-max))
    (unless nostiff (gp-stiffen-prompt)) ;(print "Out of stiffening !")
    (set-marker (process-mark process) (point)))

(defun gp-get-shell (process-name process-buffer-name cmd)
    "Explicit. Distinguishes bash/sh and [t]csh. Aimed at command gp+parameters."
  ;; We put the number of lines to 1000 so that no break will
  ;; occur when giving long comment like with "?6". We do not
  ;; want any "Return to continue", the editing job should
  ;; be done by emacs and not by gp.
  (if (member (file-name-nondirectory shell-file-name) '("bash" "sh"))
      (start-process process-name process-buffer-name
                 shell-file-name "-c"
                 (concat "(stty -echo nl; TERM=emacs; LINES=1000; PAGER=cat; COLUMNS="
                         (number-to-string (window-width))
                         "; export TERM COLUMNS LINES; " cmd ")"))
      (start-process process-name process-buffer-name
                 shell-file-name "-c"
                 (concat "stty -echo nl; env TERM=emacs PAGER=cat LINES=1000 COLUMNS="
                      (number-to-string (window-width)) " "
                    cmd))))

(defun gp-background nil
  "Same as 'gp except that it doesn't switch to the buffer `*PARI*'.
The answer is t if success, and nil otherwise."
 (save-excursion
  (if (and (processp gp-process)
           (eq 'run (process-status gp-process)))
    t ; If gp is already running, do nothing.

;; Else start up gp in the buffer.

    ;; Create the buffer `*PARI*' if required.
    (set-buffer (get-buffer-create "*PARI*"))
    (unless gp-keep-PARI-buffer-when-quitting
      (erase-buffer))
    (run-hooks 'pari-mode-hook 'gp-mode-hook)
;; Form the command line string.
    (let*((process-connection-type t) ; use PTY.
          (gp-cmd
           (concat
             gp-file-name " -s " (number-to-string gp-stack-size)
                          " -p " (number-to-string gp-prime-limit)
	     " --emacs"  ; --emacs requested by gp2.
             )))
 
;; Insert the command line string into the *PARI* buffer (for reference)
      (insert (format (gp-messager 41) gp-cmd))
;; Start gp.
      (setq gp-process (gp-get-shell "pari" "*PARI*" gp-cmd))
;; Clean up when the gp process has finished.
    (set-process-sentinel gp-process (function gp-sentinel)))
    ;; We should run the hook as the prompt may have
    ;; been changed in the .gprc:
    (run-hooks 'pari-mode-hook)
    (gp-wait-for-output (point-min))
    (setq gp-input-start (point) gp-input-end (point))
    ;; Introduce 'gp-mode
    ;; (Should be here as the prompt needs a gp-session running,
    ;; as well as the choice readline on/off):
    (unless (eq major-mode 'gp-mode) (gp-mode))
    (setq mode-line-process '(": %s"))
    (if (memq (process-status gp-process) '(signal exit))
        (setq gp-process nil) t))))

(defun gp nil
  "
   Open a buffer and a window for the execution of gp.

   The following bindings are available:
   \\{gp-map}

  The variables
  gp-file-name gp-stack-size gp-prime-limit
  determine the command line that starts gp."

  (interactive)
  (if (gp-background)
      (progn
        (gp-window-manager "*PARI*" 'gp-beginning)
        ;; Hilight first prompt:
        (goto-char (point-max))
        (gp-update-fontification))
      (message (gp-messager 8))))

(defun gp-run-in-region (beg end)
  "Run GP on the current region.  A temporary file (gp-temp-file) is
written in gp-temp-directory, but GP is run in the current directory."
;; Set gp-input-start, gp-input-end and gp-reads-this-buffer.
   (interactive "r")
   (setq gp-input-start beg gp-input-end end)
   (setq gp-reads-this-buffer (buffer-name))
   (gp-input-filter)
   (write-region beg end gp-temp-file nil nil)
   (gp)     ;; In case a GP-process was not already running, starts one.
                ;; In any case, switches to buffer "*PARI*".
   (gp-beginning-of-last-line)
   (insert (concat "\\r " gp-temp-file))
   (set-marker (process-mark gp-process) (point))
   (gp-send-input))

(defun gp-read-input (prompt default sep flag)
  "If flag is non-nil, reads string (if string is \"\" uses default).
Else, if flag is nil, set string to default.
If resulting string is not \"\" prepends sep.
As a special case, if string is \" \", return \"\"."

  (let ((string
    (if flag
;; If flag is non-nil prompt for input from mini-buffer.
      (read-input
        (concat prompt " (Default "default") "))
;; Else use the default string.
        default)))

    (if (string-equal string "")
      (if (string-equal default "") 
         ""                     ;; If string and default both "": 
         (concat sep default))  ;; If string "" and default is non empty:
      (if (string-equal string " ")
        ""                      ;; If string is a space:
        (concat sep string))))) ;; If string is non empty:

(defun gp-sentinel (proc msg)
  "Sentinel for the gp-process in buffer *PARI*."

  (gp-kill-buffer-safely "*gp-menu*")
  (gp-window-manager "*gp-help*" 'gp-remove-help-now)
      ;; We do not kill the buffer "*Completions*" as it may have
      ;; been triggered by something else.
  (gp-window-manager "*PARI*" nil)
  (gp-clear-temp-files)
  (setq gp-process nil))

(defun gp-output-filter ()
  (let ((wind (selected-window))
        (errp (save-excursion
                (goto-char (1+ gp-input-end))
                (looking-at "^  \\*\\*\\*  \\|^Unknown function"))))
      (if errp
       	(progn
	  (let ((copy (buffer-substring-no-properties (1+ gp-input-end)
                       (progn
                         (goto-char (point-max)) ;; We should already be there!
                         ;; Remove last prompt line ...
	                 (beginning-of-line)
                         ;; and final empty lines:
                         (skip-chars-backward " \t\n")
                         (point)))))
            (delete-region gp-input-end (point-max))
            (gp-store-wind-conf)
            (other-window 1)
	    (split-window-vertically)
            ;(other-window 1)
            (switch-to-buffer (get-buffer-create "*gp-messages*"))
	    (erase-buffer)
	    (insert copy)
            (shrink-window-if-larger-than-buffer)
	    (goto-char (point-min))
            (gp-info-wind-conf)
	    (select-window wind))))))

(defun gp-special-output-filter nil
  (let ((errp (save-excursion
                (goto-char (1+ gp-input-end))
                (or (looking-at "^  \\*\\*\\*   unexpected character: \\.\\.\\.")
                    (looking-at "^  \\*\\*\\*   expected character: [^\n]*\n  \\*\\*\\*   instead of: ")
                    (looking-at "^  \\*\\*\\*   expected character: [^\n]*\n  \\*\\*\\*   instead of:\n  \\*\\*\\*   \\.\\.\\.")
                    (looking-at "^  \\*\\*\\*   unknown function or error in formal parameters:\n  \\*\\*\\*   \\.\\.\\.") 
                    (looking-at "^  \\*\\*\\*   unknown function or error in formal parameters: ")
                    (looking-at "^  \\*\\*\\*   unexpected character: "))
                )))
      (if errp  ;; T if an error has been detected.
       	(progn
          (goto-char (match-end 0))
	  (let* (;; the line containing the mistake:
                 (astring (buffer-substring-no-properties (point)
                                            (progn (end-of-line) (point))))
                 ;; how many characters of astring have been sent to
                 ;; gp-latest-error:
                 (place 1)
                 ;; "location" of the mistake:
                 (which-char (+ (length astring) (- (search-forward "^")
                                                    (progn (end-of-line) (point))))))
            ;; We create gp-latest-error:
           (setq gp-latest-error (concat "\\(" (regexp-quote (substring astring 0 1))))
            (while (< place (length astring))
              (setq gp-latest-error
                    (concat gp-latest-error "[ \t\n]*\\(\\(/\\*[^\\*]*\\*/\\|\\\\\\.*$\\)[ \t\n]*\\)*"
                            (regexp-quote (substring astring place (setq place (1+ place))))))
              (if (= place which-char)
                  (setq gp-latest-error (concat gp-latest-error "\\)"))))
            (select-window (get-buffer-window gp-reads-this-buffer))
            (goto-char (point-min))
            (gp-skip-to-error))))))

(defun gp-skip-to-error nil
  "Well, it essentially does not work... No it works ! but
succeeds only every other day..."
  (interactive)

  (if (and gp-reads-this-buffer gp-latest-error
           (buffer-live-p (get-buffer gp-reads-this-buffer)))
    (progn (print "Going")
      (if (string= (buffer-name) gp-reads-this-buffer) nil
          (switch-to-buffer gp-reads-this-buffer)
          (goto-char (point-min)))
      (if (re-search-forward gp-latest-error nil t)
          (progn (goto-char (1- (match-end 1)))
                 ;; Warn the user this place is maybe not the good one !:
                 (message (gp-messager 35))
                 ;; Make the cursor blink:
                 (let ((old-color (frame-parameter nil 'cursor-color))
                       ;; Does not work... Why ? :
                       (other-color (frame-parameter nil 'background-color))
                       (how-many 6) (how-long-dark 50) (how-long-light 70) aux)
                       (setq other-color "blue")
                       (while (> how-many 0)
                         (set-cursor-color other-color)
                         (sit-for 0 how-long-light)
                         (set-cursor-color old-color)
                         (sit-for 0 how-long-dark)
                         (setq how-many (1- how-many)))
                       (set-cursor-color other-color)
                       (sit-for 0 how-long-light)
                       (set-cursor-color old-color)))
          ;; Could not locate the error:
          (message (gp-messager 34))))
    (message (gp-messager 36))))

(defun gp-run-gp nil
  "Sends a file to be run under GP."
  ;; This command is simply a compositum of 'gp-usual-start
  ;; and 'gp-meta-r. However the default file is different.
   (interactive)
   (let* ((gp-pgrm (gp-read-input (gp-messager 69)
                                  (gp-possible-file-name) "" t)))
         (if (get-buffer gp-pgrm)
             (save-excursion
               (set-buffer gp-pgrm)
               (setq gp-reads-this-buffer gp-pgrm)
               (if (buffer-modified-p) (save-buffer))
               (setq gp-input-start (point-min)
                     gp-input-end (point-max))
               (gp-input-filter)
               ;; In case 'gp-input-filter modified the buffer:
               (setq gp-pgrm (buffer-file-name))
               (if (buffer-modified-p) (save-buffer 0))))
         (gp)  ;; In case a GP-process was not already running, starts one.
               ;; In any case, switches to buffer "*PARI*". 
         (gp-beginning-of-last-line)
         (insert (concat "\\r " gp-pgrm))
         (set-marker (process-mark gp-process) (point))
         (gp-send-input)))

(defun gp-C-j nil
  (interactive)
  (insert-char ?\n 1)
  (put-text-property (1- (point)) (point) 'gp-virtual-newline t))

(defun gp-is-virtual (where)
  (get-text-property where 'gp-virtual-newline))

(defun gp-end-of-inputp nil
  ;; Beware we do not impose initial point to be at end of line !!
  (save-excursion
    (forward-char -1)
    (and (not (and (looking-at "\n")
                   (gp-is-virtual (point))))
         (not (and (looking-at "\n")
                   (save-excursion
                     (forward-char -1)
                     (looking-at "\\\\"))))
         (not (looking-at "\\\\")))))

(defun gp-match-input (limit)
  (let (rep)
    (if (and (re-search-forward gp-prompt-pattern limit t)
             (setq rep (gp-find-end-of-input limit)))
        (progn (set-match-data (list (point) (goto-char rep)))
               t)
        nil)))

(defun gp-find-end-of-input (end)
  "Gives the position of next end-of-input and nil if none."
  (save-excursion
    (while (and (re-search-forward "\n" end t)
                (not (gp-end-of-inputp))))
    (if (and (char-equal (char-after (1- (point))) ?\n)
             (gp-end-of-inputp))
        (point)
        ;; No more newlines in sight:
        (goto-char end)
        (if (gp-end-of-inputp) (point) nil))))

(defun gp-copy-input (&optional nocontrol)
  "Copy expression around point to the end of the buffer.
(Unless this is already the last expression.)
If NOCONTROL is non nil, then 'gp-complete-expression is
automatically set to t and emacs will not check whether the
expression is complete or not."
  (interactive)
;; Go back to the end of prompt, and record that point.

  (re-search-backward gp-prompt-pattern nil t)
  (goto-char (setq gp-input-start (match-end 0)))  ;; end of prompt
  (setq gp-input-start-bracketp (looking-at "[ \t]*{"))

  (let ((lastp t))  ; t if this input is the last one
                    ; (i.e. is not followed by a prompt).
    (if gp-input-start-bracketp
      (progn
        (save-excursion
          (if (re-search-forward "}" nil t)
              (setq gp-input-end (point))
              (setq gp-input-end nil)))
        (setq lastp (not (re-search-forward gp-prompt-pattern nil t)))
        (if (or (and (not lastp) gp-input-end
                     (< (match-beginning 0) gp-input-end))
                (not gp-input-end))
            ;; Bad: not the last one but well backeted, except that the
            ;; "closing" bracket in on the other side of next prompt!
            ;; Or unfinished construct (no closing }):
            (progn
              (if lastp
                ;; Repair the \n:
                (progn
                  (goto-char gp-input-start)
                  (while (search-forward "\n" nil t)
                    (put-text-property (1- (point)) (point)
                                       'gp-virtual-newline t)))
                (setq gp-input-end (match-beginning 0)))
              ;; not the last one but badly bracketed
              (setq gp-input-start-bracketp nil)))))

    (if gp-input-start-bracketp ; properly enclosed expression.
      (setq gp-complete-expression t)

      (setq gp-input-end (gp-find-end-of-input (point-max)))
      (if gp-input-end
        (setq gp-complete-expression t)
        (goto-char (point-max))
        (setq gp-input-end (point-max)
              gp-complete-expression (gp-end-of-inputp)))

      (setq lastp (equal gp-input-end (point-max)))
      (if (not lastp)
        ;; It is not the last expression:
        (setq gp-input-end (1- gp-input-end)))
      ;; Remove trailing (virtual) \n :
      (if (char-equal (char-after (1- gp-input-end)) ?\n)
        (progn
           (goto-char gp-input-end)
           (re-search-backward "[^\n]\\(\n\\)" gp-input-start t)
           (goto-char (setq gp-input-end (match-beginning 1)))
           (put-text-property (point) (1+ (point)) 'gp-virtual-newline nil)
           (setq gp-complete-expression (gp-end-of-inputp))))

      ;; We refine 'gp-complete-expression:
      (let ((ans (parse-partial-sexp gp-input-start gp-input-end)) a-pt)
         (setq gp-complete-expression
              (or nocontrol
                (and gp-complete-expression
                     (equal (nth 0 ans) 0)  ; Depth in parens is 0.
                     (not (nth 3 ans))      ; Not inside a string.
                     (or (not (nth 4 ans))  ; Not inside a comment...
                         (nth 7 ans))       ; except if it starts with \\.
                     )))))

    (goto-char (point-max))
    (if (not lastp)
     ;; It is not the last expression:
     (progn
       (insert (buffer-substring gp-input-start gp-input-end))
       (if gp-complete-expression
         nil
         (ding)
         (message (gp-messager 9)))))))

(defun gp-input-filter nil
  "Look at buffer between gp-input-start and gp-input-end.
-- If it finds a string `default(prompt,foo)', and
foo is a gp-string, try to set gp-prompt-pattern
correctly. If foo is not a string, warn the user that
something wrong may happen.
-- If a line `\\@' is found, set variable 'gp-should-wait-for-ouputp
to nil.
-- If a comment `/* */' starts by a @, the content is understood
as a Lisp command and appended to the file gp-el-temp-file. This
file is empty at the beginning. This file is loaded before execution
of the gp program."
  ;; Follow 'gp-copy-input so the input has been copied at the end
  ;; of the buffer. 'gp-input-start and 'gp-input-end are set.
  (interactive)
  (save-excursion

   ;; Take care of `/*@ foo */':
   (goto-char gp-input-start)
   (let ((first-time t))
     (while (re-search-forward "/\\*@\\(\\([^\\*]\\|\\*[^/]\\)*\\)\\*/" gp-input-end t)
       (if first-time
         (progn (setq first-time nil)
                (if (file-exists-p gp-el-temp-file)
                    ;; Remove any older version:
                    (delete-file gp-el-temp-file))))
       ;; Append the Lisp part to the file "gp-prgm":
       (write-region (match-beginning 1) (match-end 1)  gp-el-temp-file t)
       (write-region "\n" nil gp-el-temp-file t)
       (goto-char (match-end 0)))
     (if first-time nil
         ;; Load the Lisp part:
         (load-file gp-el-temp-file))))

   ;; Run filter hooks if any. It should be here since the hook
   ;; may have been defined precisely in this file.
   ;; Should not be surrounded by a save-excursion !
   (run-hooks 'gp-input-filter-hook)

 (save-excursion
 (unless gp-trust-mode
   ;; Warn the user that `default(prompt,APROMPT)' may not work properly.
   (goto-char gp-input-start)
   (while (re-search-forward "default(prompt," gp-input-end t)
     ;; Try to set the prompt if it is a simple string.
     (goto-char (match-end 0))
     (let ((start (1+ (match-end 0))))
       (if (and (looking-at "[ ]*\"")
                (re-search-forward "\")" gp-input-end t))
           (gp-make-gp-prompt-pattern
              (gp-make-prompt-pattern
                (buffer-substring-no-properties start (- (match-end 0) 2))))
           ;; Else troubles...
           (message (gp-messager 10))
           (sit-for 2))))

   ;; Take care of `\\@':
   (goto-char gp-input-start)
   (if (re-search-forward "^\\\\\\\\@$" gp-input-end t)
       (setq gp-should-wait-for-outputp nil))
    
   ;; Take care of virtual-newlines:
   (goto-char gp-input-start)
   (while (re-search-forward "[^\\\\]\\(\n\\)" gp-input-end t)
     (if (gp-is-virtual (1- (point)))
         (progn
           (replace-match "\\\n" t t nil 1)
           (setq gp-input-end (1+ gp-input-end)))))
   )

   ;; Out of (unless gp-trust-mode ...
   ;; Take care of gp-trust-mode:
   (goto-char gp-input-start)
   (let ((case-fold-search t))
      (when (re-search-forward "/\\*\\s-*Trust\\s-*=\\s-*\\(On\\|Off\\)\\s-*\\*/"
             gp-input-end t)
            (setq gp-trust-mode (not (null (member (match-string-no-properties 1)
                                               '("on" "On" "oN" "ON")))))))

   ;; Tries to understand "input" ...
   (save-excursion
   (unless (or gp-trust-mode (not gp-worryp) (not (eq major-mode 'gp-mode)))
     (goto-char gp-input-start)
     (when (re-search-forward "\\<input\\>" gp-input-end t)
       (message (gp-messager 80))
       (setq gp-trust-mode t))))
   ))

(defun gp-treat-special-inputp nil
    (cond (gp-trust-mode (setq gp-input-start (marker-position (process-mark gp-process))
                           gp-input-end (point-max) gp-complete-expression t)
           (when (save-excursion (re-search-forward gp-prompt-pattern gp-input-end t))
              (gp-copy-input t)) t)
          ((save-excursion (beginning-of-line)
                           (looking-at (concat "\\(" gp-prompt-pattern "\\)\\?\\\\")))
           (save-excursion
             (end-of-line)
             (setq gp-input-start (- (point) 2) gp-input-end (point)
                   gp-complete-expression t)))
          (t nil)))

(defun gp-send-input (&optional localp)
  "Sends input to gp. Does not send incomplete expressions
ie those starting with {, without a matching }, or those
ending with \\ .
Uses a temporary file (and \\r ) for large expressions.
If LOCALP is non nil, then it is assumed the input comes
from the *PARI* buffer, in which case if this input was a
`\r '-command, sends the output to `gp-output-filter'.
If LOCALP is nil, then if a file is being read which is
currently being displayed, sends the output to `gp-special-output-filter'.

Sub-functions are `gp-treat-special-inputp' and `gp-copy-input'
with whom it shares the variables:
`gp-input-start' `gp-input-end' `gp-complete-expression'
`gp-input-start-backetp' `gp-reads-this-buffer'."

  (if (gp-treat-special-inputp)
      nil ;; already treated.
      (gp-copy-input)) ;; does all the work!

  (if gp-complete-expression
  ;; If it is a complete expression do this:
      (progn
        (insert "\n")
        (gp-input-filter)
        (if (> (- gp-input-end gp-input-start) 1023)
  ;;  If large expression, use a temporary file.
          (progn
            (write-region gp-input-start gp-input-end gp-temp-file)
            (process-send-string gp-process (concat "\\r "gp-temp-file"\n")))
  ;;  Else use process-send-region.
          (if gp-input-start-bracketp
              (process-send-region gp-process gp-input-start gp-input-end)
            (process-send-string gp-process "{")
            (process-send-region gp-process gp-input-start gp-input-end)
            ;; a tricky one: if last line had a \\, the final } may not be seen...
            (process-send-string gp-process "\n}"))
          (process-send-string gp-process "\n"))
        (set-marker (process-mark gp-process) (point))
        (if (and gp-should-wait-for-outputp (not gp-trust-mode))
            (progn (gp-wait-for-output gp-input-end)
                   (gp-update-fontification))
            (setq gp-should-wait-for-outputp t))
        (if (and localp (not gp-trust-mode))
            ;; Sometimes the output should not be sent to the output filter:
           (progn
             (save-excursion
               (goto-char gp-input-start)
               (setq localp
                     (not (re-search-forward "\\\\r +" gp-input-end t))))
             (if (and localp gp-separate-window-for-mistakes)
                 (gp-output-filter)))
           (if (and (stringp gp-reads-this-buffer)
                    (buffer-visiblep gp-reads-this-buffer))
               ;; If an error is detected, and a buffer is visible
               ;; containing gp-reads-this-buffer, then we should move the
               ;; point to the place where the error is detected.
               (gp-special-output-filter))))

;; Else (not a complete expression) do this:
      (gp-C-j)
      (message (gp-messager 9))))

(defun gp-send-local-input nil
  "An input is declared to be 'local' if it comes from the *PARI* buffer."
  (interactive) (gp-send-input t))

(defun gp-interrupt ()
  "Interrupts gp.
This is identical to interrupt-shell-subjob in shell-mode."
  (interactive) (interrupt-process nil t))

;;---------------
;; META-COMMANDS
;;---------------

(defun gp-meta-cmd-general (cmd window-option)
  "With 'gp-beginning for window-option, it is 'gp-meta-cmd.
With nil, it is 'gp-quiet-meta-cmd."
  (progn
    (set-buffer "*PARI*")    ;; In case we use it from another buffer,
                             ;; but a gp process is running.
    (goto-char (point-max))
;; Make gp send text to the buffer end, so we can move it to the help buffer.
    (set-marker (process-mark gp-process) (point))
    (let ((temp (point)))
      ;; Send the meta command to gp.
      (process-send-string gp-process (concat cmd "\n"))
      ;; Wait for the gp-prompt to be sent.
      (gp-wait-for-output temp)

      ;; Display the output in the help buffer:
      (let ((copy (buffer-substring-no-properties temp (point-max))))
        (delete-region temp (point-max))
        (if (eq window-option 'gp-beginning)
            ;;Switch to buffer "*gp-help*":
            (gp-window-manager "*gp-help*" window-option)
            (set-buffer (get-buffer-create "*gp-help*"))
            (erase-buffer))
            
        (insert copy)
        (beginning-of-line)  ;; We remove the last prompt line.
        (delete-region (point) (point-max))
        (goto-char (point-min))))))

(defun gp-meta-cmd (cmd)
  "Send cmd to gp, and display output in help buffer"
  (save-excursion
    (let ((wind (selected-window)))
         (gp-meta-cmd-general cmd 'gp-beginning)
         (select-window wind))
         (gp-info-wind-conf)))

(defun gp-quiet-meta-cmd (cmd)
  "Send cmd to gp, and copy output in help buffer without displaying it"
  (save-excursion (gp-meta-cmd-general cmd nil)))

(defun gp-set-prompt (p)
  "Set new gp prompt (and tell both gp and emacs that you have done so)."

  (interactive "sNew prompt: ")
  (let ((my-buffer (buffer-name)) temp)
   (set-buffer "*PARI*")
   (goto-char (setq temp (point-max)))
;; New pattern matches p OR old-pattern
   (gp-make-gp-prompt-pattern (gp-make-prompt-pattern p))
;; Tell gp about the change too!
   (insert (concat "default(prompt,\"" p "\");\n"))
   (process-send-string gp-process (concat "default(prompt,\"" p "\");\n"))
   (set-marker (process-mark gp-process) (point))
   (gp-wait-for-output temp)
   (gp-update-fontification)
;; In case it is called from the menu-bar, do not write anything:
   (message "")
   (set-buffer my-buffer)))

(defun gp-make-prompt-pattern (p)
  "Make the regexp that matches the prompt p."
  ;; We use the buffer *Messages* to analyse the prompt.
  (save-excursion
    (set-buffer "*Messages*")
    (goto-char (point-max))
    (insert "\n" p) ; The "\n" is most probably useless.
    (beginning-of-line)
    (let ((where (point)) a-char)
         (setq p "")
         (while (not (eolp))
          (if (re-search-forward "%[a-zA-Z%]" nil t)
              (setq p
               (concat p
                (regexp-quote (buffer-substring-no-properties where
                                                (match-beginning 0)))
                (progn (setq a-char (buffer-substring-no-properties (1- (point)) (point)))
                       (setq where (point))
                             ;;Options from strftime:
                       (cond ((string= a-char "%") "%")
                             ((member a-char 
                               '("C" "d" "e" "H" "I" "k" "l" "m" "M" "S"
                                 "U" "V" "W" "y"))
                              "[0-9][0-9]")
                             ((member a-char '("D" "T"))
                              "[0-9][0-9]/[0-9][0-9]/[0-9][0-9]")
                             ((string= a-char "R")
                              "[0-9][0-9]:[0-9][0-9]")
                             ((member a-char '("a" "A" "b" "B"))
                              "[A-Z][a-z]*")
                             ((string= a-char "n")
                              "\n")
                             ;; If everything else fails:
                             (t (concat "%" a-char))))))
              ;; No % anymore:
              (goto-char (point-max))
              (setq p
               (concat p
                (regexp-quote (buffer-substring-no-properties where (point-max))))))))
          ;; Now p contains the regexp matching the prompt.
    ;; We erase what we have written on this buffer:
    (beginning-of-line) (backward-char 1)
    (delete-region (point) (point-max))
  ;;Return p:
  p))

(defun gp-set-simple-prompt nil
  "Set the prompt to \"? \"."
  (interactive)
  (gp-set-prompt "? "))

(defun gp-set-time-prompt nil
  "Set a prompt that gives the time."
  (interactive)
  (gp-set-prompt "(%H:%M)> "))

(defun gp-set-date-prompt nil
  "Set a prompt that gives the date."
  (interactive)
  (gp-set-prompt "%d %b %y >> "))

(defun gp-set-separator-prompt nil
  "Set a prompt with a separator "
  (interactive)
  (gp-set-prompt "-------------------------%n(%H:%M)> "))

(defun gp-meta-d ()
  "Send \\d to gp, then display output in the help buffer.
Print the gp defaults."
  (interactive)
  (gp-meta-cmd "\\d"))

(defun gp-meta-t ()
  "Send \\t to gp, then display output in the help buffer.
Print the longword format of PARI types."
  (interactive)
  (gp-meta-cmd "\\t"))

(defun gp-meta-r (file)
  "Send a \\r <file name> command to gp.
Read in gp commands from a file."
  (interactive "fRead from file: ")
  (goto-char (point-max))
  (insert (concat "\\r " (expand-file-name file)))
  (gp-send-input))

(defun gp-meta-w (file num)
  "Send a \\w<num> <file name> command to gp.
Writesgp object %<num> to <file name>."
  (interactive "FWrite to file: \nsObject number %%")
  (goto-char (point-max))
  (insert (concat "\\w"num" " (expand-file-name file)))
  (gp-send-input))

(defun gp-meta-x nil
  "Send \\x to gp, then display output in the help buffer.
Print tree of addresses and contents of last object."
  (interactive)
  (gp-meta-cmd "\\x"))

(defun gp-meta-v nil
  "If gp is running, send \\v to gp, then display output
in the help buffer. Print the version number of this
implementation of pari-gp."
  (interactive)
  (if (processp gp-process) (gp-meta-cmd "\\v")
      (message (gp-messager 11) gp-version)))

(defun gp-meta-s (num)
  "Send \\s or \\s(num) to gp, then display output in the help buffer.
Print the state of the pari stack."
  (interactive "sNumber of longwords (default 0) ")
  (if (equal num "")
    (gp-meta-cmd "\\s")
    (gp-meta-cmd (concat "\\s(" num ")" ))))

(defun gp-meta-a (num)
  "Send \\a or \\a<num> to gp, then display output in the help buffer.
Print object %<num> in raw format."
  (interactive "sPrint object (default last) %%")
  (gp-meta-cmd (concat "\\a" num)))

(defun gp-meta-b (num)
  "Send \\b or \\b<num> to gp, then display output in the help buffer.
Print object %<num> in pretty format."
  (interactive "sPrint object (default last) %%")
  (gp-meta-cmd (concat "\\b" num)))

(defun gp-meta-m (num)
  "Send \\m or \\m<num> to gp, then display output in the help buffer.
Prins object %<num> in prettymatrix format."
  (interactive "sPrint object (default last) %%")
  (gp-meta-cmd (concat "\\m" num)))

(defun gp-meta-q ()
  "Send \\q to gp. Prompt for confirmation before quiting."
  (interactive) 
  (if (y-or-n-p "Quit gp ? ") 
    (progn
     (set-buffer "*PARI*")
     (goto-char (point-max))
     (process-send-string gp-process "\\q\n")
     (setq gp-process nil) ;; Should be automatic with the previous one.
                           ;; Works better like this.
    ))
  (message ""))

(defun gp-break-long-line nil
  "gp will not accept lines longer than 1024.
gp-break-long-line breaks current line 
inserting \\ every (frame-width)-5 chars."
  (interactive)
  (let ((length (min (- (frame-width) 5) 250)))
  (move-to-column length)
  (while (not (looking-at "$"))
    (insert "\\\n")
    (move-to-column length))))

(defun gp-copy-last-input nil
  (interactive)
  (save-excursion
    (goto-char (point-max))
    (if (re-search-backward gp-prompt-pattern nil t 2)
        (progn (goto-char (match-end 0))
               (gp-copy-input)))))

(defun gp-previous-cmd nil
  "Recall previous gp command."
  (interactive)
  (gp-relative-cmd -1))

(defun gp-next-cmd nil
  "Step to gp next command line."
  (interactive)
  (gp-relative-cmd 1))

(defun gp-relative-cmd (dir)
  "Step to previous or next command line according to
the first argument being 1 or -1."
  (while (and (zerop (forward-line dir))
              (not (looking-at gp-prompt-pattern))
              (looking-at "^"))); forward-line at the end of a buffer
  (end-of-line))

(defun gp-toggle-previous-next-behavior nil
  "Change C-p/M-p C-n/M-n from previous-line and next-line to
gp-previous-cmd and gp-next-cmd and reciprocally"
  (interactive)
  (if (equal (key-binding "\C-p") 'previous-line)
      (progn
        (define-key gp-map "\M-p" 'previous-line)
        (define-key gp-map "\M-n" 'next-line)
        (define-key gp-map "\C-p" 'gp-previous-cmd)
        (define-key gp-map "\C-n" 'gp-next-cmd))
    (define-key gp-map "\C-p" 'previous-line)
    (define-key gp-map "\C-n" 'next-line)
    (define-key gp-map "\M-p" 'gp-previous-cmd)
    (define-key gp-map "\M-n" 'gp-next-cmd)))

(defun gp-toggle nil
  "Change some keys. See gp-toggle-previous-next-behavior"
  (interactive)
  (gp-toggle-previous-next-behavior)
  (message (gp-messager 14)))

(defun gp-translate (bool)
  (if bool "On" "Off"))

(defun gp-toggle-locked-mode nil
  "Toggle `gp-locked-modep'."
  (interactive)
  (message (format (gp-messager 84)
             (gp-translate (setq gp-locked-modep (not gp-locked-modep))))
  (gp-define-locked-keys)))

(defun gp-toggle-trust-mode nil
  "Toggle `gp-trust-mode'."
  (interactive)
  (message (format (gp-messager 85)
             (gp-translate (setq gp-trust-mode (not gp-trust-mode))))))

(defun gp-remove-last-output nil
  (interactive)
  (save-excursion
    (goto-char (point-max))
    (when (re-search-backward gp-prompt-pattern nil t)
       (let ((inhibit-read-only t))
         (delete-region gp-input-end (point-max))))))

(defun gp-remove-last-action nil
  (interactive)
  (save-excursion
    (goto-char (point-max))
    (if (re-search-backward gp-prompt-pattern nil t)
        (let ((where (1- (point))))
          (when (re-search-backward gp-prompt-pattern nil t)
            (let ((inhibit-read-only t))
              (delete-region (1- (point)) where)))))))

(defun gp-electric-behavior (choice)
  "Selects RET/M-RET from `sli-electric-terminate-line'
to newline and reciprocally"
  (interactive)
  (setq gp-auto-indent choice)
  (if choice
      (progn
        (define-key gp-script-map "\r"    'sli-electric-terminate-line)
        (define-key gp-script-map "\M-\r" 'newline))
  (define-key gp-script-map "\M-\r" 'sli-electric-terminate-line)
  (define-key gp-script-map "\r"    'newline)))

;;-------------------
;; PART VI : MENU-BAR
;;-------------------

;;---------------
;; MENU BUILDERS
;;---------------

(defconst gp-metakeys-gp-mode-menu
  (list
  (list (gp-messager 48) 
    (vector (gp-messager 49) 'gp-meta-r '(processp gp-process))
    (vector (gp-messager 50) 'gp-meta-w '(processp gp-process))
    "------------------------------------------------"
    (list (gp-messager 51)
      ["Pretty Format"  gp-meta-b (processp gp-process)]
      ["Matrix Pretty Format" gp-meta-m (processp gp-process)]
      ["Raw Format" gp-meta-a (processp gp-process)]
      ["Inner Structure" gp-meta-x (processp gp-process)])
    (list (gp-messager 52)
      (vector (gp-messager 53) 'gp-set-simple-prompt
              ':active t ':included '(processp gp-process) ':key-sequence nil)
      (vector (gp-messager 54) 'gp-set-time-prompt
              ':active t ':included '(processp gp-process) ':key-sequence nil)
      (vector (gp-messager 55) 'gp-set-date-prompt
              ':active t ':included '(processp gp-process) ':key-sequence nil)
      (vector (gp-messager 56) 'gp-set-separator-prompt
              ':active t ':included '(processp gp-process) ':key-sequence nil)
      (vector (gp-messager 42) 'gp-set-prompt '(processp gp-process)))
    "------------------------------------------------"
    ["PARI Types"     gp-meta-t (processp gp-process)]
    ["Default"        gp-meta-d (processp gp-process)]
    ["Version Number" gp-meta-v (processp gp-process)]
    ["Stack Info"     gp-meta-s (processp gp-process)])))

(defconst gp-metakeys-gp-script-mode-menu
  (list
  (list (gp-messager 48)
    ["PARI Types"     gp-meta-t :included (processp gp-process) :active t]
    ["Default"        gp-meta-d :included (processp gp-process) :active t]
    ["Version Number" gp-meta-v t])))

(defun gp-build-main-cmds-menu nil ""
  (nconc
   (if (equal major-mode 'gp-script-mode)
     (list
       (vector (gp-messager 57) 'gp ':active t))
     nil)
   (list (vector (gp-messager 60) 'gp-run-gp t))
   (if (eq major-mode 'gp-script-mode)
       (list (vector (gp-messager 61) 'gp-run-in-region
                                      ':active 'mark-active)) nil)
   (list (vector (gp-messager 62) 'gp-meta-q '(processp gp-process)))))

(defun gp-environment-menu nil
  (list
    (vector (concat (gp-messager 82) " (" (gp-translate gp-locked-modep) ")")
            'gp-toggle-locked-mode ':active t ':key-sequence nil)
    (vector (concat (gp-messager 83) " (" (gp-translate gp-trust-mode) ")")
         'gp-toggle-trust-mode t ':key-sequence nil)
    (vector (gp-messager 81) 'gp-customize-gp-group t ':key-sequence nil)))

(defun gp-build-utilities-menu nil ""
  (nconc
    (list 
          (vector (gp-messager 73) 'gp-skip-to-error t))
    (if (eq major-mode 'gp-mode)
        (list
        (list (gp-messager 74)
         (vector (gp-messager 75) 'gp-copy-last-input '(processp gp-process))
         (vector (gp-messager 76) 'gp-remove-last-output t)
         (vector (gp-messager 77) 'gp-remove-last-action t)))
        nil)
    (list 
      (vector (gp-messager 78) 'gp-toggle t :included '(eq major-mode 'gp-mode))
      (gp-environment-menu))
    ))

;;--------------------------------------
;; MENU-BAR ITEM USED IN GP-SCRIPT-MODE
;;--------------------------------------

(defun gp-init-script-menu-bar nil
   "Add menu-bar item GP if wanted and possible."
   (when (and gp-menu-barp
              (progn (require 'easymenu) (featurep 'easymenu))
              (eq GP-script-menu-map nil))
      (easy-menu-define GP-script-menu-map gp-script-map
       "Menu-bar item used under gp-script-mode."
       (append
         (list "GP")
         (gp-build-main-cmds-menu)                       gp-separator
         gp-metakeys-gp-script-mode-menu                 gp-separator
         (gp-build-utilities-menu)                       gp-separator
         (list (vector (gp-messager 71) 'gp-restore-wind-conf
                                        'gp-registers-list))
         ))
      (add-hook 'menu-bar-update-hook
        '(lambda nil (easy-menu-change '("GP") (gp-messager 79) (gp-environment-menu))))
      (run-hooks 'menu-bar-update-hook)))

;;-------------------------------
;; MENU-BAR ITEM USED IN GP-MODE
;;-------------------------------

(defun gp-init-menu-bar nil
  "Add menu-bar item GP if wanted and possible."
  (when (and gp-menu-barp
             (progn (require 'easymenu) (featurep 'easymenu))
             (eq GP-menu-map nil))
     (easy-menu-define GP-menu-map gp-map
      "Menu-bar item used under gp-mode."
      (append
        (list "GP")
        (gp-build-main-cmds-menu)                       gp-separator
        gp-metakeys-gp-mode-menu                        gp-separator
        (gp-build-utilities-menu)                       gp-separator
        (list (vector (gp-messager 71) 'gp-restore-wind-conf
                                       'gp-registers-list))))
     (add-hook 'menu-bar-update-hook
       '(lambda nil (easy-menu-change '("GP") (gp-messager 79) (gp-environment-menu))))
     (run-hooks 'menu-bar-update-hook)))

;;-----------------------------------------------
;; Customs def that uses the above
;; for initialisation.
;;-----------------------------------------------

(defcustom gp-auto-indent nil
"Non-nil means emacs will try to indent properly each line ended
by a carriage return. Changing its value will exchange the bindings
of \r and \M-\r."
:type 'boolean
:set (lambda (symbol val) (gp-electric-behavior val))
:initialize 'custom-initialize-set ;if you use :set, you should specify :initialize!
:group 'gp-indentation)

;;; pari.el ends here   ----------
