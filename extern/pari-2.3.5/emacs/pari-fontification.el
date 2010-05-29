;; $Id: pari-fontification.el 6643 2005-02-20 18:51:46Z kb $
;; pari-fontification.el --  part of pari.el GP/PARI editing support package.

;; fontification functions
;; See pariemacs.txt  for more details.

;;; TODO:         (please feel free to propose !)
;;  -- Introducing three levels of fontification.
;;  -- Introducing predefined combination of colors.
;;     See mupad-color-scheme and mupad-color-scheme-alist.
;;  -- An initial value of gp-fontifyp to nil ...

(provide 'pari-fontification)

;; Provides: variable:  gp-fontification-keywords
;;           functions: gp-update-fontification, gp-find-global-var, 

;; Of pari.el, it uses:
;;     functions: gp-window-manager, gp-store-wind-conf
;; Of pari-messages.el, it uses: functions: gp-messager
(eval-and-compile
(unless (fboundp 'gp-messager)
  (defun gp-messager (no) (print "Feature pari-messages is absent ..." "*Messages*")))
(unless (fboundp 'gp-window-manager)
  (defun gp-window-manager (a b) (message "Main program pari.el is absent !!")))
(unless (fboundp 'gp-info-wind-conf)
  (defun gp-info-wind-conf nil (message "Main program pari.el is absent !!")))
(unless (fboundp 'gp-store-wind-conf)
  (defun gp-store-wind-conf nil (message "Main program pari.el is absent !!"))))

(require 'font-lock)
(eval-when-compile
  (fset 'x-defined-colors nil)
  ;; for development:
  ;;(setq byte-compile-warnings (list 'free-args 'unresolved 'callargs 'redefine 'obsolete))
  ;; for users:
  (setq byte-compile-warnings (list 'unresolved 'redefine 'obsolete))
  )

(defun gp-update-fontification-buffers nil
"Update (/un)-fonctification on all the buffers
that are in gp-mode or in gp-script-mode."
  (interactive)
  (save-excursion
      (mapcar
      (lambda (abuffer)
              (set-buffer abuffer)
              (if (memq major-mode '(gp-script-mode gp-mode))
                  (if gp-fontifyp (font-lock-fontify-buffer)
                     (font-lock-unfontify-buffer))))
      (buffer-list))
  (message "")))

(defun gp-set-fontifyp (sym val)
  (set sym (and (eq window-system 'x) (x-display-color-p) val))
  (gp-update-fontification-buffers))

(defcustom gp-fontifyp t
"If this variable is nil neither fontify GP scripts nor *PARI* buffer.
Partly modified internally"
:type 'boolean
:set 'gp-set-fontifyp
:initialize 'custom-initialize-reset
:group 'gp-font-lock-and-completion)

(eval-and-compile
(mapcar (lambda (agpplace) (eval (list 'defvar (eval agpplace) agpplace)))
     (list
     ''gp-error ''gp-history ''gp-prompt  ''gp-output ''gp-input
     ''gp-help  ''gp-timer   ''gp-comment ''gp-string ''gp-control-statement
     ''gp-default-keywords   ''gp-default-set         ''gp-input-cmd
     ''gp-global-var         ''gp-function-proto      ''gp-function-args)))

(defvar pari-fontification-keywords
  (list
    '("\\<\\(buffersize\\|co\\(lors\\|mpatible\\)\\|debug\\(mem\\)?\\|echo\\|format\\|h\\(elp\\|size\\)\\|logfile\\|output\\|p\\(a\\(risize\\|th\\)\\|r\\(imelimit\\|ompt\\)\\|sfile\\)\\|\\(real\\|series\\)precision\\|simplify\\|strictmatch\\|timer\\)\\>"  (1 gp-default-keywords))
    '("\\<\\(return\\|break\\|next\\|if\\|until\\|while\\|sum\\|for\\(div\\|prime\\|step\\|vec\\|subgroup\\)?\\)\\>" (1 gp-control-statement))
    '("\\<\\(default\\)(" (1 gp-default-set))
    '("^ *\\\\[a-z].*$" . gp-default-set)
    ;; In the two following ones, in case we meet a list, the separators are also painted...
    '("\\<\\([a-zA-Z]\\w*\\)(\\([^)]*\\)) *=[^=]" (1 gp-function-proto) (2 gp-function-args))
   '("\\<global[ \t]*(\\([^)]*\\))" (1 gp-global-var t))
    )
  "Common keywords to be fontified in gp- and gp-script- mode.")

(defvar gp-fontification-keywords
  (append
  (list
    ;; Careful ! everything involving gp-prompt-pattern
    ;; should be redefined in gp-make-gp-prompt-pattern
    (list gp-prompt-pattern (list 1 'gp-prompt))
    '("^ *%[0-9]* =" .  gp-history)
    '(gp-match-output (0 gp-output))
    '("\\*\\*\\*.*$" .   gp-error)
    '("^[a-zA-Z][a-zA-Z0-9_]*([^ )]*): \\(*\\)$" (1 gp-help))
    '(gp-match-input (0 gp-input))
    '("time = \\([0-9][hmn,0-9 ]* ms\.\\)"  (1 gp-timer)))
   pari-fontification-keywords)
  "Patterns to be fontified under gp-mode.")

(defvar gp-script-fontification-keywords
  (append pari-fontification-keywords
  (list
    '(gp-find-global-var (1 gp-global-var))
    '("\\<read\\|install\\>" . gp-input-cmd)
    ))
  "Patterns to be fontified under gp-script-mode.")

(defun gp-fontification-switch nil
   (interactive)
   (gp-set-fontifyp 'gp-fontifyp (not gp-fontifyp)))

;;------------------------
;; PART V : HIGHLIGHTING
;;------------------------

(defsubst gp-default-face (gp-face default-face doc)
  (if (not (null (get gp-face 'saved-face)))
      (custom-declare-face gp-face (get 'saved-face default-face) doc)
    (copy-face default-face gp-face)
    (set-face-documentation gp-face doc)))

(defun gp-init-font-lock-faces nil
  "Define gp-faces."
  (interactive)
  (gp-default-face gp-error font-lock-warning-face
                   "*Face used in GP to highlight errors.")
  (gp-default-face gp-comment font-lock-comment-face
                   "*Face used in GP to highlight comments. Default is font-lock-comment-face.")
  (gp-default-face gp-string font-lock-string-face
                   "*Face used in GP to highlight string. Default is font-lock-string-face.")
  (gp-default-face gp-function-proto font-lock-function-name-face
                   "*Face used in GP to highlight function names in definitions.
Default is font-lock-function-name-face.")
  (gp-default-face gp-function-args font-lock-variable-name-face
                   "*Face used in GP to highlight function arguments in definitions.
Default is font-lock-varaible-name-face.")
  (gp-default-face gp-global-var font-lock-constant-face
                   "*Face used in GP Script to highlight global variables.
Default is font-lock-constant-face.")
  (gp-default-face gp-history font-lock-constant-face
                   "*Face used in GP to highlight function arguments in definitions.
Default is font-lock-constant-face.")
  (gp-default-face gp-default-keywords font-lock-builtin-face
       "*Face used in GP to highlight some keywords (buffersize, simplify ...).
Default is font-lock-builtin-face.")
  (gp-default-face gp-control-statement font-lock-keyword-face
       "*Face used in GP to highlight control-statements (for, while ...).
Default is font-lock-keywords.")
  (defface gp-prompt
    '((((class grayscale) (background light)) (:foreground "LightGray" :bold t))
      (((class grayscale) (background dark)) (:foreground "DimGray" :bold t))
      (((class color) (background light)) (:foreground "Orchid"))
      (((class color) (background dark)) (:foreground "LightSteelBlue"))
      (t (:bold t)))
    "*Face used in GP to highlight prompt.")
  (defface gp-output
    '((((class grayscale) (background light)) (:foreground "DimGray" :italic t))
      (((class grayscale) (background dark)) (:foreground "LightGray" :italic t))
      (((class color) (background light)) (:foreground "RosyBrown"))
      (((class color) (background dark)) (:foreground "LightSalmon"))
      (t (:italic t)))
    "*Face used in GP to highlight outputs.")
  (defface gp-input
    '((((class grayscale) (background light)) (:foreground "LightGray" :bold t))
      (((class grayscale) (background dark)) (:foreground "DimGray" :bold t))
      (((class color) (background light)) (:foreground "Purple"))
      (((class color) (background dark)) (:foreground "Cyan"))
      (t (:bold t)))
    "*Face used in GP to highlight inputs.")
  (defface gp-help
    '((((class grayscale) (background light)) (:foreground "LightGray" :bold t))
      (((class grayscale) (background dark)) (:foreground "DimGray" :bold t))
      (((class color) (background light)) (:foreground "Orchid"))
      (((class color) (background dark)) (:foreground "LightSteelBlue"))
      (t (:bold t)))
    "*Face used in GP to highlight help messages in *PARI* buffer.")
  (defface gp-timer
    '((((class grayscale) (background light)) (:foreground "DimGray" :italic t))
      (((class grayscale) (background dark)) (:foreground "LightGray" :italic t))
      (((class color) (background light)) (:foreground "RosyBrown"))
      (((class color) (background dark)) (:foreground "LightSalmon"))
      (t (:italic t)))
    "*Face used in GP to highlight time.")
  (defface gp-default-set
    '((((class grayscale) (background light)) (:foreground "LightGray" :bold t))
      (((class grayscale) (background dark)) (:foreground "DimGray" :bold t))
      (((class color) (background light)) (:foreground "Purple"))
      (((class color) (background dark)) (:foreground "Cyan"))
      (t (:bold t)))
     "*Face used in GP to highlight default-set.")
  (defface gp-input-cmd
    '((((class grayscale) (background light)) (:foreground "LightGray" :bold t))
      (((class grayscale) (background dark)) (:foreground "DimGray" :bold t))
      (((class color) (background light)) (:foreground "Orchid"))
      (((class color) (background dark)) (:foreground "LightSteelBlue"))
      (t (:bold t)))
    "*Face used in GP Script to highlight input-cmd (read, install)."))

(defun gp-match-output (limit)
  "Set match-data 0 to limits of next output."
  (if (re-search-forward "^ *%[0-9]* = +" limit t)
      (let ((beg (point)))
        (if (re-search-forward gp-prompt-pattern limit t)
            (progn (goto-char (match-end 0))
                   (set-match-data (list beg (1- (match-beginning 0))))
                   t)
          nil))
    nil))

(defun gp-find-global-var (limit) 
 "A parser to find global variables. Called on a gp-program outside
a function-definition, gives position via (cons start end) of
next global-variable-definition not surrounded by {} and set the
point at the end of the line. Answer nil if no global-variable is found.
The end delimiter of a function definition surrounded by {} is
'}\n' and the same holds with function definitions of the style
'fun(var)={foo}'.
LIMIT is not used."
  (let ((answer nil))
    (while (looking-at (concat comment-start-skip
                               "\\|[ \\|\t\\|\n]+\\|{\\([^}]\\|}[^\n]\\)*}\n\\|\\<[a-zA-Z]\\w*([^)]*) *={\\([^}]\\|}[^\n]\\)*}\n\\|\\<[a-zA-Z]\\w*([^)]*) *=\\([^=\\\\\"]\\|\\\\[ \t]*\\(\\\\\\\\.*$\\|/\\*\\([^\\*]\\|\\*[^/]\\)*\\*/\\)?\n\\|\"\\([^\"]*\\|\\\\\"\\)*\"\\)\\([^\\\\\n\"]\\|\"\\([^\"]*\\|\\\\\"\\)*\"\\|\\\\[ \t]*\\(\\\\\\\\.*$\\|/\\*\\([^\\*]\\|\\*[^/]\\)*\\*/\\)?\n\\)*\n\\|\\<[a-zA-Z]\\w*([^)]*)[;\n]"))
  ;; We look at a single line comment, or a long comment,
  ;; or a space/tab/newline character, or a function definition between {},
  ;; or a function definition of the type fun(var)={foo},
  ;; or a function definition not between {}, or a function call,
  ;; or any line without an equality sign.
  ;; And skip them.
      (goto-char (match-end 0)));(prin1 (point))(prin1 " "));(print "GPGLOBAL") ;(sit-for 10)
  ;; We look whether there is a global-variable being defined here:
    (if (looking-at "^\\<\\([a-zA-Z]\\w*\\)=[^=].*$")
        (progn
          (setq answer (cons (match-beginning 1) (match-end 1)))
          (goto-char (match-end 0))))
    (unless answer (goto-char limit))
    answer))

(defun gp-find-global-var (limit) (goto-char limit) nil)

(defsubst gp-search-forward-string-delimiter (lim)
  "Give the position of next \" preceded by an even number
of \\ . Move point after this point. Nil if no such place before lim."
  ;; Inspired from fontify-string-find in hilit19.el.
  (let (p)
    (while (and (setq p (search-forward "\"" lim t))
                (save-excursion
                  (forward-char -1)
                  (not (zerop (% (skip-chars-backward "\\\\") 2))))))
    p))

(defun gp-update-fontification nil "Update fontification."
  (interactive)
  (when gp-fontifyp
    (font-lock-fontify-buffer)))

(defun gp-turn-on-lazy-font-lock nil ""
  (interactive)
  (require 'lazy-lock)
  (when (featurep 'lazy-lock) (lazy-lock-mode)))

(defun gp-customize-faces nil
  (interactive)
  (if gp-tutorial-requiredp
      (let ((wind (selected-window)) s
            (msg (gp-messager 28))
            (colors-list (if (fboundp 'x-defined-colors)
                             (sort (x-defined-colors) 'string-lessp)
                             '("No colours found !!!"))))
           (gp-window-manager "*gp-help*" 'gp-beginning-temp)
           (insert msg)
           (fill-region (point-min) (point-max) 'left)
           ;; Following taken from list-colors-display of facemenu.el:
           (while colors-list
             (setq s (point))
             (insert (car colors-list))
             (indent-to 20)
             (put-text-property s (point) 'face 
                                (cons 'background-color (car colors-list)))
             (setq s (point))
             (insert "  " (car colors-list) "\n")
             (put-text-property s (point) 'face 
                                (cons 'foreground-color (car colors-list)))
             (setq colors-list (cdr colors-list)))
           (goto-char (point-min))
           (gp-info-wind-conf)
           (select-window wind)))
  (customize-apropos-faces "font-lock-.*\\|gp-.*"))

(defun gp-customize-gp-group nil
  (interactive)
  (gp-store-wind-conf)
  (customize-group "gp")
  (message (gp-messager 4)))

(defun gp-color-menu nil
  "Build the Colors menu"
  (when (and (eq window-system 'x) (x-display-color-p))
     (append
      (if (eq major-mode 'gp-script-mode)
          (list (vector (gp-messager 45) 'gp-turn-on-lazy-font-lock
                        ':active t ':key-sequence nil 
                        ':included '(and gp-fontifyp (eq major-mode gp-script-mode)))) nil)
      (list (vector (gp-messager 44) 'gp-update-fontification
                    ':active t ':included: 'gp-fontifyp)
            (vector (gp-messager 46) 'gp-fontification-switch
                    ':active t ':key-sequence nil ':included: 'gp-fontifyp)
            (vector (gp-messager 47) 'gp-update-fontification-buffers
                    ':active t ':key-sequence nil)
            (vector (gp-messager 79) 'gp-customize-faces ':active t
                    ':key-sequence nil ':included: 'gp-fontifyp)))))

(add-hook 'menu-bar-update-hook
  '(lambda nil
     (when (memq major-mode '(gp-mode gp-script-mode))
       (easy-menu-change '("GP") (gp-messager 43) (gp-color-menu) (gp-messager 72)))))

(add-hook 'pari-mode-hook
  '(lambda nil
     (define-key gp-map        "\C-l" (function gp-update-fontification))
     (define-key gp-script-map "\C-l" (function gp-update-fontification))
     (require 'font-lock)
     (gp-init-font-lock-faces)))

(add-hook 'gp-script-mode-hook
  '(lambda nil
     (make-local-variable 'font-lock-defaults)
     (make-local-variable 'font-lock-comment-face)
     (setq font-lock-comment-face (eval gp-comment))
     (make-local-variable 'font-lock-string-face)
     (setq font-lock-string-face (eval gp-string))
     (setq font-lock-defaults '(gp-script-fontification-keywords nil nil nil))))

(add-hook 'gp-mode-hook
  '(lambda nil
     (make-local-variable 'font-lock-defaults)
     (make-local-variable 'font-lock-comment-face)
     (setq font-lock-comment-face (eval gp-comment))
     (make-local-variable 'font-lock-string-face)
     (setq font-lock-string-face (eval gp-string))
     (setq font-lock-defaults '(gp-fontification-keywords nil nil nil))))

;; pari-fontification.el ends here.

