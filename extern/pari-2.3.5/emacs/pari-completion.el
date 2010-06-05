;; $Id: pari-completion.el 4609 2003-04-18 09:20:01Z ramare $
;; pari-completion.el --  part of pari.el GP/PARI editing support package.

;; documentation functions
;; See pariemacs.txt  for more details.
 
(provide 'pari-completion)

;; pari.el will use the function 'gp-quit-cpl-edit.
;; Also extends pari-mode-hook to set some key-bindings.
;; pari-help.el will use 'gp-c-array.

;; Of pari.el, it uses:
;; variables:
;;     gp-process, gp-prompt-pattern, gp-tutorial-requiredp
(defvar gp-tutorial-requiredp)
(defvar gp-prompt-pattern)
;; functions: 
;;     gp-restore-wind-conf, gp-background,
;;     gp-store-wind-conf, gp-backward-wind-conf,
;;     gp-window-manager
;; Of pari-fontification.el, it uses gp-find-global-var.
(require 'pari-messages)
;; Of pari-messages.el, it uses gp-messager.
(eval-and-compile
(unless (fboundp 'gp-messager)
  (defun gp-messager (no) (print "Feature pari-messages is absent ..." "*Messages*")))
(unless (fboundp 'gp-window-manager)
  (defun gp-window-manager (a b) (message "Main program pari.el is absent !!")))
(unless (fboundp 'gp-info-wind-conf)
  (defun gp-info-wind-conf nil (message "Main program pari.el is absent !!")))
(unless (fboundp 'gp-store-wind-conf)
  (defun gp-store-wind-conf nil (message "Main program pari.el is absent !!")))
(unless (fboundp 'gp-backward-wind-conf)
  (defun gp-backward-wind-conf nil (message "Main program pari.el is absent !!")))
(unless (fboundp 'gp-restore-wind-conf)
  (defun gp-restore-wind-conf nil (message "Main program pari.el is absent !!")))
(unless (fboundp 'gp-background)
  (defun gp-background nil (message "Main program pari.el is absent !!")))
(unless (fboundp 'gp-find-global-var)
  (defun gp-find-global-var (l) (print "Feature pari-fontification is absent ..." "*Messages*"))))

(defvar gp-process nil "Defined in pari.el")
(defcustom gp-readline-enabledp t
"t if readline is enabled. Emacs will try to set it properly
whenever a gp-session is started."
:type 'boolean    :group 'gp-font-lock-and-completion)

(defcustom gp-additional-cpl-file ""
"Name (string) of a completion file used in supplement for completion.
This file should have the format of 'gp-menu files."
:type 'file       :group 'gp-font-lock-and-completion)

(defconst gp-c-array (make-vector 800 0)
  "obarray used for completing gp command names.")
;; pari-2.0.19-beta contains 514 function names.
;; We extend it by 286 for local ones.

(defvar gp-cpl-lists-list
    '(gp-c-array)
 "List of the lists/arrays to be used for completion on top of the
completion already delivered by readline if present and by the general
'gp-c-array which has to be the first element of this list.")

(defconst gp-function-proto-pstart "\\(^ *\\|{\\)\\([a-zA-Z][_a-zA-Z0-9]*\\)([^)]*) *=[^=]"
"Regexp matching the beginning of a function-definition")

(defvar gp-completion-input-start 0)

(defun gp-choose-complete nil
  "Try to see whether readline is enabled or not
and select proper completion function. To be used
when the buffer *PARI* is selected."
  (save-excursion
    (goto-char (point-min))
    (if (re-search-forward "readline .*\\(dis\\|en\\)abled" nil t)
      (progn
        (forward-char -6)
        (setq gp-readline-enabledp (looking-at "n")))
    ;; Else use default:
      (message (gp-messager 1)))))

(defun gp-proper-name (filename)
  "We replace the dots in filename by -."
  (mapconcat
     (lambda (achar) (char-to-string (if (= achar ?.) ?- achar)))
     (string-to-list filename) ""))

(defun name-extension (filename)
  "Return the extension suffix of filename, if any"
  (if (> (length filename) (length (file-name-sans-extension filename)))
      (substring filename (1+ (length (file-name-sans-extension filename))))
      ""))

;;-------------------------
;; GP COMPLETION FUNCTIONS
;;-------------------------

(defun gp-mouse-2 (event)
  "A kind of hook for 'mouse-choose-completion."
  (interactive "e")
  (funcall 'mouse-choose-completion event)
  ;; 'mouse-choose-completion comes from the standard file "mouse.el".
  (gp-restore-wind-conf) (forward-word 1))

(defun gp-clear-list (lst)
  "Remove the lists `(\"\")' from LST."
  (let ((newlist nil))
    (mapcar (lambda (liststring)
              (or (string= (car liststring) "")
                  (setq newlist (cons liststring newlist))))
            lst)
    newlist))

(defun gp-clear-list2 (lst)
  "Remove the empty words from LST."
  (let ((newlist nil))
    (mapcar (lambda (astring)
              (or (string= astring "")
                  (setq newlist (cons astring newlist))))
            lst)
    newlist))

(defun gp-make-cpl-list (abuffer)
"Take a buffer in the format of pari.menu, and create the list
whose name is the concatenation of \"gp-cpl-\" and the buffer-name
and which contains all the non-commented lines of the buffer.
The file must have at least one comment line, starting with #, All
lines before the first comment line are IGNORED. Finally add this list
name to 'gp-cpl-lists-list."
  (save-excursion
   (let ((lst nil) lst-aux astring)
     (set-buffer abuffer)
     (save-restriction
       (widen)
       (goto-char (point-min))
       (re-search-forward "#")
       (while (not (eobp))
         (forward-line 1)
         (or (looking-at "#")
             (add-to-list 'lst
                (list
                 (buffer-substring-no-properties (point)
                                   (line-end-position))))))
       (setq astring
             (concat "gp-cpl-"
                     (gp-proper-name (buffer-name))))
       (make-symbol astring)
       (set (intern astring) (gp-clear-list lst))
       (setq lst-aux (list (intern astring)))
       (setcdr lst-aux (cdr gp-cpl-lists-list))
       (setcdr gp-cpl-lists-list lst-aux)
       (kill-buffer abuffer)
      ))))

(defun gp-cpl-file (afile)
  "Same as `gp-make-cpl-list' except that we start with a file."
  (interactive "fFile of command names: ")
    (gp-make-cpl-list (find-file-noselect afile)))

(defun gp-add-symbol  (name)
  "Add a name to the obarray, if it is not already there."
  (make-symbol name)
  (intern name gp-c-array))

(defun gp-find-word-to-complete nil
  (save-excursion
   (let ((pt (point)) ans)
      (if (char-equal (preceding-char) ?() (forward-char -1))
      (if (not (bolp))
          (progn
            (forward-char -1)           
            (if (looking-at "\\w")
                (progn (forward-char 1) (forward-word -1))
              (forward-char 1))))
      ;; In case it is a command-word:
      (if (= (preceding-char) ?\\) (forward-char -1))
      (setq ans (if (= (point) pt) " "
                  (buffer-substring-no-properties (point) pt)))
      ;(princ "\n") (princ (list "(gp-find-word-to-complete) word:" ans))
      ans)))

(defun gp-string-to-list (astring)
  "ASTRING is a succession of gp-words separated by spaces or newlines.
The answer is the list of these words."
  (let ((lst nil) (beg 0) (end 1))
    (while (<= end (length astring))
      (cond ((member (aref astring (1- end)) '(?\  ?\n))
             (if (not (= beg (1- end)))
                 (setq lst (nconc lst
                                  (list (substring astring beg (1- end))))))
             (setq beg end end (1+ end)))
            (t (setq end (1+ end)))))
    ;; taking care of the last one:
    (if (not (= beg (1- end)))
        (setq lst (nconc lst (list (substring astring beg (1- end))))))
    lst))

(defun gp-sort-and-minimise (list1 list2)
  "Take two lists of strings and build the list of all their
elements with no repetition and sorted out."
   (let ((lst (sort (nconc list1
                           (mapcar
                             (lambda (elt) (if (member elt list1) "" elt))
                             list2))
                    'string-lessp)))
    (if (string= (car lst) "") (cdr lst) lst)))


(defun gp-make-standard-word (word)
  "If WORD has a final \"()\", remove it."
 ;; When asking for completion and there is a unique completion, readline
 ;; adds sometimes `()' at the end of the completion.
  (if (and (> (length word) 1)
           (string= (substring word (- (length word) 2)) "()"))
      (substring word 0 (- (length word) 2))
      word))

(defsubst standard-string= (word1 word2)
  (string= (gp-make-standard-word word1)
           (gp-make-standard-word word2)))

(defsubst gp-standard-lst (word comp)
   (cond ((and (string= (car comp) "") (null (nth 1 comp)))
          (list ""))
         ((null (nth 1 comp))
          (list (concat word (car comp))))
         (t (nth 1 comp))))

(defun gp-merge-cpls (word comp1 comp2)
  (let* ((lst1 (gp-standard-lst word comp1))
         (lst2 (gp-standard-lst word comp2))
         (a-local-cpl-list (mapcar 'list (gp-sort-and-minimise lst1 lst2))))
     (gp-ask-cpl-via-list word 'a-local-cpl-list)))

(defun gp-ask-cpl-via-list (word lst)
  "Careful! LST is a symbol whose value is a list of completion type,
ie a list of lists whose cars are strings used for completion."
  ;; LST can be an array also.
  (setq lst (symbol-value lst))
  (let ((comp (try-completion word lst))
        to-insert fun-list)
    (cond ((equal comp nil)         ; No completion.
           (list "" nil))
          ((equal comp t)           ; Already complete.
           (list "" nil))
          ((> (length comp) (length word)) ; Some completion with a kernel.
           (setq to-insert (substring comp (length word)))
           (setq fun-list
                 (all-completions comp lst))
           (if (< (length fun-list) 2)
               (list to-insert nil)  ; Unique completion.
               (list to-insert fun-list)))
          (t (setq fun-list 
                   (all-completions comp lst))
             (if (< (length fun-list) 2)
                 (list "" nil)       ; Unique completion.
                 (list "" fun-list))))))

(defun gp-ask-cpl-via-readline (context)
  (let ((to-insert nil) (fun-list ""))
                    
   (if (gp-background)
    (save-excursion
      (set-buffer "*PARI*")
      (goto-char (point-max))
      (set-marker (process-mark gp-process) (point))
      (let ((temp (point)) (last nil))

  ;; ask for all completions (readline command)
        (process-send-string gp-process (concat context "\t" ))
        (let ((notdone t))
	  (while notdone 
	    (accept-process-output gp-process);(sleep-for 1)(princ " .. ")
	    (let ((p (point)))
	      (if (or
		    (not (and (processp gp-process) 
			    (eq 'run (process-status gp-process))))
		  (search-backward "@E_N_D" (1+ temp) t))
	;; If gp is not running, or @E_N_D  has appeared, stop.
	      (progn 
		(message (gp-messager 6))
		(setq notdone nil last (point)))
	;; Else wait a bit longer.
	      (message (gp-messager 15)) (goto-char p)))))

      ;; Get end of completed-part:
      (search-backward "@" nil t)
      (setq to-insert (buffer-substring-no-properties temp (point)))
      (forward-char 1) ;; In order to skip the "@".
      ;; Possible further completions:
      (if (< (point) last)
	(setq fun-list (buffer-substring-no-properties (point) (1- last))))
      (delete-region temp (point-max))
      ;; clear line in the gp-process:
      (process-send-string gp-process "\C-A\C-K"))))

   (list to-insert (gp-string-to-list fun-list))))

(defun gp-general-complete (completion-function word)
  "Answer a list whose car is an extension of WORD, and whose cdr
is a list of list of possible matching words."
    (let ((ans (funcall completion-function word)))
    ;; 'gp-find-word-to-complete puts the point at
    ;; the end of the word to complete.

    ;; Insert the beginning of the completion
    ;; BEFORE any window change :    
    (if (not (string= (car ans) ""))
        (progn
          (insert (car ans))
          ;; In case of a direct completion via readline:
          (if (char-equal (preceding-char) ?)) (forward-char -1))))

    (if (equal (nth 1 ans) nil)
    ;; at most one match:
	(if (and (get-buffer "*Completions*")
                 (get-buffer-window "*Completions*"))
            ;; Occurs whenever an earlier completion has
            ;; been asked for.
            (progn
              (gp-restore-wind-conf)
              (forward-word 1)
              ;; In case of a completion via readline:

              (if (and (char-after (point))
                       (char-equal (char-after (point)) ?()) (forward-char 1))
              (if (char-equal (preceding-char) ?)) (forward-char -1))))
    ;; more than two matches:
    (if (string= (car ans) "")
      ;; We do not display anything if a partial completion was possible:
      (progn
        (if (not (and (get-buffer "*Completions*")
                      (get-buffer-window "*Completions*")))
            ;; No use storing wind-conf if some completion is in
            ;; progress.
            (gp-store-wind-conf))
        (with-output-to-temp-buffer "*Completions*"
	  (display-completion-list (nth 1 ans))))))))

(defun gp-ask-cpl-via-readline-and-emacs (word)
  (interactive)
  (let ((lst
         (if (or (and gp-readline-enabledp (string= word ""))
                 (and gp-readline-enabledp gp-process
                      (equal (process-buffer gp-process) (current-buffer))))
             ;; Do not use general completion (let readline work !):
             '("" nil)
             ;; Ask for general completion:
           (gp-ask-cpl-via-list word (car gp-cpl-lists-list)))))
    (mapcar
      (lambda (a-cpl-list)
        (setq lst
          (gp-merge-cpls
            word
            (gp-ask-cpl-via-list word a-cpl-list)
            lst)))
      (cdr gp-cpl-lists-list))

    (cond ((and gp-readline-enabledp gp-process
                (equal (process-buffer gp-process) (current-buffer)))
           (save-excursion
              (if (re-search-backward gp-prompt-pattern nil t)
                  (setq gp-completion-input-start (match-end 0))  ;; end of prompt
                (setq gp-completion-input-start (point-min))))
           ;(princ "\n")
           ;(princ "((gp-ask-cpl-via-readline-and-emacs) readline enabled and inside gp buffer)")
           (gp-merge-cpls
              word lst
              (gp-ask-cpl-via-readline
                (buffer-substring-no-properties gp-completion-input-start (point)))))
          (gp-readline-enabledp
           ;(princ "\n")
           ;(princ "((gp-ask-cpl-via-readline-and-emacs) readline enabled and outside gp buffer)")
           (gp-merge-cpls
              word lst
              (gp-ask-cpl-via-readline
                (buffer-substring-no-properties (line-beginning-position) (point)))))
          (t lst))))

(defun gp-complete nil
  (interactive)
  (gp-general-complete 'gp-ask-cpl-via-readline-and-emacs
                       (gp-find-word-to-complete)))

;;------------------
;; COMPLETION FILES
;;------------------

(defsubst gp-cpl-stamp (my-cpl-file)
  "Put a completion-file-stamp on a buffer."
  ;; Do not convert that in any other langage ! See gp-actualise-stamp.
  (insert (format "\nCompletion File Name: %s\n\n" my-cpl-file))
  (insert
    "----------------------------------------------------------------\n"
    "                     Created: " (current-time-string) "\n"
    "                     By:      " (user-full-name) "\n\n"
    "                     Last Modification: " (current-time-string) "\n"
    "----------------------------------------------------------------\n"
    "\n### Function Names : (one per line)\n"))

(defsubst gp-actualise-stamp nil
  "Actualise the completion-file-stamp of a buffer."
   (goto-char (point-min))
   (if (re-search-forward "Last Modification: " nil t)
   ;; We have found this string and update what's behind:
       (let ((kill-whole-line nil)) ;; local value of global parameter.
         (backward-char 1) ;;so that we are sure something is on this line. 
         (kill-line)
         (insert " " (current-time-string)))))

;; Edition of completion file. We follow a loose way of working
;; in case the user edits other buffers in between.

(defun gp-edit-cpl-file (my-cpl-file)
  "Edit my-cpl-file."
  (interactive
    (list (gp-read-input (gp-messager 33)
                   (concat (gp-possible-file-name) ".cpl") "" t)))
 
    (gp-store-wind-conf)
    (or (file-exists-p (expand-file-name my-cpl-file))
        ;; If the file does not exist, create it (the list may exists though):
        (gp-prepare-cpl-file t))
    (switch-to-buffer-other-window
       (find-file-noselect my-cpl-file))
    (goto-char (point-min))
    (if (eobp) (gp-cpl-stamp my-cpl-file)
               (re-search-forward "#.*$" nil t)
               (goto-char (match-end 0))
               (if (eobp) (insert "\n") (forward-char 1)))
    (message (gp-messager 16)))

(defsubst gp-cpl-bufferp (abuffer)
  (string= (name-extension abuffer) "cpl"))

(defun gp-quit-cpl-edit nil
  (interactive)
  (if (gp-cpl-bufferp (buffer-name))
      (progn
        ;; After entering 'gp-edit-cpl-file,
        ;; the user may have edited another completion file...
        ;; We don't bother since nothing bad will happen. The
        ;; behaviour of emacs may simply daze the user.
        (gp-actualise-stamp)
        (save-buffer 0)  ;; No backup file.
        (gp-backward-wind-conf)
       )))

(defsubst gp-make-cpl-help (file)
  (if gp-tutorial-requiredp
      (let ((wind (selected-window)))
       (gp-window-manager "*gp-help*" 'gp-beginning-temp)
       (insert (format (gp-messager 30) file file file file))
       (fill-region (point-min) (point-max) 'left)
       (select-window wind))))

(defun gp-show-help (astring)
   (gp-window-manager "*gp-help*" 'gp-beginning-temp)
   (insert astring)
   (setq fill-column (1- (window-width (get-buffer-window "*gp-help*"))))
   (fill-individual-paragraphs (point-min) (point-max) 'left)
   ;; Remove help window :
   (gp-window-manager "*gp-help*" 'gp-remove-help-old-config)
   (gp-restore-wind-conf))

(defun gp-cpl-file-info nil
  (interactive)
  (gp-show-help (gp-messager 29)))

(defmacro gp-cpl-file-has (astring)
  "t if the edited completion file has the string ASTRING
at a beginning of line followed by \n or a  space or a #.
Also t when ASTRING is the empty string."
 (`
  (if (string= (, astring) "") t
    (save-excursion
     (goto-char (point-min))
     (re-search-forward "#")  ; There exists such a line.
     (end-of-line)
     (if (eobp)
         nil  ;; Return value is nil.
        (forward-line 1)
        (re-search-forward
         (concat "^" (regexp-quote (, astring)) "[\n| |#]") nil t))))))

(defun gp-prepare-cpl-file (option)
"  Write in the file 'buffername.cpl' which has the format of a completion
file (i.e. a gp-menu file) the names of the functions and of the  global
variables of the visited file. OPTION is t means save the buffer on file,
nil means don't do that if the file wasn't existing already."
  (let* ((file (buffer-name)) (my-cpl-file (concat file ".cpl")))
   (or option (gp-make-cpl-help file))
   ;; Prepare buffer:
   (save-excursion
    (if (or option
            (file-exists-p (expand-file-name my-cpl-file)))
      (set-buffer (find-file-noselect my-cpl-file))
      (set-buffer (get-buffer-create my-cpl-file)))
    (if (file-exists-p (expand-file-name my-cpl-file))
        (progn
           ;; Assume it has the format of a completion-file:
           (re-search-forward "#" nil t)
           (end-of-line)
           (if (eobp) (insert "\n")
               (forward-line 1) (beginning-of-line)
               (kill-region (point) (point-max))))
        (if option (gp-cpl-stamp my-cpl-file)
            (insert "\n### Function Names : (one per line)\n")))
    ;; Add function names:
    (set-buffer file)
    (goto-char (point-min))
    (let ((thelist nil))
    (while (re-search-forward gp-function-proto-pstart nil t)
           (add-to-list 'thelist (match-string 2)))
    (setq thelist (sort thelist (function string-lessp))) ; We order things.
    (set-buffer my-cpl-file)
    (mapcar (lambda (fn) (insert fn "\n")) (gp-clear-list2 thelist)))

    ;; Prepare buffer for names of global variables:
    (insert (gp-messager 25) "\n")

    ;; Add global-variable-names:
    (when (fboundp 'gp-find-global-var)
      (set-buffer file)
      (goto-char (point-min))
      (let (theplace (thelist nil))
        (while (setq theplace (gp-find-global-var nil))
          (add-to-list 'thelist
                       (buffer-substring-no-properties (car theplace) (cdr theplace))))           
        (setq thelist (sort thelist (function string-lessp))) ; We order things.
        (set-buffer my-cpl-file)
        (mapcar (lambda (fn) (insert fn "\n")) (gp-clear-list2 thelist))))
    
    (if (or option (file-exists-p my-cpl-file))
        (progn
          ;; Prepare buffer for closing, no backup-file:
          (gp-actualise-stamp)
          (save-buffer 0)))
    ;; Add it to the possible completions:
    (gp-make-cpl-list (buffer-name)))

   ;; Remove help window
   (if (and (not option) gp-tutorial-requiredp)
     (progn
       (gp-window-manager "*gp-help*" 'gp-remove-help-old-config)
       (gp-restore-wind-conf))
   )))

(defun gp-make-cpl-file nil
  (interactive)
  (gp-prepare-cpl-file nil))

(defsubst gp-cpl-file-menu nil ""
  (nconc
   (list
    (vector (gp-messager 65) 'gp-cpl-file ':active t ':key-sequence nil)
    (vector (gp-messager 66) 'gp-edit-cpl-file
            ':active t ':key-sequence nil))
   (if (eq major-mode 'gp-script-mode)
       (list (vector (gp-messager 67) 'gp-make-cpl-file
                     ':active t ':key-sequence nil)
             ["Info" gp-cpl-file-info :active t :key-sequence nil
              :included gp-tutorial-requiredp])
     nil)))

(add-hook 'pari-mode-hook
  '(lambda nil
     (if (file-really-exists-p (concat (buffer-name) ".cpl"))
         ;; The local completion for this file.
         (gp-cpl-file (concat (buffer-name) ".cpl")))
     (if (file-really-exists-p gp-additional-cpl-file)
         ;; Add this file to the usual completion array.
         (gp-cpl-file gp-additional-cpl-file))
     (define-key gp-map "\t"          (function gp-complete))
     (define-key gp-map "\M-i"        (function gp-complete))
     (define-key gp-map "\M-\t"       (function gp-complete)) ;; C-i in fact !!!
     (define-key gp-script-map "\M-i" (function gp-complete))))

(add-hook 'gp-mode-hook
  '(lambda nil (gp-choose-complete)))

(add-hook 'menu-bar-update-hook
  '(lambda nil
     (when (memq major-mode '(gp-mode gp-script-mode))
       (easy-menu-change '("GP") (gp-messager 64) (gp-cpl-file-menu) (gp-messager 71))
       (when GP-script-menu-map
         (easy-menu-add-item GP-script-menu-map nil
            (vector (gp-messager 72) 'gp-complete t) (gp-messager 73)))
       (when GP-menu-map
         (easy-menu-add-item GP-menu-map nil
            (vector (gp-messager 72) 'gp-complete t) (gp-messager 73))))))


;; pari-completion.el ends here ---------------
