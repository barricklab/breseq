;; $Id: sli-tools.el 6643 2005-02-20 18:51:46Z kb $
;; sli-tools.el --- structured languages indentation package

;; It works out some tools for indentation of structured programs.
;; It has been written for mupad.el and pari.el but should apply to
;; any other structured language like Pascal.
;; See sli-tools and sli-structures below.

;; The way it works inside:
;; sli-tell-indent is the main engine. They are two cases
;; either we want to indent the line the cursor is on,
;; or we want determine the indent of the next line.
;; See also sli-forward-sexp.

;;   BASICS FROM SLI-STRUCTURES:
;;     you should read the information concerning this variable,
;;     but some basics are required to go further.
;;     In the construct
;;             if toto then tata
;;             else
;;               titi 
;;             end_if
;;     "if" is called a HEAD or a head-key,
;;     "else" is called a STRONG, 
;;     "end_if" is called an END.
;;     Basically, the "else" is aligned on the "if" and the
;;     "end_if" on previous "else"/"elif"/"if" if the HEREDITY applies.
;;     HEREDITY applies unless otherwise specified.
;;     The key "then" is called a SOFT-key: it implies special
;;     indentation afterwards but is not aimed at being under
;;     the "if".
;;     Keys can also be termed
;;          FIXED (usually global stuff),
;;          BEACON (like "do" in a while-construct),
;;          RELATION (math),
;;          SEPARATOR,
;;          CONSTRUCTOR,
;;          SPECIAL-HEAD (initial declarators like "local", "var", "remember").

;;     The same END can be used for several HEADs, a word can be
;;     a HEAD and a SPECIAL-HEAD, but if so, its corresponding HEAD name
;;     cannot be its own.
;;     HEADs, RELATIONs, BEACONs, SEPARATORs, SOFTs, ENDs should all be different,
;;     and SPECIAL-HEADs can only be also HEADs.
;;     SOFTs, STRONGs or ENDs can be used in fifferent constructs.

;;     SYNTACTICALLY speaking, chars used in these strings should be word-constituents,
;;     symbols, open-parenthesis, close-parenthesis or generic-parenthesis.
;;     If sli-case-fold is t, upper/lowercase letters are irrelevant *but*
;;     sli-structures and all should use lowercase letters.

;;   INDENT OF THIS LINE:
;;     we look if the first word on this line is a fixed/strong/end/soft
;;        if yes --> fixed keys are easy
;;               --> soft keys: find its ancestor (a strong or a head)
;;                   this ancestor is necessarily on another line,
;;                   so compute the indent required after this key.
;;               --> strong/end keys: find its ancestor and align
;;                   our key on the ancestor (strong or head), with possible offset.
;;                   If the attribute is 'absolute, apply this indent.
;;                   Else, apply it except if this key belongs to sli-no-heredity-list,
;;                   in which case the alignment is on the head.
;;                
;;        if no  --> use indentation of previous line
;;   INDENT OF NEXT LINE:
;;     see if previous line has an unclosed head/strong/soft.
;;        if yes --> use its indentation.
;;        if no  --> use indentation of previous line.
;;   SEE sli-tools for more info.

;;  REGION scanned: the region scanned is extremely important for lengthy programms,
;;  since no unclosed constructs may be found before the very beginning of the file.
;;  So we provide the variable `sli-safe-place-regexp' which indicates where one
;;  can start: after the end of the first grouping. For inctance
;;  "^\(\\\\--\)$" means that a line containing only "\\--" indicates a place
;;  outside any construct. One can start after the string "--" or before the "\\".

;;  COMMENTS: nothing much has been done for indenting comments just now.

;; Use of properties:
;;  -- 'sli-type can be
;;          'head 'special-head 'strong 'soft 'end 'math-relation 'beacon
;;          'block-comment-start 'block-comment-end 'string
;;  -- 'sli-ancestor if present is a buffer location:
;;          for 'end  it is point at beginning of opening 'head or an intermediate 'strong
;;          for 'strong it is point at beginning of corresponding 'head
;;          for 'special-head it is point at beginning of previous 'special-head or 'head
;;          for 'block-comment-end it is point at beginning of corresponding 'block-comment-start
;; -- 'sli-reverse-ancestor if present is a buffer location:
;;          for 'head it is point at beginning of closing 'end              *Not Always Present*
;;          for 'strong it is point at beginning of next 'strong or 'end    *Not Always Present*
;;          for 'special-head it is point at beginning of closing separator *Not Always Present*
;;          for 'block-comment-start it is point at beginning of corresponding 'block-comment-end
;; -- 'sli-time if present is an integer representing the time when
;;          the sli-properties were last set.
;; These properties are lazily computed: everytime we can deduce such a property,
;; we do it, but we do not go out of our way to do so. So the absence of a property
;; only means it has not been computed, and *not* it doesn't exist.

;; Maintainer: Olivier Ramare <ramare@agat.univ-lille1.fr>

;; BUGS:
;; (1) If I remember well, strings spreading over several lines may
;;     raise some troubles.
;; (2) sli-tutor has some troubles if used in the middle of already
;;     complete structures.
;; (3) Due to lazy computations of text properties, sli-show-sexp may
;;     show wrong things. Wait a bit and things will become ok.
;;     See `sli-prop-do-not-recompute-time'.
;; Use of sli-special-head-heads-alist ??

(provide 'sli-tools)

;;------------------------------------------------------
;; Variables that defines how indentation should occur.
;; See mupad.el for an example.
;;------------------------------------------------------

;; We use "" and  \" for strings.

(defgroup sli nil
"sli customization group"
:group 'editing :prefix "sli-")

(defcustom sli-handles-sexp nil "A true value advises forward/backward/scan-sexp/s"
:type 'boolean :group 'sli)

;; These values are modified in sli-tools:
(defvar sli-verbose nil "A true value gives (debugging) infos")
(defvar sli-prop-verbose nil "A true value gives (debugging) infos on text properties")

(eval-and-compile
;; The next variables are here to pacify the compiler !
;; Do *not* assign any value to them or they may override ....
(defvar block-comment-end)
(defvar block-comment-start))

(defvar sli-structures nil
  "List of lists. Each item is a vector or a list which we call a STRUCTURE
in this explanation. There are several kind of structures :

([HEAD-STRING head INDENT-HEAD]
 [SOFT-STRING1 soft INDENT-SOFT1]
 ([STRONG-STRING1 strong INDENT-STRONG1]
  [SOFT-STRING2 soft INDENT-SOFT2])
 ([STRONG-STRING2 strong INDENT-STRONG2])
 [END-STRING end])
is the usual structure, like in 'if/then/(elif/then)/(else)/end_if'.  Between
the 'head' and the 'soft', INDENT-HEAD is used on subsequent lines to offset the
new line with respect to the beginning of HEAD-STRING. When the 'soft' is found,
INDENT-SOFT1 is used still with respect to the 'head'.  The next part is
optional.  The STRONG-STRING is aligned on its 'head' and INDENT-STRONG is used
after that, with respect to the STRONG-STRING. Finally the END-STRING is aligned
on the previous STRONG-STRING (the 'heredity principle'). If you want to change
this alignement, use `sli-shift-alist' below.  Note that an INDENT-* value can
be either an integer or a cons pair whose first element is the symbol 'absolute
and the second one is an integer: it means that the indentation is not relative
but absolute with respect to the left margin. It applies also to the next
strong/end key.  In this construct, you can also use [SPECIAL-HEAD-STRING
special-head INDENT-SPECIAL-HEAD SEPARATORS]. This key is closed by SEPARATORS
which is either a separator which belongs to `sli-separators' or a list of
separators all in `sli-separators' in which case the first one is the one used
by sli-maid. No other construct should happen between the special-head and its
separator except comments and keys termed CONSTRUCTORs; for instance the
'proc/(option)/begin/end_proc' construct of MuPAD is
a head/special-head/strong/end. You can use several [END-STRING end]. The first
one is going to be used by the maid. Furthermore you can use the same END-STR
for several constructs. It then applies to the first 'head' that appears
(going backward). Concerning SPECIAL-HEAD, the syntax could make believe that
a string could be used after a HEAD with some separators and after another one
with some other separators: in fact they are merge internally so the union
of all appearing separators for this SPECIAL-HEAD is being used.

([BEACON-STRING beacon INDENT-BEACON]) specifies a special string that can be
found between a 'head' or a 'strong' and its corresponding 'soft'. The typical
example being 'for t from 1 to 2 do' and has pattern
'head/beacon/beacon/soft'. If a newline is asked after the 'from' but before the
'to', indentation is done with respect to the beginning of 'from' and
INDENT-BEACON is added except if this newline is asked just after the beacon
key, in which case indentation is done like from before the beacon but
'math-relation's are ignored. Simply because 'math-relation' are supposedly
closed by the appearance of a beacon, whether a separator has occured or not.

([RELATION-STRING math-relation INDENT-RELATION]) specifies a mathematical type
of relation (like '='). Such operators acts either as beacons (example 'while
t=3D55 do' with pattern 'strong/math-relation/soft') or else are closed by
someone in `sli-separators'. They may contain further structures in between like
in 'foo = if ok then gonethrough=t ; 3 else 5 end_if'.  INDENT-RELATION is used
before the appearance of the proper separator.

HEAD-STRINGs, MATH-RELATION-STRINGs, BEACON-STRINGs, SEPARATORs should all be
different, except one case for HEAD-STRINGs indicated below.  SOFT-STRINGs and
STRONG-STRINGs are different from any of the above, but a same soft or strong
key can be used in different constructs. Usual examples are 'then' and 'do' and
the 'elif' in 'if/elif/end_if' and '%if/elif/end_if'.  But because of the way
things are, the corresponding INDENT should be the same throughout. Note that
longest match is always taken, so that if 'while(' is a head (like in gp) and
'(' is also a head (almost everywhere), indentation after 'while(' is the one it
should. Same applies for the two constructs '%if' and 'if' in mupad.

Concerning HEAD-STRINGs, all starting heads are to be distincts, but inside a
construct, an existing head can be used as a special head. The typical case in
MuPAD is 'category' which is normally a head but can be used like a special head
inside a 'domain' statement.

CONSTRUCTORs are treated in a special way and keys declared as head or end
or whatever can also be termed constructor. Usual example: ( is a head and
is also declared as a constructor.

Cdr's are to be evaled.

If downcase/uppercase is relevant is controled by the variable `sli-case-fold'.
If sli-case-fold is t, sli-structures should use lowercase letters.

Technical note: the first element of this list *has to* contain a 'head'. ")

(defvar sli-case-fold nil
"The strings used as separators, relations, and all. Not yet used.
If set to t, all keywords in sli-stryctures, sli-shift-alist ...
should be in lowercase.")

(defvar sli-escape-key-string ""
"The strings used as separators, relations, and all. Not yet used.")

(defvar sli-shift-alist nil
"Usual 'strong/end' are aligned on the previous
occurence of a corresponding head/strong.
You can add an offset between two keys.
This is also valid in case of an absolute indent.
Elements of this list have format ([key1 key2] . offset).
Cdr's are to be evaled.")

(defvar sli-no-heredity-list nil
"Usual 'strong/end' are aligned on the previous
occurence of a corresponding head/strong except
if mentionned in this list.
Elements of this list have format [head-key key].")

(defvar sli-separators nil "Do not forget `sli-is-a-separatorp'.")

(defvar sli-is-a-separatorp-fn 'sli-is-a-separatorp-default
  "Function called to decide if character after POINT
is a separator. This function takes an optional argument
which is the value of POINT and should be surrounded by
save-excursion and save-match-data, see `sli-is-a-separatorp-default'.")

(defun sli-is-a-separatorp-default (&optional pt)
  (save-excursion
    (when pt (goto-char pt))
    (save-match-data
      (if sli-separators
          (let ((case-fold-search sli-case-fold))
            (looking-at (regexp-opt sli-separators)))
        nil))))

(defun sli-is-a-separatorp (&optional pt)
  (funcall sli-is-a-separatorp-fn pt))

(defvar sli-put-newline-fn 'sli-put-newline-default
"Function used to insert a newline. Takes no argument.")

(defun sli-put-newline-default nil (insert-char ?\n 1))

(defun sli-put-newline nil
"Indirection. Puts a newline according to `sli-put-newline-fn'
and takes care not to write anything on read-only parts."
  (unless (get-text-property (point) 'read-only)
    (funcall sli-put-newline-fn)))

(defvar sli-safe-place-regexp "^\\(//--+\\|/\\*-+-\\*/\\)$"
"Marker used to tell emacs this point is outside a commented area, a string or a sexp. The safe place starts at beginning of match-group 1 and ends at end of match-group 1.")

(defvar sli-fixed-keys-alist '()
"Some keys should be placed at a fixed place with respect to the
indentation of previous line when following a RELATION sign. See
`sli-relation-keys'. This is the corresponding alist.
List of (STRING . INDENTATION).")

(defvar sli-keys-with-newline nil
"When `sli-maid' tries to further your constructs, some keys should be
followed by a newline before completion is added.")

(defvar sli-keys-without-newline nil
"When `sli-maid' tries to further your constructs, some keys should never be
followed by a newline.")

(defvar sli-maid-correction-alist nil "See `sli-maid'")

(defvar sli-add-to-key-alist nil "See `sli-maid'.")

(defvar sli-more-maidp t "See `sli-maid'.")

(defvar sli-tab-always-indent t "See `sli-electric-tab'.")

(defvar sli-comment-starts '()
"A list of possible starters of one-line comments.
That is to say an extension of `comment-start' in this special case.")

(defvar sli-block-comment-middle-offset -1
"Indentation of block comments: they start with block-comment-start and then
either some whitespace and a word on the same line, on which case next lines
are aligned on this first word. Or the text starts on next line in which case
they start at column-of-end-of-block-comment-start + this-offset.
Exception for the last line if it contains only one word ending with
'block-comment-end in which case this word is where placed at
column-of-end-of-block-comment-start+sli-block-comment-end-offset spaces
from the margin.")

(defvar sli-block-comment-end-offset -1
"See `sli-block-comment-middle-offset'.")


;;;--------------------------------------------------------------------------
;;; Inner variables
;;;--------------------------------------------------------------------------

(defvar sli-head-keys nil)
(defvar sli-special-head-keys nil)
(defvar sli-soft-keys nil)
(defvar sli-beacon-keys nil)
(defvar sli-math-relation-keys nil)
(defvar sli-relation-keys nil)
(defvar sli-constructor-keys nil)
(defvar sli-keys-nomrelations nil) ; nomrelations means no-math-relations
(defvar sli-strong-keys nil)
(defvar sli-end-keys nil)
(defvar sli-keys nil)
(defvar sli-max-keys-length 0
"An integer: the maximum length of a keyword in sli-structures.
Used in `sli-anchored-posix-search-backward', a fix for `posix-search-backward'. ")
(defvar sli-all-keys-nomrelations-noseparators-regexp nil)
(defvar sli-all-keys-regexp nil) ; including string quotes and all kind of comments.
(defvar sli-all-end-strong-regexp nil)
(defvar sli-fixed-regexp nil)
(defvar sli-head-regexp nil)
(defvar sli-strong-regexp nil)
(defvar sli-all-keys-and-constructors-regexp nil)

(defvar sli-head-end-alist nil "The alist ((end . head) ...).")
(defvar sli-ends-head-alist nil "The alist ((head . (end1 end2 ...) ...).")
(defvar sli-heads-strong-alist nil "The alist ((strong . (head1 head2 ...)) ...).")
(defvar sli-special-head-alist nil "The alist ((special-head . (separator1 separator2 ...)) ...).")
(defvar sli-special-head-heads-alist nil
  "The alist ((special-head . heads) ...) for those special heads that are also heads.")
(defvar sli-special-head-previous-keys-alist nil
  "The alist ((special-head . keys) ...) for special-heads that can be heads.
keys are the keys that can be before special-head.")
(defvar sli-companion-strong-keys-alist nil
  "The alist  ((strong/head . (strongs that could be after)) ...).
The car should be a member of the cdr if the car is a strong.")
(defvar sli-soft-alist nil 
  "The alist ((ambiguous-soft . (head-or-strong1 head-or-strong2 ...)) ...).")
(defvar sli-soft-head-or-strong-alist nil "The alist ((head-or-strong . soft) ...)")
(defvar sli-first-offset-alist nil)  ; to apply before the soft
        ; it applies to head/strong keys that are followed by a soft with no
        ; head or strong in between. Morally speaking this soft "closes" the head/strong.
(defvar sli-relevant-alist nil
"An alist. Put all head/strong/end's in one bundle. say two keys are linked if
they occur in a same constructs. Close this relation transitively.
this is the alist ((key . (keys in the same class)) ...).")
(defvar sli-ancestors-alist nil
"The alist ((end/strong-key . (head/strong1 head/strong2 ...)) ...)
of keys that can occur before the first key.")

(defvar sli-second-offset-alist nil "Alist (key . offset) where
OFFSET is the one to apply after the soft key if it exist, after
KEY if it doesn't have any soft. KEY can be a head/end/strong/soft.")  ; to apply after the soft
(defvar sli-special-head-offset-alist nil "Alist (special-head . offset).")
(defvar sli-relation-offset-alist nil)

(defvar sli-maid-alist nil)
(defvar sli-ambiguous-keys nil
  "List of keys that may ask for a different following key according
to context. They *should be* soft or strong keys.")

;; Only to shut up compiler. These two variables should be defined when the
;; correct buffer is set ! Used by sli-show-sexp.
(defvar sli-overlay-beg nil "overlay set by `sli-show-sexp' and showing the head key.")
(defvar sli-overlay-end nil "overlay set by `sli-show-sexp' and showing the end key.")

(defvar sli-prop-do-not-recompute-time 10
"Time span in milliseconds under which it is not necessary to recompute
text properties alloted by sli-tools.")
(defvar sli-prop-used 0
"Number of times text-properties have been used.")
(defvar sli-key-is-a-special-headp nil
  "Set by `sli-get-corresponding-key' and `sli-get-first-non-end-key'.")

(mapc 'make-variable-buffer-local 
'(sli-verbose sli-prop-verbose sli-handles-sexp sli-overlay-beg sli-overlay-end
sli-prop-do-not-recompute-time sli-structures sli-shift-alist sli-separators
sli-is-a-separatorp-fn sli-more-maidp sli-add-to-key-alist
sli-math-relation-keys sli-max-keys-length sli-no-heredity-list sli-head-keys
sli-special-head-keys sli-soft-keys sli-beacon-keys sli-relation-keys
sli-keys-nomrelations sli-strong-keys sli-end-keys sli-keys sli-prop-used
sli-all-keys-nomrelations-noseparators-regexp sli-all-keys-regexp sli-all-end-strong-regexp
sli-soft-head-or-strong-alist sli-head-end-alist sli-heads-strong-alist
sli-special-head-alist sli-special-head-heads-alist
sli-special-head-previous-keys-alist sli-ends-head-alist sli-head-regexp
sli-strong-regexp sli-relevant-alist sli-ancestors-alist sli-fixed-keys-alist
sli-fixed-regexp sli-companion-strong-keys-alist sli-soft-alist
sli-first-offset-alist sli-second-offset-alist sli-relation-offset-alist
sli-maid-alist sli-ambiguous-keys sli-constructor-keys sli-all-keys-and-constructors-regexp
sli-block-comment-middle-offset sli-block-comment-end-offset sli-key-is-a-special-headp
sli-special-head-offset-alist))

;;;-----------------------------------------------------------------------------
;;; This section is devoted to some precomputations from sli-structures.
;;; Lots of work is done several time, but I prefer this modularity
;;; since it is easier to modify.
;;;-----------------------------------------------------------------------------

(defun sli-split-list (lst)
  (let ((wordother '()) (otherword '()) (wordword '()) (otherother '()) ls)
    (mapc
     (lambda (wd)
       (setq ls (string-to-list wd))
      (cond
       ((and (= (char-syntax (car ls)) ?w) (= (char-syntax (car (last ls))) ?w))
        (add-to-list 'wordword wd))
       ((= (char-syntax (car ls)) ?w)
        (add-to-list 'wordother wd))
       ((= (char-syntax (car (last ls))) ?w)
        (add-to-list 'otherword wd))
       (t (add-to-list 'otherother wd))))
     lst)
    (list wordword wordother otherword otherother)))

(defun sli-regexp-opt (lst)
  (let ((qlst (sli-split-list lst)))
    (if (null (elt qlst 0))
        (if (null (elt qlst 1))
            (if (null (elt qlst 2))
                (if (null (elt qlst 3))
                    "\\<\\>"
                  (regexp-opt (elt qlst 3) t)) ; grouping required for posix
              (concat
               (regexp-opt (elt qlst 2) t) "\\>"
               (if (null (elt qlst 3))
                   ""
                 (concat "\\|" (regexp-opt (elt qlst 3) t)))))
          (concat
            "\\<" (regexp-opt (elt qlst 1) t)
            (if (null (elt qlst 2))
                (if (null (elt qlst 3))
                    ""
                  (concat "\\|" (regexp-opt (elt qlst 3) t)))
              (concat
               "\\|" (regexp-opt (elt qlst 2) t) "\\>"
               (if (null (elt qlst 3))
                   ""
                 (concat "\\|" (regexp-opt (elt qlst 3) t)))))))
      (concat
       "\\<" (regexp-opt (elt qlst 0) t) "\\>"
       (if (null (elt qlst 1))
            (if (null (elt qlst 2))
                (if (null (elt qlst 3))
                    ""
                  (concat "\\|" (regexp-opt (elt qlst 3) t)))
              (concat
               "\\|" (regexp-opt (elt qlst 2) t) "\\>"
               (if (null (elt qlst 3))
                   ""
                 (concat "\\|" (regexp-opt (elt qlst 3) t)))))
          (concat
            "\\|\\<" (regexp-opt (elt qlst 1) t)
            (if (null (elt qlst 2))
                (if (null (elt qlst 3))
                    ""
                  (concat "\\|" (regexp-opt (elt qlst 3) t)))
              (concat
               "\\|" (regexp-opt (elt qlst 2) t) "\\>"
               (if (null (elt qlst 3))
                   ""
                 (concat "\\|" (regexp-opt (elt qlst 3) t)))))))))))

(defun sli-flatten (ls)
  (let ((res '()))
    (mapc
      (lambda (ph)
        (cond
          ((listp ph) (setq res (append res (sli-flatten ph))))
          (t (setq res (append res (list ph))))))
      ls)
    res))

(defun sli-scan-structures-locally (stru symbol)
  (let ((res '()))
    (mapc (lambda (ph)
              (setq res
                (append res
                  (cond
                    ((listp ph) (sli-scan-structures-locally ph symbol))
                    ((equal (elt ph 1) symbol) (list (elt ph 0)))
                    (t '())))))
            stru)
    res))

(defsubst sli-compact-list (lst)
  ; remove same consecutive occurences.
  (let* ((old (car lst)) (nlst (list old))  (lst (cdr lst)))
    (while lst
      (if (equal (car lst) old)
          (setq lst (cdr lst))
          (setq nlst (cons (setq old (car lst)) nlst) lst (cdr lst))))
    (nreverse nlst)))

(defun sli-scan-structures (symbol)
  (let ((res '()))
    (mapc
      (lambda (st)
        (when (equal (elt st 1) symbol)
          (add-to-list 'res (elt st 0))))
      (sli-flatten sli-structures))
  res))

(defun sli-get-ends-head-alist nil
  (let ((res '()) all-ends) ; forme la liste (head-key . (end1 end2 ...))
   (mapc
     (lambda (ph)
       (when (equal (elt (elt ph 0) 1) 'head)
         (setq all-ends '())
         (mapc
           (lambda (s)
             (when (and (vectorp s) (equal (elt s 1) 'end))
               (setq all-ends (append all-ends (list (elt s 0))))))
          ph)
         (add-to-list 'res (cons (elt (elt ph 0) 0) all-ends))))
     sli-structures)
   res))

(defun sli-get-head-end-alist nil
  (let ((res '()) all-heads) ; forme la liste (end-key . (head1 head2 ...))
   (mapc
     (lambda (end)
       (setq all-heads '())
       (mapc
	(lambda (s)
	  (if (member end (cdr s))
	      (add-to-list 'all-heads (car s))))
	sli-ends-head-alist)
       (add-to-list 'res (cons end all-heads)))
   sli-end-keys)
   res))

(defun sli-get-strong (ph)
  (let ((res '()))
    (mapc
      (lambda (st)
        (when (equal (elt st 1) 'strong)
          (add-to-list 'res (elt st 0))))
      ph)
  res))

(defun sli-get-heads-strong-alist nil
  (let ((res '()) (aux '()) possible-heads) ; forme la liste des (strong-key . (head-key1 head-key2 ...))
   ; Peut-etre plusieurs strong pour chaque head.
   (mapc
     (lambda (ph)
       (if (equal (elt (elt ph 0) 1) 'head)
           (let ((strongs (sli-get-strong (sli-flatten ph))))
              (unless (null strongs)
                 (mapc (lambda (st)
                         (setq aux (add-to-list 'aux
                                                (cons st (elt (elt ph 0) 0)))))
                       strongs)))))
     sli-structures)
   ; Une strong peut etre liee a plusieurs heads. Il faut les reunir:
   (mapc
    (lambda (strong)
      (setq possible-heads '())
      (mapc
       (lambda (ajoint)
         (when (equal (car ajoint) strong)
           (setq possible-heads (append possible-heads (list (cdr ajoint))))))
       aux)
      (when (> (length possible-heads) 1)
        (add-to-list 'sli-ambiguous-keys strong))
      (setq res (append res (list (cons strong possible-heads)))))
    (sli-compact-list (sort (mapcar 'car aux) 'string-lessp)))
   res))

(defun sli-get-soft-alist nil ; forme la liste (soft . (head of strong using it))
  (let ((resaux '()) loc (res '()) astrong-list (asoft-list '()))
   (mapc
    (lambda (ph)
      (setq astrong-list '())
      (mapc
       (lambda (ve)
	 (cond
          ((equal (elt ve 1) 'soft) (unless (null astrong-list)
                                      (add-to-list 'resaux (cons (elt ve 0) astrong-list))
                                      (add-to-list 'asoft-list (elt ve 0))))
          ((member (elt ve 1) '(strong head)) (setq astrong-list (list (elt ve 0))))))
       (sli-flatten ph)))
    sli-structures)
    ;; now gather identical soft:
    (mapc
      (lambda (asoft)
        (setq loc '())
        (mapc
	  (lambda (dd)
	    (when (string-equal asoft (car dd))
              (setq loc (append loc (cdr dd)))))
          resaux)
        (add-to-list 'res (cons asoft (sli-compact-list (sort loc 'string-lessp)))))
      asoft-list)
    res
   ))

(defun sli-common-pointp (l1 l2)
  "t if l1 and l2 have a common element. Test is done through member."
  (let ((ok nil))
    (mapc (lambda (c) (setq ok (or ok (member c l1)))) l2)
    ok))

(defun sli-get-companion-alist nil ; case ?? It was not there.
  (let ((res '()))
    ; on prend les car de sli-heads-strong-alist on leur
    ; associe la liste des car qui ont au moins une tete en commun :
    (mapc
      (lambda (co)
        (let ((end (cdr co)) (companions '()))
          (mapc
            (lambda (coo)
               (when (sli-common-pointp (cdr coo) end)
                 (setq companions (add-to-list 'companions (car coo)))))
            sli-heads-strong-alist)
          (setq res (append res (list (cons (car co) companions))))))
      sli-heads-strong-alist)
    ; on prend les cdr de sli-heads-strong-alist on leur
    ; associe la liste des car possibles :
    (mapc
      (lambda (head)
        (let ((companions '()))
          (mapc
            (lambda (coo)
               (when (member head (cdr coo))
                 (setq companions (add-to-list 'companions (car coo)))))
            sli-heads-strong-alist)
          (setq res (add-to-list 'res (cons head companions)))))
      (sli-compact-list (sort (sli-flatten (mapcar 'cdr sli-heads-strong-alist)) 'string-lessp)))
    res))

(defun sli-get-soft-head-or-strong-alist nil
  (let ((res '()) asoft astrong-list)
    (mapc
     (lambda (ass)
       (setq asoft (car ass))
       (setq res (append res (mapcar (lambda (st) (cons st asoft)) (cdr ass)))))
     sli-soft-alist)
    res))

(defun sli-equivalence-classes-local (lst)
  (cond
   ((null lst) lst)
   (t (let (lstbis (done nil) (l1 (car lst)))
	(setq lstbis
	      (mapcar
	       (lambda (c)
		 (if (sli-common-pointp l1 c)
		     (progn
		       (setq done t)
		       (sli-compact-list (sort (append l1 c) 'string-lessp)))
		   c))
	       (sli-equivalence-classes-local (cdr lst))))
	(unless done
	  (setq lstbis (append lstbis (list l1))))
	lstbis))))

(defun sli-equivalence-classes (lst)
  (while (> (length lst) (length (setq lst (sli-equivalence-classes-local lst)))))
  lst)

(defun sli-get-relevant-alist nil
  (let (key-lst (res '()))
    ;; relevant keys are head/strong or end keys.
    (mapc
     (lambda (class)
       (mapc
	(lambda (el)
	  (add-to-list 'res (cons el class)))
	class))
     (sli-equivalence-classes
      (delq nil ; nil had better not be the first one ...
	    (mapcar
	     (lambda (ph)
	       (setq key-lst '())
	       (mapcar
		(lambda (co)
		  (when (member (elt co 1) '(head strong end))
		    (add-to-list 'key-lst (elt co 0))))
		ph)
	       key-lst)
	     (mapcar 'sli-flatten sli-structures)))))
    res))

(defun sli-get-ancestors-alist nil
  (append
   ;; Ancestors for end-keys:
   (mapcar
    (lambda (end)
      (cons end
	    (sli-flatten
	     (mapcar
	      (lambda (head)
		(or (assoc head sli-companion-strong-keys-alist) ; works only if a strong is present
		    (cdr (assoc end sli-head-end-alist))))
	      (cdr (assoc end sli-head-end-alist))))))
    sli-end-keys)
   ;; Ancestors for strong-keys:
   (mapcar
    (lambda (strong)
      (cons strong
            (append (cdr (assoc strong sli-heads-strong-alist))
                    ;; The next one is bad: for "begin" it associates "begin" which
                    ;; can not be an anscestor ...
                    (cdr (assoc strong sli-companion-strong-keys-alist)))))
    sli-strong-keys)))

(defun sli-get-first-offset-alist nil
  (let ((res '()) last-head-or-strong stru pl)
    (mapc
     (lambda (ph)
       (setq last-head-or-strong nil stru (sli-flatten ph))
       (while (not (null stru))
         (setq pl (car stru))
         (cond
           ((member (elt pl 1) '(head strong)) (setq last-head-or-strong pl))
           ((equal (elt pl 1) 'soft)
            (when last-head-or-strong
              (setq res (append res (list (cons (elt last-head-or-strong 0)
                                                (elt last-head-or-strong 2))))
                    last-head-or-strong nil))))
           (setq stru (cdr stru))))
     sli-structures)
    res))

(defun sli-get-second-offset-alist nil
  (let ((res '()) last-cand stru pl)
    (mapc
     (lambda (ph)
       (setq last-cand nil stru (sli-flatten ph))
       (while (not (null stru))
         (setq pl (car stru))
         (cond
           ((equal (elt pl 1) 'head)
            (setq last-cand pl))
           ((and (member (elt pl 1) '(end strong))
		 (not (assoc (elt pl 0) sli-special-head-heads-alist))) ;; ???
            (when last-cand ;; no soft after last-cand.
              (setq res (append res (list (cons (elt last-cand 0)
                                                (elt last-cand 2))))))
            (if (equal (elt pl 1) 'end)
                (setq last-cand nil)
              (setq last-cand pl)))
           ((equal (elt pl 1) 'soft)
            (when last-cand ;; last-cand is followed by a soft
              (setq res (append res (list (cons (elt last-cand 0)
                                                (elt pl 2))))
                    last-cand nil))))
	 (setq stru (cdr stru))))
     sli-structures)
    res))

(defun sli-get-relation-offset-alist nil
  (let ((res '()))
    (mapc
      (lambda (ph)
        (mapc
	 (lambda (pl)
	   (cond
	    ((member (elt pl 1) '(math-relation beacon))
	     (add-to-list 'res (cons (elt pl 0) (elt pl 2))))))
	 ph))
      (mapcar 'sli-flatten sli-structures))
    res))

(defun sli-get-special-head-offset-alist nil
  (let ((res '()))
    (mapc
      (lambda (ph)
        (mapc
	 (lambda (pl)
	   (cond
	    ((member (elt pl 1) '(special-head))
	     (add-to-list 'res (cons (elt pl 0) (elt pl 2))))))
	 ph))
      (mapcar 'sli-flatten sli-structures))
    res))

(defun sli-get-maid-alist-locally (ph lst)
  (let ((res '()) aux resaux (nlst '()))
    (cond
      ((null ph))
      ((listp (car ph))
       (setq ; process the internal with no 'lst' since it is optional:
	     aux (sli-get-maid-alist-locally (car ph) '())
             ; Then process the remainder with both candidates 'lst' and (cadr aux):
             resaux (sli-get-maid-alist-locally (cdr ph) (append (cadr aux) lst))
             ; glue things together:
             res (list (append aux (car resaux)) (cadr resaux))))
      (t (setq aux (elt (car ph) 0) ; the new 'last-word (lst=(last-word))
               ph (cdr ph))
         ; Link 'lst' to the new compulsory:
         (mapc (lambda (s) (add-to-list 'res (cons s aux))) lst)
	 (while (and (not (null ph)) (listp (car ph)))
           ; (car ph) is an optional construct. Scan it with no 'lst'
           (setq resaux (sli-get-maid-alist-locally (car ph) '())
                 ; gather all 'last-words':
                 nlst (append nlst (cadr resaux))
                 ; gather all bindings :
                 res (append res (car resaux))
		 ph (cdr ph)))
	 (when (car ph) ; aux is linked to the new guy:
	   (add-to-list 'res (cons aux (elt (car ph) 0)))
           ; the new guy is linked with all the 'last-words':
	   (mapc(lambda (s) (add-to-list 'res (cons s (elt (car ph) 0)))) nlst))
         ; process things farther:
         (setq resaux (sli-get-maid-alist-locally ph '())
               res (list (append (car resaux) res)
                         (if (null (cadr resaux)) (append (list aux) nlst)
                             (cadr resaux))))))
     res))

(defsubst sli-full-stuff (key alist fn1 fn2)
  (let ((res '()) aux)
    (while alist
      (when (setq aux (funcall fn1 (funcall fn2 key alist)))
        (add-to-list 'res aux))
      (setq alist (cdr alist)))
    res))

(defsubst sli-full-assoc (key alist)
  "The list of cdrs in alist whose car is key."
  (sli-full-stuff key alist 'cdr 'assoc))

(defsubst sli-full-rassoc (key alist)
  "The list of cars in alist whose cdr is key."
  (sli-full-stuff key alist 'car 'rassoc))

(defun sli-get-automatic-maid-alist nil
;; sli-ambiguous-keys is also created here.
  ;(setq sli-ambiguous-keys nil)
  (let ((res '()))
    (mapc
      (lambda (ph)
        (setq res (append res (car (sli-get-maid-alist-locally ph '())))))
      sli-structures)  ;(princ "\n") (princ (list "sli-get-automatic-maid-alist" res))
    (add-to-list 'res (cons block-comment-start block-comment-end))
    ; well, soft keys may correspond to different strong keys...
    (mapcar (lambda (co) (let ((to (sli-full-assoc co res)))
                           (cons co (if (null (cdr to)) (car to)
                                      (progn
                                        (add-to-list 'sli-ambiguous-keys co)  to)))))
          (sli-compact-list (sort (mapcar 'car res) 'string-lessp)))))

(defun sli-get-maid-alist nil
  ;; First, create the list automatically:
  (setq sli-maid-alist (sli-get-automatic-maid-alist))
  ;(princ "\n") (princ (list "sli-get-maid-alist" sli-maid-alist))
  ;; But now users may want something else. A typical example is
  ;; for-from-do-end_for where the proposed completion of "for"
  ;; is "do" because "from" is only a beacon.
  ;; Correction is done is two steps: first the elements who have
  ;; a car is sli-maid-correction-alist are removed from
  ;; from sli-maid-alist and then sli-maid-correction-alist
  ;; is added.
  (let ((new-lst '()) (correction-words (mapcar 'car sli-maid-correction-alist)))
    (while sli-maid-alist  ;(princ "\n")  (car sli-maid-alist)
      (unless (member (caar sli-maid-alist) correction-words)
        (setq new-lst (append new-lst (list (car sli-maid-alist)))))
      (setq sli-maid-alist (cdr sli-maid-alist)))
    (append new-lst sli-maid-correction-alist)))

(defun sli-get-special-head-alist nil
  (let ((res '()) aux)
    (mapc
     (lambda (ph)
       (if (equal (elt ph 1) 'special-head)
           (progn
             (if (setq aux (assoc (elt ph 0) res))
                 ;; This special-head has already been used, but maybe with
                 ;; different separators. Merge everything ... Sorry !
                 (progn
                   (setq res (delq aux res));(print res)
                   (setq aux (cdr aux))
                   (mapc (lambda (wd) (add-to-list 'aux wd))
                         (if (listp (elt ph 3)) (elt ph 3)(list (elt ph 3))))
                   (add-to-list 'res (cons (elt ph 0) aux)))
               (add-to-list 'res (cons (elt ph 0)
                                       (if (listp (elt ph 3))
                                           (elt ph 3)
                                         (list (elt ph 3)))))))))
     (sli-flatten sli-structures))
    res))

(defun sli-agglomerate (lst)
  "LST is a list of list (beg end).
If beg1 = beg2= ... = begN, we answer (beg1 end1 end2 ... endN)."
  (let ((res '()) beg (listend '()))
    (mapc
     (lambda (ph)
       (unless (assoc (setq beg (elt ph 0)) res) ;; already done
         (setq listend '())
         (mapc 
          (lambda (nph)
            (when (equal (elt nph 0) beg)
              (add-to-list 'listend (elt nph 1))))
          lst)
         (setq res (append res (list (append (list beg) listend))))))
     lst)
    res))

(defun sli-get-special-head-head-alist nil
  (let ((res '()) previous-head (previous-keys '()))
    (mapc
     (lambda (ph)
       (cond
	((equal (elt ph 1) 'head)
	 (setq previous-head (list (elt ph 0)) previous-keys (list (elt ph 0))))
	((and (equal (elt ph 1) 'special-head) (member (elt ph 0) sli-head-keys))
	 (add-to-list 'res (cons (elt ph 0) previous-head)); (print (list (elt ph 0) previous-keys))
	 (add-to-list 'sli-special-head-previous-keys-alist (cons (elt ph 0) previous-keys)))
	(t (add-to-list 'previous-keys (elt ph 0)))))
     (sli-flatten sli-structures))
    ;; Some work for sli-special-head-previous-keys-alist and res:
    ;;   some special-head are linked to different things.
    (setq sli-special-head-previous-keys-alist (sli-agglomerate sli-special-head-previous-keys-alist))
    (sli-agglomerate res)))

(defun sli-get-max-keys-length (lst)
  (let ((res 0))
    (mapc (lambda (to) (setq res (max res to)))
          (mapcar 'length lst))
    res))

(defun sli-precomputations nil
  ;; variables:
  ;(princ "\nPrecomputations: variables")
  (setq sli-head-keys (sli-scan-structures 'head)
        sli-special-head-keys (sli-scan-structures 'special-head)
        sli-soft-keys (sli-scan-structures 'soft)
        sli-beacon-keys (sli-scan-structures 'beacon)
        sli-math-relation-keys (sli-scan-structures 'math-relation)
        sli-relation-keys (append sli-beacon-keys sli-math-relation-keys)
        sli-strong-keys (sli-scan-structures 'strong)
        sli-end-keys (sli-scan-structures 'end)
        sli-constructor-keys (sli-scan-structures 'constructor)
        sli-keys-nomrelations (append sli-head-keys sli-soft-keys sli-strong-keys sli-beacon-keys
				      sli-special-head-keys ;; momentanous !!
                                      sli-end-keys)
	sli-keys (append sli-keys-nomrelations sli-relation-keys)
        sli-max-keys-length (sli-get-max-keys-length sli-keys))
  ;(princ "...done.\n")
  ;;regexps:
  ;(princ "\nPrecomputations: regexps")
  (setq sli-all-end-strong-regexp (sli-regexp-opt (append sli-end-keys sli-strong-keys))
        sli-fixed-regexp (sli-regexp-opt (mapcar 'car sli-fixed-keys-alist))
        sli-head-regexp (sli-regexp-opt sli-head-keys)
        sli-strong-regexp (sli-regexp-opt sli-strong-keys)
        sli-all-keys-nomrelations-noseparators-regexp
          (sli-regexp-opt (append sli-keys-nomrelations sli-comment-starts
                                 (list "\"" block-comment-start block-comment-end)))
        sli-all-keys-regexp
          (sli-regexp-opt (append sli-keys sli-separators sli-comment-starts
                                  (list "\"" block-comment-start block-comment-end)))
        sli-all-keys-and-constructors-regexp
          (sli-regexp-opt (append sli-keys sli-separators sli-comment-starts
                                  sli-constructor-keys
                                  (list "\"" block-comment-start block-comment-end))))
  ;(princ "...done.\n")
  ;; association lists:
  ;(princ "\nPrecomputations: alists")
  (setq sli-ends-head-alist (sli-get-ends-head-alist)
	sli-head-end-alist (sli-get-head-end-alist)
	sli-heads-strong-alist (sli-get-heads-strong-alist) ; sli-ambiguous-keys also is partly created there.
        sli-companion-strong-keys-alist (sli-get-companion-alist)
        sli-soft-alist (sli-get-soft-alist)
        sli-soft-head-or-strong-alist (sli-get-soft-head-or-strong-alist)
	sli-special-head-alist (sli-get-special-head-alist)
	sli-special-head-heads-alist (sli-get-special-head-head-alist) ;; sli-special-head-previous-keys-alist is also created here
                   sli-relevant-alist (sli-get-relevant-alist)
        sli-ancestors-alist (sli-get-ancestors-alist)
	;; offsets :
        sli-first-offset-alist (sli-get-first-offset-alist)
        sli-second-offset-alist (sli-get-second-offset-alist)
        sli-relation-offset-alist (sli-get-relation-offset-alist)
        sli-special-head-offset-alist (sli-get-special-head-offset-alist)
        ;; the maid :
        sli-maid-alist (sli-get-maid-alist) ; sli-ambiguous-keys also is partly created there.
        )
  ;(princ "...done.\n")
  )

;;;--------------------------------------------------------------------------------------
;;; End of the section devoted to precomputations from sli-structures.
;;;--------------------------------------------------------------------------------------

;;;--------------------------------------------------------------------------------------
;;; This section is devoted to some simple functions extracting informations
;;; from the variables defined above.
;;;--------------------------------------------------------------------------------------

  ;; A full-key is a cons (STRING . PT) where PT is the
  ;; value of point at the beginning of STRING.

(defsubst sli-keyword (el)
  (if sli-case-fold (downcase el) el))

(defsubst sli-member (el lst)
  (if sli-case-fold (member (downcase el) lst) (member el lst)))

(defsubst sli-following-key (key)
  (cdr (assoc (sli-keyword key) sli-maid-alist)))

(defun sli-indent-after (key &optional before-soft)
  ;; answer is an integer or a cons ('absolute . integer)
  (setq key (sli-keyword key))
  (eval
   (cond
    ;; See how special-heads are handled: if specified by sli-key-is-a-special-headp
    ;; put to t they take precedence, otherwise the head-offset has precedence.
    ;; If no head exist then the offset as a special-head is finally used.
    (sli-key-is-a-special-headp
     (cdr (assoc key sli-special-head-offset-alist)))
    ((and before-soft (sli-member key (append sli-head-keys sli-strong-keys)))
     (cdr (assoc key sli-first-offset-alist)))
    ((sli-member key (append sli-head-keys sli-strong-keys))
     (cdr (assoc key sli-second-offset-alist)))
    ((sli-member key sli-relation-keys)
     (cdr (assoc key sli-relation-offset-alist)))
    ((sli-member key sli-soft-keys)
     (cdr (assoc key sli-second-offset-alist)))
    ((sli-member key sli-special-head-keys)
     (cdr (assoc key sli-special-head-offset-alist)))
    (t 0))))

(defsubst sli-get-shift (beg end)
  (or (eval (cdr (assoc (vector (sli-keyword beg) (sli-keyword end))
                        sli-shift-alist))) 0))

(defsubst sli-get-strongs-from-strong-or-head (strong)
  (cdr (assoc (sli-keyword strong) sli-companion-strong-keys-alist)))

(defsubst sli-get-heads-from-end (end)
  (cdr (assoc (sli-keyword end) sli-head-end-alist)))

(defsubst sli-get-heads-from-strong (strong)
  (cdr (assoc (sli-keyword strong) sli-heads-strong-alist)))

(defsubst sli-get-ends-from-head (head)
  (cdr (assoc (sli-keyword head) sli-ends-head-alist)))

(defsubst sli-get-head-and-strong-from-soft (soft)
  (cdr (assoc (sli-keyword soft) sli-soft-alist)))

(defsubst sli-get-ends-from-strong (strong)
  (sli-flatten
   (mapcar 'sli-get-ends-from-head
           (sli-get-heads-from-strong strong))))

(defsubst sli-get-relevant (key)
  (cdr (assoc (sli-keyword key) sli-relevant-alist)))

(defsubst sli-get-special-head-previous-keys (key)
  (cdr (assoc (sli-keyword key) sli-special-head-previous-keys-alist)))

(defsubst sli-get-special-head-previous-heads (key)
  (cdr (assoc (sli-keyword key) sli-special-head-heads-alist)))

(defsubst sli-possible-ancestors (key)
  (cdr (assoc (sli-keyword key) sli-ancestors-alist)))

;;;-------------------------------------------------------------------------------------------
;;; Some general primitives.
;;;-------------------------------------------------------------------------------------------

(defsubst sli-remove-trailing-spaces nil
  (if (and (looking-at "\\s-+\\($\\|\\'\\)")
           (not (text-property-any (match-beginning 0) (match-end 0) 'read-only t)))
      (delete-horizontal-space)))

(defsubst sli-remove-trailing-spaces-previous-line nil
  (save-excursion
    (forward-line -1)
    (end-of-line)
    (save-restriction
      (condition-case err
          (unwind-protect
              (progn
                (narrow-to-region (line-beginning-position) (point))
                (while (and (progn
                              (forward-char -1)
                              (looking-at "\\s-"))
                            (not (text-property-any (match-beginning 0) (match-end 0) 'read-only t))))
                (unless (looking-at "\\s-") (forward-char 1)); in case we are not at bol
                (when sli-verbose
                  (princ "\n")
                  (princ (list "(sli-remove-trailing-spaces-previous-line) removing spaces from/to: "
                               (point) (line-end-position))))
                (delete-char (- (line-end-position) (point))))
            (widen))
        (error (when sli-verbose (princ "\n(sli-remove-trailing-spaces-previous-line): ") (princ err)) nil)))))

(defsubst sli-only-spacep (&optional pt)
  ;; t if the line contains only spaces.
  (unless pt (setq pt (point)))
  (let ((only-spacep t))
    (mapc (lambda (ch) (setq only-spacep
                             (and only-spacep (= (char-syntax ch) ?\ ))))
          (string-to-list
           (buffer-substring-no-properties (line-beginning-position) pt)))
    only-spacep))

(defun sli-only-spaces-on-line-before nil
  "t if point is between beginning-of-line
and first non-whitespace character, nil else.
nil if point is at beginning of line."
  (let (res)
    (save-excursion
      (save-restriction
        (narrow-to-region (line-beginning-position) (line-end-position))
        (skip-syntax-forward " ") ; beware: linefeed/newline are whitespaces
        (setq res
              (if (= 0 (current-column))
                  nil
                (= (current-indentation) (current-column)))))
      (widen))
    res))

(defun sli-backward-to-indentation nil
  (interactive)
  (if (not (sli-only-spaces-on-line-before))
      (delete-char -1)
    (let ((foundp nil) (cc (current-indentation)) ncc)
      ;;(if sli-verbose
      ;;  (print (list "(sli-backward-to-indentation) Current indentation: " cc)))
      (save-excursion
        (while (and (not (bobp)) (not foundp))
	  (forward-line -1)
          (beginning-of-line) ; for the bobp to work
          (setq foundp (> cc (setq ncc (current-indentation))))))
      (save-restriction
        (narrow-to-region (line-beginning-position) (line-end-position))
        (skip-syntax-forward " ")
        (if (not foundp)
	    (backward-delete-char-untabify cc)
          (backward-delete-char-untabify (- cc ncc)))
        (widen)))))

(defsubst sli-point-to-indent (pt)
  (save-excursion
    (progn (goto-char pt) (current-column))))

(defsubst sli-indent-at (full-key)  ;; used only here
  ;; A full-key is a cons (STRING . PT) where PT is the
  ;; value of point at the beginning of STRING. PT alone is also accepted.
  (sli-point-to-indent (if (consp full-key) (car full-key) full-key)))

(defsubst sli-in-one-line-comment nil
  (and sli-comment-starts ; if sli-comment-starts is nil, answer is nil
      (re-search-backward (regexp-opt sli-comment-starts) (line-beginning-position) t)))

(defsubst sli-get-safe-backward-place nil
  (save-excursion
    (when (eobp) (forward-char -1))
    (if (re-search-backward sli-safe-place-regexp nil t)
        (match-end 1) (point-min))))

(defsubst sli-get-safe-forward-place nil
  (save-excursion
    (when (bobp) (forward-char 1))
    (if (re-search-forward sli-safe-place-regexp nil t)
        (match-beginning 1) (point-max))))

(defsubst sli-within-long-comment nil
  (let*((aux (sli-get-safe-backward-place))
	(res (parse-partial-sexp aux (point)))) ;(princ (list " Yol " (nth 4 res) (not (nth 7 res))))
    (if (and (nth 4 res) (not (nth 7 res)))
        (nth 8 res)
      nil)))

(defun sli-anchored-posix-search-backward (regexp lim &optional no-error)
;;; ??? DOES NOT SEEM TO WORK:  (posix-search-backward regexp lim no-error))
  (let ((case-fold-search sli-case-fold))
    (and (re-search-backward regexp lim no-error)
         (let*((end-pt (match-end 0))
               (beg (- end-pt sli-max-keys-length)))
           ;(princ "\n") (princ (list "Anchored posix. Candidate: " (match-beginning 0) (match-end 0)  " beg=" beg))
           ;;(princ (save-excursion (goto-char beg) (posix-search-forward regexp end-pt t)))
           (while (save-excursion
                    (goto-char beg)
                    (posix-search-forward regexp end-pt t)
                    (< (match-end 0) end-pt))
             ;;(princ "\n") (princ (list "Inside anchored posix: " (match-beginning 0) " beg=" beg))
             (setq beg (1+ beg)))
           ;(princ "\n") (princ (list "Out of anchored posix: " (match-beginning 0) " beg=" beg))
           (goto-char (match-beginning 0))))))

;;;---------------------------------------------------------------------------------
;;;  Handling text properties
;;;---------------------------------------------------------------------------------

(defsubst sli-prop-should-remove (beg props)
  (let ((lola 0) (res t))
    (or (and (get-text-property beg 'sli-time)
             (> (- (cadr (current-time)) (get-text-property beg 'sli-time))
                sli-prop-do-not-recompute-time))
        (progn
          (while (< lola (/ (length props) 2))
            (setq res (and res (get-text-property beg (elt props (* lola 2))))
                  lola (+ 1 lola)))
          ;; res is nil if one of the properties that PROPS wants to set
          ;; is not already set.
          (not res)))))

(defsubst sli-prop-word (beg)
  (buffer-substring-no-properties beg (next-property-change beg)))

(defsubst sli-prop-full-key (beg)
  (cons (buffer-substring-no-properties beg (next-property-change beg)) beg))

(defsubst sli-prop-region (beg)
  (cons beg (next-property-change beg)))

(defun sli-prop-renew (beg end props)
  "PROPS is '(sli-type head sli-ancestor 66) for instance."
  (let ((old-buff-modp (buffer-modified-p)))
    (when (sli-prop-should-remove beg props)
      (remove-text-properties beg end '(sli-type nil sli-ancestor nil sli-reverse-ancestor nil sli-time nil))
      (when sli-prop-verbose
        (princ "\n((sli-prop-renew) propertying ")(princ beg))
      (add-text-properties beg end props)
      (add-text-properties beg end (list 'sli-time (cadr (current-time))))
      (set-buffer-modified-p old-buff-modp))))

(defsubst sli-prop-renew2 (full-key props)
  "Same as sli-prop-renew except that full-key replaces BEG END"
  (sli-prop-renew (cdr full-key) (+ (cdr full-key) (length (car full-key))) props))

(defsubst sli-prop-has-type (beg)
  "Answer sli-type at BEG if it exists and is not stale.
Answer is nil otherwise."
  (if (sli-prop-should-remove beg '(sli-time 0)) nil
    (setq sli-prop-used (+ 1 sli-prop-used))
    (get-text-property beg 'sli-type)))

;;;---------------------------------------------------------------------------------
;;; The real stuff starts here.
;;;---------------------------------------------------------------------------------
;;;
;;; Functions to get pairs ....
;;;

(defun sli-reduce-skel (skel &optional full)
  ; (cdr skel) is reduced if FULL is nil. With a t value,
  ; (cdr skel goes through reduction.
  (if (null skel) nil
   (let*((word (car skel)) end-lst strong-lst
         (found-strongp nil) (found-endp nil)
         (skel (if full (sli-reduce-skel (cdr skel) t) (cdr skel))))
     (cond
       ((sli-member word sli-end-keys) ; don't do a thing !
        (append (list word) skel))
       ((sli-member word sli-head-keys)
        ;; its end should be below or it is the key we seek. Erase this closed part.
        (setq end-lst (sli-get-ends-from-head word))
        ;(princ "\n") (princ (list "(sli-reduce-skel): end-lst is " end-lst))
        (while (and skel (not (sli-member (car skel) end-lst)))
          (setq skel (cdr skel)))
        ;(princ "\n") (princ (list "(sli-reduce-skel): last skel is " skel))
        (if (null skel) (list word) (cdr skel))) ; the answer.
       ((sli-member word sli-strong-keys)
        ;; its end should be below or it is the key we seek.
        (setq end-lst (sli-get-ends-from-strong word)
              strong-lst (sli-get-strongs-from-strong-or-head word))
        (mapc (lambda (s)
                (setq found-endp (or found-endp (sli-member s end-lst))
                      found-strongp (or found-strongp (sli-member s strong-lst))))
              skel)
        (cond
         (found-endp
          (while (and skel (not (sli-member (car skel) end-lst)))
            (setq skel (cdr skel))))
         ;; So word is a strong key with no end below.
         (found-strongp
          (while (and skel (not (sli-member (car skel) strong-lst)))
            (setq skel (cdr skel)))
          (when (and (cdr skel) (sli-member (cadr skel) strong-lst))
            (setq skel (cdr skel)))))
        (append (list word) skel))))))
 
(defun sli-find-matching-key (pt whatwewant relevant &optional givekey forspecialhead) ; goes backward
"PT is supposedly at beginning of an end/strong-key, out of comment or
string and we look for the first element of WHATWEWANT which is not
in a complete expression. RELEVANT is the list of keys that may
intervene. If GIVEKEY, then full-key is given else key only.
That's a kind of backward-sexp...
If FORSPECIALHEAD is t, then if we find a special-head before PT,
we stop and answer t.
Supports imbedded comments. Answer nil if not found."
  (save-excursion
    (goto-char pt)
    ;(princ "\n") (princ (list "(sli-find-matching-key) getting in with " pt whatwewant relevant))
    (if (and (sli-prop-has-type (point))
             (get-text-property (point) 'sli-ancestor)
             (sli-member (sli-prop-word (get-text-property (point) 'sli-ancestor)) whatwewant))
        (sli-prop-full-key (get-text-property (point) 'sli-ancestor))
      (when (and (sli-prop-has-type (point))
                 (get-text-property (point) 'sli-ancestor))
        ;; but the ancestor is not the good one. Still go till there :
        (setq pt (get-text-property (point) 'sli-ancestor)))
      (let ((level-comment1 0) (skel '())
            (foundp nil) (ans nil) (case-fold-search sli-case-fold)
            word start (in-stringp nil) ancestor
            (aregexp (sli-regexp-opt
                      (append relevant
                              (list "\"" block-comment-start block-comment-end)))))
        (while (and (not foundp) (not (bobp)))
        ;(princ "\n") (princ (list "(sli-find-matching-key) word " word "skel" skel))
          (if (sli-anchored-posix-search-backward aregexp nil 1)
              (cond
               ((string= (setq word (match-string-no-properties 0)) "\"")
                (if (= (preceding-char) ?\\)
                    (setq in-stringp t) ; it should already be.
                  (setq in-stringp (not in-stringp))))
               (in-stringp)
               ; Out of strings:
               ((string= word block-comment-end)
                (sli-prop-renew (match-beginning 0) (match-end 0) '(sli-type block-comment-end))
                (setq level-comment1 (+ 1 level-comment1)))
               ((string= word block-comment-start)
               ; in case the string we look for is a block-comment-start
                (sli-prop-renew (match-beginning 0) (match-end 0) '(sli-type block-comment-start))
                (setq level-comment1 (1- level-comment1))
                (when (and (< level-comment1 0)
                           (equal (list block-comment-start) whatwewant))
                           ; in case the string we look for is a block-comment-start
                           ;(princ (list "Found !" (point)))
                  (setq ans (if givekey (cons word (point)) (point))
                        foundp t)))
               ((sli-member word sli-comment-starts)) ; within a one-line-comment
               ((> level-comment1 0)); within a multiline-comment
               ;; Out of imbedded comments. Now word is in RELEVANT.
               ((not (sli-member word relevant)) ; should not happen!!
                (setq foundp t ans nil))
               ((and forspecialhead
                     (sli-member word whatwewant)) 
                ;; Avoid crossed recursivity of next point.
                (setq foundp t ans (if givekey (cons word (point)) (point))))
               ((setq ancestor (sli-is-a-special-head (point) word))
                ;; crossed recursivity ... But point is going backward !
                (sli-prop-renew2 (cons word (point))
                                (list 'sli-type 'special-head 'sli-ancestor (cdr ancestor)))
                (if (or (sli-separator-directly-afterp pt word)
                        (sli-in-one-line-comment))
                    (goto-char (+ (cdr ancestor) (length (car ancestor))))
                  (setq ans (if givekey (cons word (point)) (point)) foundp t)))
               ((save-excursion (sli-in-one-line-comment)))
               (t (setq skel (sli-reduce-skel (append (list word) skel))
                        forspecialhead nil)
                  (when (and (= 1 (length skel)) (sli-member (car skel) whatwewant))
                    (setq ans (if givekey (cons word (point)) (point))
                          foundp t)))) ; end of cond
            )) ; end of while
            ;(princ "\n") (princ (list "(sli-find-matching-key) out with " ans))
        ans))))

(defsubst sli-special-head-headp (word)
  "Answer not nil if WORD is a special-head that can be a head."
  (assoc (sli-keyword word) sli-special-head-heads-alist))

(defun sli-is-a-special-head (pt word)
  "Answer nil if WORD located at PT is not a special-head.  WORD should not be
in comment, and PT is before WORD.  If WORD is a special-head that can be a
head, answer is nil if it acts like a head; else answer is
(previousword . previouspt) where previousword is the one that showed that word
was a special-head: it is thus a special-head or a head located before (word . pt). "
  (save-match-data
    (cond 
     ((sli-special-head-headp word)
      (cond
        ((and (eq (sli-prop-has-type pt) 'special-head)
              (get-text-property pt 'sli-ancestor))
         (sli-prop-full-key (get-text-property pt 'sli-ancestor)))
        ;; An easy trick: if a separator is not after, it can't be a special-head !
        ((not (sli-separator-directly-afterp (point-max) word))
         (sli-prop-renew pt (+ pt (length word)) '(sli-type head))
         nil)
        (t (let ((appui (sli-find-matching-key
                         pt (sli-get-special-head-previous-heads word)
                         (sli-get-relevant word) t t)))
             (if (consp appui)
                 (sli-prop-renew pt (+ pt (length word))
                                 (list 'sli-type 'special-head 'sli-ancestor (cdr appui)))
               (if appui (sli-prop-renew pt (+ pt (length word)) '(sli-type special-head))
                 (sli-prop-renew pt (+ pt (length word)) '(sli-type head))))
             appui))))
     (t (sli-member word sli-special-head-keys)))))
  
(defun sli-get-corresponding-key (pt whatwewant)
  ;; answer is (block-comment-start . point)
  ;; if PT is within a multiline-comment.
  ;; PT is at the beginning of the word we want to match.
  ;; This function skips
  ;; head/end blocks by using sli-find-matching-key.
  ;; Answers the first element of what we want that is not
  ;; enclosed in a construct.
  (save-excursion
    (goto-char pt)
    (let ((level-comment1 0) (foundp nil) beg aux
          word start (in-stringp nil) (case-fold-search sli-case-fold)
          (relevant (append whatwewant
                            sli-comment-starts
                            (list "\"" block-comment-start block-comment-end)))
          aregexp)
      (dolist (wd whatwewant)
        (dolist (x (cdr (assoc (sli-keyword wd) sli-relevant-alist)))
          (when (sli-member x sli-end-keys) (add-to-list 'relevant x))))
      (setq aregexp (sli-regexp-opt relevant) sli-key-is-a-special-headp nil)
      ;(princ "\n") (princ (list "(sli-get-corresponding-key) getting in " relevant))
      (while (and (not foundp) (not (bobp)))
        (if (sli-anchored-posix-search-backward aregexp nil 1)
          (cond
            ((string= (setq word (match-string-no-properties 0)) "\"")
             (if (= (preceding-char) ?\\)
                 (setq in-stringp t) ; it should already be.
               (setq in-stringp (not in-stringp))))
            (in-stringp)
            ; Out of strings:
            ((string= word block-comment-end)
             (sli-prop-renew (point) (+ (point) (length word)) '(sli-type block-comment-end))
             (setq level-comment1 (1+ level-comment1)))
            ((string= word block-comment-start)
             (sli-prop-renew (point) (+ (point) (length word)) '(sli-type block-comment-start))
             (if (= level-comment1 0)
                 (setq foundp t)
               (setq level-comment1 (1- level-comment1))))
            ((sli-member word sli-comment-starts)) ; within a one-line-comment
            ((> level-comment1 0)); within a multiline-comment
            ;; Out of imbedded comments:
            ((sli-member word sli-end-keys)
             (setq start (point))
             (unless (sli-in-one-line-comment)
               (if (setq beg (sli-find-matching-key
                              start (sli-get-heads-from-end word) (sli-get-relevant word) t))
                   (progn
                     (goto-char (cdr beg))
                     (sli-prop-renew start (+ start (length word))
                                     (list 'sli-type 'end 'sli-ancestor (cdr beg)))
                     (sli-prop-renew2 beg (list 'sli-type 'head 'sli-reverse-ancestor start)))
                 (goto-char (point-min)))))
            ((not (sli-member word whatwewant)))
            ((sli-special-head-headp word) ;; special heads that can be heads
	     (when sli-verbose
	       (princ "\n")
	       (princ
		(list "(sli-get-corresponding-key) Found special-head that could be a head: "
		      word "...")))
	     (if (setq aux (sli-is-a-special-head (point) word))
		 ;; acts like a special head:
		 (unless (or (sli-separator-directly-afterp pt word)
			     (sli-in-one-line-comment))
                   (sli-prop-renew (point) (+ (point) (length word))
                                   (list 'sli-type 'special-head 'sli-ancestor (cdr aux)))
		   (setq foundp t sli-key-is-a-special-headp t))
	       ;; acts like a head:
	       (when sli-verbose (princ "\n(                            ... and is indeed one !)"))
               (sli-prop-renew (point) (+ (point) (length word)) '(sli-type head))
	       (setq foundp (sli-member word whatwewant))))
            ((sli-member word sli-special-head-keys)
             (unless (or (sli-separator-directly-afterp pt word)
                         (sli-in-one-line-comment))
               (setq foundp t)))
            ((sli-member word whatwewant)
             (setq start (point))
             (unless (sli-in-one-line-comment)
                     (setq foundp t)))
            (t nil))
           ))
       ;(princ "\n") (princ (list "(sli-get-corresponding-key) out with " (if foundp (cons word (point)) nil)))
      (if foundp (cons word (point)) nil))))

(defsubst sli-get-key-for-soft (pt soft)
  (sli-get-corresponding-key pt (sli-get-head-and-strong-from-soft soft)))

(defun sli-get-key-for-strong (pt strong)
  (sli-get-corresponding-key pt (sli-get-heads-from-strong strong)))

(defun sli-get-key-for-end (pt end)
  "Looking for head of (END.PT)."
  (sli-get-corresponding-key pt (sli-get-heads-from-end end)))

(defsubst sli-get-head-from-ambiguous (pt key)
  (let (auxkey)
    (cond
     ((sli-member key sli-strong-keys)
      (sli-get-key-for-strong pt key))
     ((sli-member key sli-soft-keys)
      (unless (sli-member (car (setq auxkey (sli-get-key-for-soft (point) key))) sli-head-keys)
        (setq auxkey (sli-get-key-for-strong pt (car auxkey))))
      (if auxkey auxkey 'sli-fail))
     (t 'sli-fail))))

(defun sli-separator-directly-afterp (end word)
  "t if there is SEPARATOR between (1+ point) and end
which is not within a comment or a string and such that
no keyword appear in between  except maybe someone in
sli-constructor-keys."
  (save-excursion
    (forward-char 1)
    ;(princ "\n") (princ (list "Getting in sli-separator-directly-afterp with " (point) end word))
    (let ((level-comment1 0) (level 0) (foundp nil)
           wd (in-stringp nil) (directlyp nil)
           (separators (cdr (assoc (sli-keyword word) sli-special-head-alist)))
           (case-fold-search sli-case-fold))
      (while (and (not foundp) (< (point) end))
        (when (posix-search-forward sli-all-keys-and-constructors-regexp end 1)
          (cond
            ((string= (setq wd (match-string-no-properties 0)) "\"")
             (if (= (preceding-char) ?\\)
                 (setq in-stringp t) ; it should already be.
               (setq in-stringp (not in-stringp))))
            (in-stringp)
            ; Out of strings:
            ((string= wd block-comment-end)
             (setq level-comment1 (1- level-comment1)))
            ((string= wd block-comment-start)
             (setq level-comment1 (1+ level-comment1)))
            ((sli-member wd sli-comment-starts) (forward-line 1)) ; within a one-line-comment
            ((> level-comment1 0)); within a multiline-comment
            ;; Out of imbedded comments:
            ((and (member wd separators) (sli-is-a-separatorp (1- (point))))
             (setq foundp t directlyp t))
            ((sli-member wd sli-constructor-keys))
            (t (setq foundp t))
          )))
      ;(princ "\n") (princ (list "Out of sli-separator-directly-afterp. directlyp =  " directlyp))
      directlyp)))

;;;----------------------------------------------------------------------------
;;;--- beginning of forward/backward/scan-sexp/s
;;;----------------------------------------------------------------------------
 
(defsubst sli-move-a-bit-before nil
  (let ((p (point))(case-fold-search sli-case-fold))
    (save-restriction
      (unwind-protect
          (progn
            (narrow-to-region
             (progn (re-search-backward "\\s-" nil 1)
                    (when (and (not (eobp))
                               (not (member (char-syntax (char-after)) '(?w ?_ ?\( ?\) ?$))))  
                      (forward-char 1)); at beob
                    (point))
             (progn (re-search-forward "\\s-" nil 1) 
                    (when (and (not (bobp))
                               (not (member (char-syntax (preceding-char)) '(?w ?_ ?\( ?\) ?$)))) 
                      (forward-char -1)); at eob
                    (point)))
            (goto-char p)
            (when (member (char-syntax (preceding-char)) '(?w ?_ ?\( ?\) ?$))
              (skip-syntax-backward "w_()$"))
            (while (and (<= (point) p);(princ (list "sli-move-a-bit-before" (point)))
                        (posix-search-forward sli-all-keys-regexp nil t))); we have gone too far.
            (sli-anchored-posix-search-backward sli-all-keys-regexp nil 1))
	(widen)))
    (when sli-verbose (print (list p (point))))
    (if (> (point) p) (progn (goto-char p) nil) t)))

(defun sli-skip-to-beginning-of-keyword nil
  (sli-move-a-bit-before))

(defun sli-find-full-key-at-point (&optional move)
  (save-excursion
    (if (or (sli-move-a-bit-before) move)
        (progn
          ;(princ "\n")(princ (list "sli-find-full-key-at-point"(match-string-no-properties 0) (point)))
          (cons (match-string-no-properties 0) (point)))
      nil)))

(defun sli-backward-sexp (&optional arg)
  "A backward-sexp. If point is after an end or a strong,
go to its head. If point is in the middle of the text,
use backward-word. If ARG, repeat that many times.
Answer POINT of where to go" 
  (save-restriction
    (condition-case err
        (progn
          (if (and arg (< arg 0))
              (sli-forward-sexp (- arg))
            (let ((n (or arg 1)) first-stuff beg pt (modifiedp (buffer-modified-p))
                  (case-fold-search sli-case-fold))
              (while (> n 0)
                (setq first-stuff (sli-find-full-key-at-point t))
                (goto-char (setq pt (cdr first-stuff)))
               (when sli-verbose
                  (princ "\n") (princ (list "(sli-backward-sexp) to be matched: " first-stuff)))
                (if (or (null first-stuff)
                        (search-backward " " pt t)
                        (not (sli-member (car first-stuff) (append sli-end-keys sli-strong-keys))))
                    ;; The previous word is not an end or a strong:
                    (progn
                      (when sli-verbose
                        (princ "\n") (princ (list "(sli-backward-sexp) nothing special")))
                      ;; Do *not* use backward-sexp, it is advised !!!
                      (forward-word -1))
                  (cond
                   ((and (sli-prop-has-type (cdr first-stuff))
                         (get-text-property (cdr first-stuff) 'sli-ancestor)
                         (sli-member (sli-prop-word (get-text-property (cdr first-stuff) 'sli-ancestor))
                                 sli-head-keys))
                    (goto-char (get-text-property (cdr first-stuff) 'sli-ancestor)))
                   ((sli-member (car first-stuff) sli-end-keys)
                    (setq beg (sli-get-key-for-end
                               (if (and (sli-prop-has-type (cdr first-stuff))
                                        (get-text-property (cdr first-stuff) 'sli-ancestor))
                                   ;; An ancestor exists. It is a strong. Still it is better than nothing.
                                   (get-text-property (cdr first-stuff) 'sli-ancestor)
                                 (cdr first-stuff))
                               (car first-stuff)))
                    (when sli-verbose
                      (princ "\n") (princ (list "(sli-backward-sexp) match: " beg)))
                    (cond
                     ((and (consp beg) (equal (car beg) block-comment-start));;un'
                      (sli-prop-renew2 beg '(sli-type block-comment-start))
                      (goto-char (cdr beg)))
                     ((consp beg)
                      (sli-prop-renew2 first-stuff (list 'sli-type 'end 'sli-ancestor (cdr beg)))
                      (sli-prop-renew2  beg (list 'sli-type 'head 'sli-reverse-ancestor (cdr first-stuff)))
                      (goto-char (cdr beg)))
                     (t nil)))
                   ((sli-member (car first-stuff) sli-strong-keys)
                    (setq beg (sli-get-key-for-strong (cdr first-stuff) (car first-stuff)))
                    (when sli-verbose (princ "\n") (princ (list "(sli-backward-sexp) match: " beg)))
                    (cond
                     ((and (consp beg) (equal (car beg) block-comment-start));;un'
                      (sli-prop-renew2  beg '(sli-type block-comment-start))
                      (goto-char (cdr beg)))
                     ((consp beg)
                      (sli-prop-renew2 first-stuff (list 'sli-type 'strong 'sli-ancestor (cdr beg)))
                      (sli-prop-renew2 beg (list 'sli-type 'head 'sli-reverse-ancestor (cdr first-stuff)))
                      (goto-char (cdr beg)))
                     (t nil)))
                   (t (when sli-verbose (princ "\n(sli-backward-sexp) Should not be here!)")))))
                (setq n (- n 1)))
              (set-buffer-modified-p modifiedp)))
          (when sli-verbose (princ "\n") (princ (list "(sli-backward-sexp) answer: " (point))))
          (point))
      (error (princ "\n(sli-backward-sexp): ") (princ err) nil))))

(defun sli-find-end-forward (pt word)
  "WORD is a head or a strong. PT is at beginning of WORD.
Answer is (endword . endpoint)."
  (let ((whatwewant-regexp (if (sli-member word sli-head-keys)
                               (sli-regexp-opt (sli-get-ends-from-head word))
                             (sli-regexp-opt (sli-get-ends-from-strong word))))
        foundp end his-head (case-fold-search sli-case-fold))
    (if (and (sli-prop-has-type pt)
             (get-text-property pt 'sli-reverse-ancestor))
        (sli-prop-full-key (get-text-property pt 'sli-reverse-ancestor))
      ;; Start the swallow/unswallow process :  
      (save-restriction
        (unwind-protect
            (progn 
              (narrow-to-region (sli-get-safe-backward-place) (sli-get-safe-forward-place))
              (while (and (re-search-forward whatwewant-regexp nil t)
                          (not foundp))
                (goto-char (match-beginning 0))
                (setq end (cons (match-string-no-properties 0) (match-beginning 0)))
                (when sli-verbose
                  (princ "\n")
                  (princ (list "(sli-find-end-forward) Potential end:" end)))
                (setq his-head (sli-get-key-for-end (point) (car end))
                      foundp (or (null his-head) ; meaning we don't understand a thing!
                                 (and (consp his-head) (<= (cdr his-head) pt))))
                (when sli-verbose
                  (princ "\n")
                  (princ (list "(sli-find-end-forward) His head:" his-head)))
                (when (consp his-head)
                  (sli-prop-renew2 end (list 'sli-type 'end 'sli-ancestor (cdr his-head)))
                  (sli-prop-renew2
                   his-head (list 'sli-type 'head 'sli-reverse-ancestor (cdr end))))
                ;; In case the end found was closing something in between, continue from after:
                (goto-char (+ (cdr end) (length (car end))))
                ))
          (widen)))
      (if foundp end nil))))

(defun sli-forward-sexp (&optional arg)
  "A forward-sexp. If point is before a head or a strong,
go to its end. If point is in the middle of the text,
use forward-word. If ARG, repeat that many times.
Answer POINT of where to go." 
  (save-restriction
    (condition-case err
        (progn
          (if (and arg (< arg 0))
              (sli-backward-sexp (- arg))
            (let ((n (or arg 1)) end beg aux (modifiedp (buffer-modified-p))
                  (case-fold-search sli-case-fold))
              (while (> n 0)
                (sli-skip-to-beginning-of-keyword)
                (cond
                 ((posix-looking-at (regexp-opt (if (boundp 'block-comment-start)
                                                    (append sli-comment-starts (list block-comment-start))
                                                    sli-comment-starts)))
                  ;; In comment: use text forward-sexp.
                  (when sli-verbose (princ "\n((sli-forward-sexp) comments)"))
                  ;; Do *not* use forward-sexp !!!
                  (forward-word 2))
                 ((or (setq aux (member (sli-prop-has-type (point)) '(head strong)))
                      (and (posix-looking-at sli-head-regexp)
                           (not (sli-is-a-special-head (match-beginning 0) (match-string-no-properties 0))))
                      (posix-looking-at sli-strong-regexp))
                  (if aux
                      (setq beg (sli-prop-full-key (point))
                            end (sli-find-end-forward (point) (car beg)))
                    (setq beg (cons (match-string-no-properties 0) (match-beginning 0))
                          end (sli-find-end-forward (point) (match-string-no-properties 0))))
                  (when sli-verbose
                    (princ "\n") (princ (list "(sli-forward-sexp) to be matched: " beg))
                    (princ "\n") (princ (list "(sli-forward-sexp) match: " end)))
                  (cond
                   ((and (consp end) (equal (car end) block-comment-end));;un'
                    (sli-prop-renew2 beg (list 'sli-type (if (sli-member (car beg) sli-head-keys)
                                                             'head 'strong)))
                    (sli-prop-renew2 end '(sli-type block-comment-end))
                    (goto-char (+ (length block-comment-end) (cdr end))))
                   ((consp end) 
                    (sli-prop-renew2 
                     beg (list 'sli-type (if (sli-member (car beg) sli-head-keys) 'head 'strong)
                               'sli-reverse-ancestor (cdr end)))
                    (sli-prop-renew2 end (list 'sli-type 'end 'sli-ancestor (cdr beg)))
                    (goto-char (+ (length (car end)) (cdr end))))
                   (t nil)))
                 (t (when sli-verbose (princ "\n((sli-forward-sexp) nothing found)"))
                    (forward-word 2)))
                (setq n (- n 1)))
          (set-buffer-modified-p modifiedp)))
          (when sli-verbose (princ "\n") (princ (list "(sli-forward-sexp) answer: " (point))))
          (point))
      (error (princ "\n(sli-forward-sexp): ") (princ err) nil))))

(defun sli-scan-sexps (pt count)
  (goto-char pt)
  (when sli-verbose (princ "\n((sli-scan-sexps))"))
  (if (< count 0)
      (sli-backward-sexp count)
    (sli-forward-sexp count)))
    
(defvar sli-select-end-of-overlay-fn
  'sli-select-end-of-overlay-fn-default
"Function used to give the end of the overlay.
Takes two arguments KEY and PT.
Default value is `sli-select-end-of-overlay-fn-default'.")

(defun sli-select-end-of-overlay-fn-default (key pt)
  (+ pt (length key)))

(defun sli-select-end-of-overlay (key pt)
  (funcall sli-select-end-of-overlay-fn key pt))

(defun sli-show-sexp (&optional arg)
  "POINT is on a head or end key.
This key is highlighted as well as its corresponding end/head.
Color used is `show-paren-match-face'. Nothing is highlighted
if no corresponding key is found. 
  When used with prefix C-u, remove stale text properties and
recompute things by setting `sli-prop-do-not-recompute-time' to 0."
  (interactive "P")
  (save-excursion
    (save-restriction
      (let ((old-sli-prop-do-not-recompute-time sli-prop-do-not-recompute-time))
        (unwind-protect
            (let ((full-key (sli-find-full-key-at-point)) pt
                  (modifiedp (buffer-modified-p)))
              (when sli-verbose
                (princ "\n")
                (princ (list "(sli-show-sexp) full-key-at-point: " full-key)))
              (when full-key (setq pt (goto-char (cdr full-key))))
              (when (and arg (= (car arg) 4)) ;; call prefixed by C-u
                (setq sli-prop-do-not-recompute-time 0))
              (setq sli-prop-used 0)
              (narrow-to-region (sli-get-safe-backward-place) (sli-get-safe-forward-place))
              (cond
               ((and full-key (sli-member (car full-key) sli-head-keys)
                     (not (sli-is-a-special-head (cdr full-key) (car full-key))));(print full-key)
                     (move-overlay sli-overlay-beg (cdr full-key) 
                                   (sli-select-end-of-overlay (car full-key) (cdr full-key)))
                     (if (and (sli-forward-sexp)
                              (equal (get-text-property (1- (point)) 'sli-type) 'end)); ?(1- point) ??
                         (progn 
                           (overlay-put sli-overlay-beg 'face 'show-paren-match-face)
                           (overlay-put sli-overlay-end 'face 'show-paren-match-face)
                           (when sli-prop-verbose
                             (princ "\n")
                             (princ (list "(sli-show-sexp) overlay-end:"
                                          (get-text-property pt 'sli-reverse-ancestor) (point))))
                           (goto-char (1- (point)))
                           (setq full-key (sli-find-full-key-at-point));(print full-key)
                           (move-overlay sli-overlay-end (get-text-property pt 'sli-reverse-ancestor)
                                         (sli-select-end-of-overlay (car full-key) (cdr full-key))))
                       (overlay-put sli-overlay-beg 'face 'show-paren-mismatch-face)
                       (move-overlay sli-overlay-end (point-min) (point-min))))
               ((and full-key (sli-member (car full-key) sli-end-keys))
                (move-overlay sli-overlay-end (cdr full-key)
                              (sli-select-end-of-overlay (car full-key) (cdr full-key)))
                (if (and (sli-backward-sexp)
                         (equal (get-text-property (point) 'sli-type) 'head))
                    (progn 
                      (overlay-put sli-overlay-beg 'face 'show-paren-match-face)
                      (overlay-put sli-overlay-end 'face 'show-paren-match-face)
                      (when sli-prop-verbose
                        (princ "\n")
                        (princ (list "(sli-show-sexp) overlay-beg:"
                                     (get-text-property pt 'sli-ancestor) (next-property-change (point))))
                        (princ "\n")
                        (princ (list "(sli-show-sexp) number of text-properties used:" sli-prop-used)))
                      (goto-char (get-text-property pt 'sli-ancestor))
                      (setq full-key (sli-find-full-key-at-point))
                      (move-overlay sli-overlay-beg (cdr full-key)
                                    (sli-select-end-of-overlay (car full-key) (cdr full-key))))
                  (overlay-put sli-overlay-end 'face 'show-paren-mismatch-face)
                  (move-overlay sli-overlay-beg (point-min) (point-min))))
               (t ;; Erase overlays:
                (when sli-prop-verbose
                  (princ (list "\n(sli-show-sexp) Erasing overlays")))
                (move-overlay sli-overlay-beg (point-min) (point-min))
                (move-overlay sli-overlay-end (point-min) (point-min))))
              (set-buffer-modified-p modifiedp)
              (widen)
              (setq sli-prop-do-not-recompute-time old-sli-prop-do-not-recompute-time)))))))

(defvar sli-show-sexp-idle-timer nil)

(defun sli-show-sexp-semi-mode (arg)
  "When ARG>0 corresponding head/end keys are automatically
shown with an idle timer. When ARG=0, sli-show-sexp is bound
to f8. When ARG is anything else, remove `sli-overlay-beg' and
`sli-overlay-end'."
  (when sli-show-sexp-idle-timer
      (cancel-timer sli-show-sexp-idle-timer))
  (cond
   ((< 0 arg)
    (setq sli-show-sexp-idle-timer
          (run-with-idle-timer (if (featurep 'lisp-float-type) (/ (float 1) (float 8)) 1)
                               t 'sli-show-sexp)))
   ((= 0 arg)
    (move-overlay sli-overlay-beg (point-min) (point-min))
    (move-overlay sli-overlay-end (point-min) (point-min))
    (local-set-key [f8] 'sli-show-sexp))
   (t 
    (move-overlay sli-overlay-beg (point-min) (point-min))
    (move-overlay sli-overlay-end (point-min) (point-min)))))

(defadvice forward-sexp (around sli-handles-forward-sexp (&optional arg))
  (interactive)
  (if (bound-and-true-p sli-handles-sexp) (sli-forward-sexp arg) ad-do-it))

(defadvice backward-sexp (around sli-handles-backward-sexp (&optional arg))
  (interactive)
  (if (bound-and-true-p sli-handles-sexp) (sli-backward-sexp arg) ad-do-it))

(require 'advice)
(ad-activate 'forward-sexp  'around)
(ad-activate 'backward-sexp 'around)

;;;----------------------------------------------------------------------------
;;;--- end of forward/backward/scan-sexp/s
;;;----------------------------------------------------------------------------
;;;
;;; Indentation
;;;

(defun sli-get-first-fixed-or-strong-or-end-or-soft (pt)
  ; Go to first non whitespace char on line on which PT lies and before PT.
  ; Then nil if within comment or first word is not a fixed/end/strong/soft key,
  ; the cons (KEY . point-at-its-beginning) otherwise.
  (save-excursion
    (save-restriction
      (unwind-protect
	  (let (aux (case-fold-search sli-case-fold))
	    (narrow-to-region (progn (beginning-of-line) (point)) pt)
	    (skip-chars-forward " \t")
            ;(princ "\n") (princ (list "(sli-get-first-fixed-or-strong-or-end-or-soft)" (point)))
	    (cond ((setq aux (sli-prop-has-type (point)))
                   (cond ((member aux '(block-comment-end block-comment-start))
                          (cons (eval aux) (point)))
                         ((and (or (member aux '(end strong soft))
                                   (assoc (sli-keyword (sli-prop-word (point))) sli-fixed-keys-alist))
                               (<= (next-property-change (point)) pt))
                          (sli-prop-full-key (point)))
                         (t nil)))
                  ((posix-looking-at (regexp-opt (append sli-comment-starts (list block-comment-start))))
                   (sli-prop-renew (match-beginning 0) (match-end 0) (list 'sli-type 'block-comment-start))
		   (cons block-comment-start (point)))
		  ((posix-looking-at (regexp-opt (list block-comment-end)))
		   (sli-prop-renew (match-beginning 0) (match-end 0) (list 'sli-type 'block-comment-end))
		   (cons block-comment-end (point)))
                  ((posix-looking-at (sli-regexp-opt sli-soft-keys))
                   (sli-prop-renew (match-beginning 0) (match-end 0) (list 'sli-type 'soft))
                   (cons (match-string-no-properties 0) (point)))
		  ((or (posix-looking-at sli-fixed-regexp)
		       (posix-looking-at sli-all-end-strong-regexp))
		   (cons (match-string-no-properties 0) (point)))
		  (t nil)))
	(widen)))))

(defun sli-get-first-non-end-key (pt &optional nomrelation) ; goes backward
"Find first non-end-key before PT outside comment
or string which is not matched by an end-key.
Imbedded comments are supported.
If NOMRELATION is t, then this key is not a math-relation
either. Answer is a full-key (KEY . POINT)
where POINT indicates the beginning of the occurence
of KEY we're interested in.
Answer is (block-comment-start . point)
if PT is within a multiline-comment."
  (save-excursion
    (goto-char pt)
    (let ((level-comment1 0) (foundp nil) beg
          (accessible-separator (sli-member (char-to-string (preceding-char)) sli-separators))
          word start (in-stringp nil) (case-fold-search sli-case-fold)
          (aregexp
             (if nomrelation sli-all-keys-nomrelations-noseparators-regexp sli-all-keys-regexp)))
      (setq sli-key-is-a-special-headp nil)
      (while (and (not foundp) (not (bobp)))
        (if (sli-anchored-posix-search-backward aregexp nil 1)
          (progn ;(princ "\n") 
        	 ;(princ (list "(sli-get-first-non-end-key). word = " (match-string-no-properties 0) (point)))
            (cond
             ((string= (setq word (match-string-no-properties 0)) "\"")
              (if (= (preceding-char) ?\\)
                  (setq in-stringp t) ; it should already be.
                (setq in-stringp (not in-stringp))))
             (in-stringp)
             ;; Out of strings:
             ((string= word block-comment-end)
              (setq start (point))
              ;(princ "\n") (princ (list "(sli-get-first-non-end-key) In block-comment."))
              (unless (sli-in-one-line-comment)
                (if (setq beg (sli-find-matching-key
                               start (list block-comment-start) (list block-comment-start) t))
                    (progn 
                      (goto-char (cdr beg))
                      (sli-prop-renew start (+ start (length word))
                                      (list 'sli-type 'block-comment-end 'sli-ancestor (cdr beg)))
                      (sli-prop-renew2
                        beg (list 'sli-type 'block-comment-start 'sli-reverse-ancestor start)))
                  (setq level-comment1 (1+ level-comment1))
                  (goto-char  (point-min)))))
             ((string= word block-comment-start)
              (sli-prop-renew start (+ start (length word)) '(sli-type block-comment-start))
              (if (= level-comment1 0)
                  (setq foundp t)
                (setq level-comment1 (1- level-comment1))))
             ((sli-member word sli-comment-starts)) ; within a one-line-comment
             ((> level-comment1 0)); within a multiline-comment
             ;; Out of imbedded comments:
             ((sli-is-a-separatorp) ; only if NOMRELATION is t.
              (setq start (point))
              (unless (sli-in-one-line-comment)
                (goto-char start) (setq accessible-separator t)))
             ((sli-member word sli-math-relation-keys) ; only if NOMRELATION is t.
              (unless accessible-separator
                (setq start (point))
                (unless (sli-in-one-line-comment)
                  (goto-char start) (setq foundp t))))
             ((sli-member word sli-end-keys)
              (setq start (point))
              (unless (sli-in-one-line-comment)
                (if (setq beg (sli-find-matching-key 
                               start (sli-get-heads-from-end word) (sli-get-relevant word) t))
                    (progn
                      (goto-char (cdr beg))
                      (sli-prop-renew start (+ start (length word))
                                      (list 'sli-type 'end 'sli-ancestor (cdr beg)))
                      (sli-prop-renew2 beg (list 'sli-type 'head 'sli-reverse-ancestor start)))
                  (goto-char (point-min)))))
             ((sli-special-head-headp word) ;; special heads that can be heads
              (when sli-verbose
                (princ "\n")
                (princ
                 (list "(sli-get-first-non-end-key) Found a special head that could be a head: "
                       word " at " (point) "...")))
              (if (sli-is-a-special-head (point) word)
                  ;; acts like a special head:
                  (unless (or (sli-separator-directly-afterp pt word)
                              (sli-in-one-line-comment))
                    (sli-prop-renew (point) (+ (point) (length word)) '(sli-type special-head))
                    (setq foundp t sli-key-is-a-special-headp t))
                ;; acts like a head:
                (when sli-verbose (princ "\n((sli-get-first-non-end-key) ... and is indeed one !)"))
                (sli-prop-renew (point) (+ (point) (length word)) '(sli-type head))
                (setq foundp t)))
             ((sli-member word sli-special-head-keys);(princ " lyo ")
              (unless (or (sli-separator-directly-afterp pt word)
                          (sli-in-one-line-comment))
                (setq foundp t)))
             ((sli-member word sli-separators))      ;; momentanous
             (t (setq foundp (not (sli-in-one-line-comment))))))
          ))
          ;(princ "\n")
          ;(princ (list "Out of sli-get-first-non-end-key with "
          ;		   (if foundp (cons word (point)) nil) accessible-separator))
      (if foundp (cons word (point)) nil))))


(defsubst sli-compute-indent-after (full-key &optional before-soft)
  (let ((the-indent (sli-indent-after (car full-key) before-soft))) ;(princ full-key)
    ;(princ (list "Yummy!!" the-indent)) 
    (throw 'indent (if (consp the-indent)
		       (cdr the-indent) ; absolute indent
		     (+ (sli-point-to-indent (cdr full-key))
			the-indent)))))

(defsubst sli-on-same-linep (pt1 pt2)
  ;(princ "\n") (princ (list "(sli-on-same-linep)" pt1 pt2 ?\n
  ;                          (string-to-list (buffer-substring-no-properties pt1 pt2))))
  (if (member ?\n (string-to-list (buffer-substring-no-properties pt1 pt2)))
      nil t))

(defun sli-tell-indent-within-long-comment (afterp pos-beg-comment)
  (when sli-verbose
    (princ "\n")
    (princ (list "(sli-tell-indent-within-long-comment) getting in with afterp = " afterp
                 " and pos-beg-comment = "pos-beg-comment)))
  ;; AFTERP like in sli-tell-indent.
  ;; If pos-beg-comment and (point) are on the same line, do nothing:
  (when (and (not afterp) (sli-on-same-linep pos-beg-comment (point)))
    (when sli-verbose
      (princ "\n")
      (princ "((sli-tell-indent-within-long-comment) On same line as beginning of comment : no indent.)"))
    (throw 'indent 0))
  (let*((pos-first-char (save-excursion
                          (goto-char (+ pos-beg-comment (length block-comment-start)))
                          (skip-syntax-forward "^w") (point)))
        (on-same-linep (and (or (not afterp) (< pos-first-char (point)))
                            ; because if afterp is true, a \n will be inserted just before (point)
                            (sli-on-same-linep pos-beg-comment pos-first-char)))
        (pos-end-comment (save-excursion
                           (goto-char (+ pos-beg-comment (length block-comment-start)))
                           (search-forward block-comment-end nil t)))
        (end (line-end-position))
        (special-last-linep (and pos-end-comment
                                 (= pos-end-comment
                                    (save-excursion
                                      (beginning-of-line);(princ (point))(princ " ")
                                      (skip-syntax-forward "-" end);(princ (point))(princ " ")
                                      (skip-syntax-forward "^-" end);(princ (point))(princ " ")
                                      (point))))))
  
  ;; check whether heredity should apply:
  ;(princ (count-lines pos-beg-comment (point)))
  (when (and (not afterp)
             (not special-last-linep)
             (> (count-lines pos-beg-comment (point)) 2))
    (throw 'indent (save-excursion
                     (forward-line -1)
                     (current-indentation))))
  ;; Else align on the start :
  (when sli-verbose
    (princ "\n")
    (princ (list "(sli-tell-indent-within-long-comment) align on first line?"
                 (and on-same-linep (not special-last-linep)))))
  (if (and on-same-linep (not special-last-linep))
      (throw 'indent (sli-point-to-indent pos-first-char))
    ;; Special treatment of last line of comment:
    (when sli-verbose
      (princ "\n")
      (princ (list "(sli-tell-indent-within-long-comment) last line?" special-last-linep)))
    (if special-last-linep
        ;; only one word on this line ending with block-comment-end.
        ;; For instance "**/"
        (throw 'indent (+ (sli-point-to-indent pos-beg-comment)
                          (length block-comment-start)
                          sli-block-comment-end-offset
                          ))
      (throw 'indent (+ (sli-point-to-indent pos-beg-comment)
                        (length block-comment-start)
                        sli-block-comment-middle-offset
                        ))))))

(defun sli-tell-indent (&optional afterp nomrelation point-is-the-end) ;; used only here
  "Gives the indentation of line on which point lies.
Or on line after if AFTERP is t."
  ;; This indentation depends on what is on the previous
  ;; line except that the first word of the line could be
  ;; a strong or end key in which case it is to be aligned
  ;; on the previous head/strong of the same block.
  ;; The only thing we don't do is if a string spreads across lines.
  (sli-remove-trailing-spaces); for current-indentation
  (catch 'indent
  (let ((pos-beg-comment (if afterp (sli-within-long-comment)
                           (save-excursion
                             (beginning-of-line)
                             (sli-within-long-comment)))))
    (when pos-beg-comment
      (when sli-verbose
        (princ "\n") (princ (list "(sli-tell-indent) looking on next line ?" afterp))
        (princ "\n")
        (princ (list "(sli-indent-line) Within long comment starting at " pos-beg-comment)))
      (sli-tell-indent-within-long-comment afterp pos-beg-comment)))
  
  (unless (or afterp point-is-the-end) (end-of-line))
  
  (let*((pt (point)) wd-lst beg-str full-key appui head opp the-indent
        (first-stuff (and (not afterp) (sli-get-first-fixed-or-strong-or-end-or-soft pt)))
        is-a-fixed-keyp)
    (when sli-verbose
      (princ "\n") (princ (list "(sli-tell-indent) looking on next line ?" afterp))
      (princ "\n") (princ (list "(sli-tell-indent) first-stuff on line = " first-stuff)))
    ; Zeroth case, indentation of this line and (car first-stuff) is a block-comment-end:
    (when (and (not (null first-stuff))
               (string= (car first-stuff) block-comment-end))
      (when sli-verbose
        (princ "\n") (princ (list "(sli-tell-indent) first-stuff is block-comment-end")))
      (throw 'indent 0))
    ; First case, indentation of this line and (car first-stuff) is a fixed key:
    (when (and (not (null first-stuff))
               (setq opp (assoc (sli-keyword (car first-stuff)) sli-fixed-keys-alist)))
      (when sli-verbose
	(princ "\n") (princ (list "(sli-tell-indent) first-stuff is in sli-fixed-keys-alist")))
      (setq is-a-fixed-keyp t)
      ;; Old treatment:
      ;(throw 'indent (+ (save-excursion (forward-line -1) (current-indentation))
      ;                  (eval (cdr opp)))))
      )
    ; Second case, line starts by a soft key:
    ; it has to be done in case of "if 2<3 \n then" since the "then"
    ; has been aligned with respect to the math-relation and not to the "if"
    (when (and first-stuff (sli-member (car first-stuff) sli-soft-keys))
      (setq appui (sli-get-key-for-soft (cdr first-stuff) (car first-stuff)))
      (when sli-verbose
	(princ "\n") (princ (list "(sli-tell-indent) first-stuff is in sli-soft-keys")))
      (sli-compute-indent-after appui))
    ; Third case, indentation of this line
    ; and (car first-stuff) is not a fixed key or a comment or a soft-key:
    (when (and first-stuff (not (string= (car first-stuff) block-comment-start))
               (not is-a-fixed-keyp))
      ; line starts by a strong/end key. We select the key from which to
      ; compute the indent. Usually we align it on the previous head/strong
      ; key and add possible offset. That's the heredity principle. But we can also
      ; align strong/end-keys on the head if this head is in sli-no-heredity-list.
      ; Another case is when the previous corresponding strong/head had the
      ; attribute 'absolute, in which case its indentation applies.
      (setq appui
            (sli-find-matching-key   ; backward
             (cdr first-stuff) ; where to start the search.
             (sli-possible-ancestors (car first-stuff))
             (sli-get-relevant (car first-stuff)) t))
      ; see whether the absolute attribute is present:
      (when (and (not (null appui))
		 (consp (setq the-indent (sli-indent-after (car appui))))
		 (eq (car the-indent) 'absolute))
        (sli-prop-renew2 first-stuff
                         (list 'sli-type (if (sli-member (car first-stuff) sli-end-keys) 'end 'strong)
                               'sli-ancestor (cdr appui)))
        (sli-prop-renew2 appui
                         (list 'sli-type (if (sli-member (car appui) sli-head-keys) 'head 'strong)
                               'sli-reverse-ancestor (cdr first-stuff)))
        (when sli-verbose
          (princ "\n") (princ (list "(sli-indent) Absolute indent. Indent resting on: " (car appui))))
	(throw 'indent (+ (cdr the-indent)
			  (sli-get-shift (car appui) (car first-stuff)))))
      ; see whether heredity applies:
      (unless (or (null appui) (sli-member (car appui) sli-head-keys))
        ; select head from appui and not from full-key because
        ; (1) it is shorter (2) (car head) *is* a strong key.
        (setq head (sli-get-head-from-ambiguous (cdr appui) (car appui)))
        ;(princ "\n") (princ (list "heredity ? for " (vector (car head) (car first-stuff))))
        (sli-prop-renew2 first-stuff
                         (list 'sli-type (if (sli-member (car first-stuff) sli-end-keys) 'end 'strong)
                               'sli-ancestor (cdr appui)))
        (if (eq head 'sli-fail)
            (sli-prop-renew2 appui (list 'sli-type 'strong 'sli-reverse-ancestor (cdr first-stuff)))
          (sli-prop-renew2 appui
                           (list 'sli-type 'strong
                                 'sli-reverse-ancestor (cdr first-stuff) 'sli-ancestor (cdr head)))
          (sli-prop-renew2 head (list 'sli-type 'head 'sli-reverse-ancestor (cdr appui)))
          (when (sli-member (vector (car head) (car first-stuff)) sli-no-heredity-list)
            (setq appui head))))
      (when sli-verbose
	(princ "\n((sli-tell-indent) indentation of this line and not in comment)")
	(princ "\n") (princ (list "                  Resting on: " (car appui) (cdr appui))))
      (throw 'indent (if (null appui) 0
                       (+ (sli-get-shift (car appui) (car first-stuff))
                          (sli-indent-at (cdr appui))))))
    ; Fourth case, indentation of this line and (car first-stuff) is a comment:
    (when (and first-stuff (string= (car first-stuff) block-comment-start))
      ; PT is within multi-line-comment.
      (sli-prop-renew2 first-stuff '(sli-type block-comment-start))
      (when sli-verbose
	(princ "\n((sli-tell-indent) indentation of this line and in comment)"))
      (throw 'indent (current-indentation)))

    (unless afterp
      ; ; Fifth case : line doesn't start by a strong/end/soft key:
      (save-excursion
        (if (= -1 (forward-line -1))
            ; we are already on the first line:
            (if first-stuff (throw 'indent (current-indentation))
                (throw 'indent 0)))
	(when sli-verbose
	  (princ "\n((sli-tell-indent) line doesn't start by a strong/end/soft key)"))
        (end-of-line)
        (setq pt (point))))

    ;; This point can be reached only if AFTERP is t OR first-stuff is nothing special
    ;; (which could be a fixed key).
    (setq first-stuff (sli-get-first-non-end-key pt nomrelation)) ; backward search
    ;; sli-key-is-a-special-headp is set.
    (when sli-verbose
      (princ "\n") (princ (list "(sli-tell-indent) indentation of line after?" afterp))
      (princ "\n") (princ (list "(sli-tell-indent) key deciding of indent = " first-stuff)))
       
    (cond
      ((null first-stuff)
       ;; no construct active or within comment. Don't do a thing:
       (when sli-verbose
	  (princ "\n((sli-tell-indent) no construct active or within comment)"))
       (throw 'indent (current-indentation)))
      ((string= (car first-stuff) block-comment-start)
       (sli-prop-renew2 first-stuff '(sli-type block-comment-start))
       (when sli-verbose
	  (princ "\n") (princ (list "(sli-tell-indent) within comment")))
       (throw 'indent (current-indentation)))
      (sli-key-is-a-special-headp ;; a special head;
       (when sli-verbose (princ "\n((sli-tell-indent) within special-head.)"))
       (sli-compute-indent-after first-stuff t))
      ((and (sli-member (car first-stuff) (append sli-head-keys sli-strong-keys))
            (not (assoc (sli-keyword (car first-stuff)) sli-soft-head-or-strong-alist)))
       ;; head/strong without soft:
       (when sli-verbose
	  (princ "\n")
          (princ (list "(sli-tell-indent) within a head/strong construct never followed by a soft")))
       (sli-prop-renew2 first-stuff (list 'sli-type (if (sli-member (car first-stuff) sli-head-keys) 
                                                        'head 'strong)))
       (sli-compute-indent-after first-stuff))
      ((sli-member (car first-stuff)
                   (append sli-head-keys sli-strong-keys sli-special-head-keys))
       ;; head/strong with soft missing or special-head:
       (sli-prop-renew2 first-stuff (list 'sli-type 
                                          (cond ((sli-member (car first-stuff) sli-head-keys) 'head)
                                                ((sli-member (car first-stuff) sli-strong-keys) 'strong)
                                                (t 'special-head))))
       (when sli-verbose
	  (princ "\n((sli-tell-indent) within special-head or head/strong sometimes")
          (princ "\n                   followed by currently missing soft)"))
       (sli-compute-indent-after first-stuff t))
      ((and is-a-fixed-keyp
            (sli-member (car first-stuff) sli-relation-keys))
       (throw 'indent
              (+ (eval (cdr opp))
                 (save-excursion
                   (goto-char (cdr first-stuff))
                   (beginning-of-line)
                   (skip-syntax-forward "-" (cdr first-stuff))
                   (sli-point-to-indent (point))))))
      ((sli-member (car first-stuff) sli-relation-keys)
       ; relation: if it is just before point ignore it:
       ;  (but can you tell me why????)
       (if (save-excursion
	     (save-restriction
	       (unwind-protect
		   (progn
		     (narrow-to-region (goto-char (cdr first-stuff)) pt)
		     (posix-looking-at (concat (car first-stuff) " *$")))
		 (widen))))
           (save-excursion
             (goto-char (cdr first-stuff))
             (sli-tell-indent t t point-is-the-end))
	 (when sli-verbose
	   (princ "\n") (princ (list "(sli-tell-indent) last non-end-key is in sli-relation-keys")))
	 (sli-compute-indent-after first-stuff)))
      ((sli-member (car first-stuff) sli-soft-keys)
       ; a soft key. Find its head/strong and align things on it.
       (setq full-key (sli-get-key-for-soft (cdr first-stuff) (car first-stuff)))
       (when sli-verbose
	 (princ "\n") (princ (list "(sli-tell-indent) last non-end-key is in sli-soft-keys")))
       (sli-compute-indent-after full-key))))))

;;;-----------------------------------------------------------------------
;;;  Functions that are used outside. Avoid using the two first ones
;;;  as they are not nicely surrounded by a condition-case !
;;;-----------------------------------------------------------------------

(defsubst sli-safe-insert (wd)
  (unless (get-text-property (point) 'read-only)
    (insert wd)))

(defsubst sli-insert-indent (ind)
  (or (null ind)
    (let ((beg (point)) (last (current-column)) move-p (cc (current-indentation))
          (old-buff-modp (buffer-modified-p)))
      (when sli-verbose
	(princ "\n") (princ (list "(sli-insert-indent) indent for: " (point))))
      ;;(princ "\n") (princ (list "(sli-insert-indent) buffer-modifiedp: " old-buff-modp))
      (save-excursion
        (setq move-p (re-search-backward "[^ \t]" (line-beginning-position) t))
        (beginning-of-line)
        (if (get-text-property (point) 'read-only)
            (setq move-p t)
          (delete-horizontal-space) ; Simply because I Hate \t chars.
          (indent-to ind)) ;(insert-char ?  ind)
        )                  ;(princ "\nInserting indent: done.")
      ;; If ind is cc on unmodified buffer, declare the buffer as unmodified:
      (set-buffer-modified-p (or old-buff-modp (not (= cc ind))))
      ;; if point was inside the removed spaces,
      ;; then now it is at the beginning of the line.
      ;; Not what we wanted.
      ;(princ "\n") (princ (list "Deplacement Automatique ?" move-p))
      (unless move-p ; point has been moved automatically
        (move-to-column ind))
      )))

(defun sli-indent-line nil
  (save-restriction
    (condition-case err
        (save-excursion
          (sli-insert-indent (sli-tell-indent)))
      (error (princ "\n(sli-indent-line): ") (princ err) nil))))

(defun sli-indent-region (beg end)
  (interactive "r")
  (save-restriction
    (condition-case err
        (save-excursion
          (setq end (progn (goto-char end) (end-of-line) (point)))
          (narrow-to-region (progn (goto-char beg) (sli-get-safe-backward-place))
                            (progn (goto-char end) (sli-get-safe-forward-place)))
          (when sli-verbose
            (princ "\n")
            (princ (list "(sli-indent-region) Narrowing to: " (point-min) (point-max))))
          ;; Use text-properties as much as possible:
          (let ((sli-prop-do-not-recompute-time 10000) (modifiedp (buffer-modified-p)))
            (remove-text-properties beg end '(sli-type nil))
            (goto-char beg)
            (while (progn (sli-indent-line)
                          (and (re-search-forward "$" end t)
                               (not (= end (point)))))
              (forward-line 1))
            (set-buffer-modified-p modifiedp)))
      (error (princ "\n(sli-indent-region): ") (princ err) nil))))

(defun sli-electric-tab nil ;; linked to 'indent-line-function
  "The interactive counterpart of 'sli-indent-line.
Does a number of other things:
 -- if there are nothing but spaces between beginning-of-line
    and (point), then indents the line and sends (point)
    to the first non space ot tab character of the line.
 -- else if sli-tab-always-indent then indents the line
    the cursor being 'relatively' fixed.
In a program, use `sli-indent-line'."
  (interactive)
  (save-restriction
    (condition-case err
        (unwind-protect
            (progn
              (setq sli-prop-used 0)
              (narrow-to-region (sli-get-safe-backward-place) (sli-get-safe-forward-place))
              (if (sli-only-spacep)
                  (progn
                    (sli-indent-line)
                    (skip-chars-forward " \t"))
                (when sli-tab-always-indent (sli-indent-line)))
              (when sli-verbose
                (princ "\n")
                (princ (list "(sli-electric-tab) number of text-properties used:" sli-prop-used))))
          (widen))
      (error (princ "\n(sli-electric-tab): ") (princ err) nil))))

(defun sli-electric-terminate-line (&optional beg)
  "Terminate line and indent next line."
  (interactive)
  (save-restriction
    (condition-case err
        (unwind-protect
            ;(if (sli-within-long-comment)
            ;    (sli-put-newline)
              (setq sli-prop-used 0)
	      (when sli-verbose
		(princ "\n")
		(princ (list "(sli-electric-terminate-line) narrowing to "
			     (sli-get-safe-backward-place) (sli-get-safe-forward-place))))
              (narrow-to-region (sli-get-safe-backward-place) (sli-get-safe-forward-place))
              (let (this-indent next-indent only-spacep)
                (sli-remove-trailing-spaces)
                (setq only-spacep (sli-only-spacep))
				; (princ "\n") (princ (list "only-spacep = " only-spacep))
                (sli-insert-indent (setq this-indent (sli-tell-indent nil nil t)))
                (unless only-spacep (sli-safe-insert " "))
                                        ;--> in case of thendo with point between then and do.
                (setq next-indent (sli-tell-indent t nil t))
                (when sli-verbose
		  (princ "\n") (princ (list "(sli-electric-terminate-line) indent before:" this-indent))
		  (princ "\n") (princ (list "(sli-electric-terminate-line) indent after:" next-indent)))
                (unless only-spacep (if (= (char-syntax (preceding-char)) ?\ )(delete-char -1)))
		;(princ "\n") (princ (list "(sli-electric-terminate-line) inserting a newline at: " (point)))
                (sli-put-newline)
                (sli-remove-trailing-spaces-previous-line)
		;(princ "\n") (princ (list "(sli-electric-terminate-line) inserting indent at: " (point)))
                (sli-insert-indent next-indent))
              (when sli-verbose
                (princ "\n")
                (princ (list "(sli-electric-terminate-line) number of text-properties USED:" sli-prop-used)))
          (widen))
      (error (princ "\n(sli-electric-terminate-line): ") (princ err) nil))))

(defun sli-newline (&optional beg)
  "Insert a newline without indenting current line.
Next line is properly indented."
  (interactive)
  (save-restriction
    (condition-case err
        (unwind-protect
            ;(if (sli-within-long-comment)
            ;    (sli-put-newline)
              (narrow-to-region (sli-get-safe-backward-place) (sli-get-safe-forward-place))
              (sli-remove-trailing-spaces)
              (sli-put-newline)
              (sli-remove-trailing-spaces-previous-line)
              (sli-insert-indent (sli-tell-indent nil nil t))
          (widen))
      (error (princ "\n(sli-newline): ") (princ err) nil))))

(defun sli-maid (&optional arg on-listp)
  "Closes constructs for you, puts the children to bed and
may order a pizza if you know how to ask.
 Usually, adds the corresponding part of `sli-add-to-key-alist'
except when the call is prefixed by C-u. If the variable
`sli-more-maidp' is nil, this behaviour is reversed.
The word to pursue the structure is taken from `sli-maid-alist'.
This list is created automatically but can be corrected
by specifying special furtherings in `sli-maid-correction-alist'"
  (interactive "P")
  (save-restriction
    (condition-case err
        (unwind-protect
            ;; *Before* any narrowing, check the possibility of inserting !!
            (unless (get-text-property (point) 'read-only)
              (narrow-to-region (sli-get-safe-backward-place) (sli-get-safe-forward-place))
              (let*((full-key (sli-get-first-non-end-key (point) t)) (key nil) (head nil) smore
                    (where-to-write '()) is-a-special-head-head-keyp has-answered)
		(when sli-verbose
		  (princ "\n")
		  (princ (list "(sli-maid) Key to be continued: " full-key)))
                (sli-remove-trailing-spaces)
                ;; Sort ambiguity arising from ambiguous-keys:
                (when (and full-key (sli-member (car full-key) sli-ambiguous-keys))
                  (if (eq (setq smore (sli-get-head-from-ambiguous (cdr full-key) (car full-key))) 'sli-fail)
                      (setq head 'sli-fail)
                    (setq head (car smore)))
		  (when sli-verbose
		    (princ "\n")
		    (princ (list "(sli-maid) The previous key was soft/strong and ambiguous. Its head is : " head))))
                ;; Sort ambiguity head-special-head-keys:
                (when (and full-key (sli-special-head-headp (car full-key)))
		  ;(print (list "yes" (cdr full-key) (sli-get-special-head-previous-keys (car full-key))
		  ;             (sli-get-relevant (car full-key))))
                  (setq head
			(sli-find-matching-key
			 (cdr full-key) 
                         (sli-get-special-head-previous-keys (car full-key)) (sli-get-relevant (car full-key)) t))
		  (setq is-a-special-head-head-keyp t)
		  (when sli-verbose
		    (princ "\n")
		    (princ (list "(sli-maid) The previous key was a head-special-head. Its head is : " head))))
                ;; Go out of one-line-comment:
                (when (save-excursion (sli-in-one-line-comment))
                  (if on-listp
                      (setq where-to-write (append where-to-write (list 'newline)))
                    (sli-electric-terminate-line)))
                ;; add a newline before insertion if required:
                (unless (sli-only-spacep)
                  (when (or (and (not is-a-special-head-head-keyp)
                                 full-key (sli-member (car full-key) sli-keys-with-newline))
			    (and is-a-special-head-head-keyp (not head)
				 (sli-member head sli-keys-with-newline)))
                    (if on-listp
                        (setq where-to-write (append where-to-write (list 'newline)))
                      (sli-electric-terminate-line))))
                ;(princ "\n") (princ (list "Inside mupad-maid. full-key/head = " full-key head))
                ;; find or insert closing-key:
                (cond
                 ((eq head 'sli-fail) (message "Could not resolve ambiguity"))
                 ((null full-key)
                     ;; No construct to be closed.
                  (setq key (buffer-substring-no-properties
                             (save-excursion (forward-word -1) (point)) (point))))
                 ((equal (car full-key) block-comment-start)
                  (if on-listp
                      (setq where-to-write (append where-to-write (list (setq key block-comment-end))))
                    (setq has-answered t)
                    (sli-safe-insert (setq key block-comment-end))))
                 ((and (sli-member (car full-key) sli-separators)
                                        ; Beware !! this key could be **very far**
                       (= (count-lines (cdr full-key) (point)) 0))
                  (setq key nil)) ; We shall put a newline, see below.
		 (is-a-special-head-head-keyp ; a special head possibly a head
		  (if head 
		      ;; it is a special head:
		      (setq key (cadr (assoc (sli-keyword (car full-key)) sli-special-head-alist)))
		    ;; it is a  head:
		    (setq key (sli-following-key (car full-key))))
		  ;(print (list "yes" key head))
		  (unless (and (not (null key))
                                 (or (not (member (char-syntax (string-to-char key)) '(?w ?_ ?\( ?\) ?$)))
                                     (= (char-syntax (preceding-char)) ?\ )))
                      (if on-listp
                          (setq where-to-write (append where-to-write '(" ")))
                        (setq has-answered t)
                        (sli-safe-insert " ")))
		  (if on-listp
                      (setq where-to-write (append where-to-write (list key)))
                    (setq has-answered t)
                    (sli-safe-insert key)))
                 ((and (sli-member (car full-key) sli-special-head-keys)
                       (not (sli-separator-directly-afterp (point-max) (car full-key))))
                  (if on-listp
                      (setq where-to-write
                            (append where-to-write
                                    (list (cadr (assoc (sli-keyword (car full-key)) sli-special-head-alist)))))
                    (setq has-answered t)
                    (sli-safe-insert (cadr (assoc (sli-keyword (car full-key)) sli-special-head-alist)))))
                 (t (setq key (if head ;  completion of an ambiguous-key:
                                  (car (sli-get-ends-from-head head))
                                (sli-following-key (car full-key))))
                               ;(princ " Yol ")         ; add a space if required:
                    (unless (and (not (null key))
                                 (or (not (member (char-syntax (string-to-char key)) '(?w ?_ ?\( ?\) ?$)))
                                     (= (char-syntax (preceding-char)) ?\ )))
                      (if on-listp
                          (setq where-to-write (append where-to-write '(" ")))
                        (setq has-answered t)
                       (sli-safe-insert " ")))
                    (or (null key)
                        (if on-listp
                            (setq where-to-write (append where-to-write (list key)))
                          (setq has-answered t)
                          (sli-safe-insert key)))))
                ;(princ "\n") (princ (list "Inside mupad-maid. key = " key))
                ;; add things if required:
                (unless (if sli-more-maidp
                            (and arg (= (car arg) 4)) ; call is  prefixed by C-u
                          (not (and arg (= (car arg) 4)))) ; call is not prefixed by C-u
                  (cond
                   ((or (null key) (eq head 'sli-fail)))
                   ((setq smore (assoc (sli-keyword key) sli-add-to-key-alist))
                    (if on-listp
                        (setq where-to-write (append where-to-write (list (cdr smore))))
                      (setq has-answered t)
                      (sli-safe-insert (cdr smore))))))
                ;; Add a newline if required:
		;(princ "\n(sli-maid) looking if a newline is required")
                (cond
                 ((or (sli-member key sli-keys-without-newline) (eq head 'sli-fail)))
                 ((eobp) (if on-listp
                             (setq where-to-write (append where-to-write (list 'newline)))
                           (sli-electric-terminate-line)))
                 ((or (null key)
                      (< 2 (count-lines (point)
                                        (save-excursion (skip-chars-forward " \t\n") (point)))))
                  (if on-listp
                      (setq where-to-write (append where-to-write (list 'indent 'forward-line 'indent)))
		    ;(princ "\n(sli-maid) indentation plus going to next line")
                    (sli-indent-line) (forward-line 1) (indent-to (sli-tell-indent nil nil t))))
                                        ; beware if it is only an empty line.
                 (t (if on-listp
                        (setq where-to-write (append where-to-write (list 'indent)))
		      ;(princ "\n(sli-maid) indentation but not going to next line")
                      (sli-indent-line))))
                (unless has-answered (message "Nothing to do"))
                where-to-write))
          (widen))
        (error (princ "\nsli-maid can't understand what to do: ")(princ err) nil))))

(defun sli-tutor nil
  "*Adds what all you should add to end your construct."
  ;; Not so good if used in the middle of a mess ...
  ;; in mupad, try "while foo do" with point before "do".
  (interactive)
  (condition-case err
      (let ((some-more '()) what-to-do)
        (while (and (setq some-more (sli-maid nil t))
                    (not (member some-more '((indent) (newline)))))
                    ; (princ "\n") (princ (list "Tutor:" some-more (point)))
          (while some-more
            (cond
             ((equal (setq what-to-do (car some-more)) 'newline)
              (sli-electric-terminate-line))
             ((equal what-to-do 'indent)
              (sli-indent-line))
             ((equal what-to-do 'forward-line)
              (forward-line 1))
             (t (sli-safe-insert what-to-do)))
            (setq some-more (cdr some-more)))))
    (error (princ "\nsli-tutor can't understand what to do: ")(princ err) nil)))

(defun sli-tools
  (struct shift sep sepp fixed safe keyn keynn mkey comm noher
          &optional newl corr showsexpp case-fold eoov)
"Once these tools are loaded, you should have
`sli-newline' and `sli-electric-terminate-line'
which behave like `newline-and-indent' and
`reindent-then-newline-and-indent'. Also
`indent-line-function' is `sli-electric-tab'
and
`indent-region-function' is `sli-indent-region'.
Finally `sli-backward-to-indentation' is a good
function to bind [backspace] to.

When `sli-handles-sexp' is t then forward-sexp,
backward-sexp and scan-sexps are advised so that
for instance C-M-f on a head sends cursor on its end.

`sli-show-sexp' works like show-paren-mode. Two
ways: either showsexpp is t, either showsexpp is nil
in which case one should press [f8] to see the
corresponding key.
C-u[f8] forces to recompute text-properties locally.

C-M-f/C-M-b run forward-sexp/backward-sexp in a special
way: heads will be atuned to ends and strongs to either
one.
Finally, `sli-maid' tries to further constructs for you
while `sli-tutor' strives to end all constructs.

For these tools to work, the parameters are
`sli-structures'
`sli-shift-alist'
`sli-separators'
`sli-is-a-separatorp-fn'
`sli-fixed-keys-alist'
`sli-safe-place-regexp' ; safe place starts at the end of first grouping
`sli-keys-with-newline'
`sli-keys-without-newline'
`sli-add-to-key-alist'
`sli-comment-starts'
`sli-no-heredity-list'
`sli-put-newline-fn'
`sli-maid-correction-alist'
showsexpp
`sli-case-fold'
`sli-select-end-of-overlay-fn'
and you should also set
`block-comment-start'      `block-comment-end'
`sli-more-maidp'           `sli-tab-always-indent'
and the syntax table should be ok.
Beware that `block-comment-start' and `block-comment-end'
are NOT regexp but simple strings."
  (interactive)
  (condition-case err
      (progn
        (setq sli-structures struct           sli-shift-alist shift
              sli-separators sep              sli-fixed-keys-alist fixed
              sli-case-fold case-fold         sli-keys-with-newline keyn
              sli-keys-without-newline keynn  sli-add-to-key-alist mkey
              sli-comment-starts comm         sli-no-heredity-list noher
              sli-maid-correction-alist corr)
        (if  safe
            (setq sli-safe-place-regexp safe)
          (setq sli-safe-place-regexp "\\(\\'\\|\\`\\)")) ;beginning/end of buffer !
        (when sepp
          (setq sli-is-a-separatorp-fn sepp))
        (when newl
          (setq sli-put-newline-fn newl))
        (when eoov
          (setq sli-select-end-of-overlay-fn eoov))
        (set (make-local-variable 'indent-line-function) 'sli-electric-tab)
        (set (make-local-variable 'indent-region-function) 'sli-indent-region)
        (setq sli-handles-sexp t sli-verbose nil sli-prop-verbose nil)
        (setq sli-overlay-beg (make-overlay (point-min) (point-min)))
        (setq sli-overlay-end (make-overlay (point-min) (point-min)))
        (overlay-put sli-overlay-beg 'face 'show-paren-match-face)
        (overlay-put sli-overlay-end 'face 'show-paren-match-face)
        (overlay-put sli-overlay-beg 'priority 0)
        (overlay-put sli-overlay-end 'priority 0)
        (sli-show-sexp-semi-mode (if showsexpp 1 0))
        (sli-precomputations))
    (error (princ "\nSomething went wrong in sli-tools: ")(princ err) nil)))

;;------------------ sli-tools ends here. 2671 lines ??
