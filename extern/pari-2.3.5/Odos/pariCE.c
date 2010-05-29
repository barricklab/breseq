#include <stdlib.h>
#include <windows.h>


int isspace(int c) {
	if ((c >= 0x09) && (c <= 0x0D)) return 1;
	else if (c == 0x20) return 1;
	return 0;
}

int isdigit(int c) {
	if ((c >= '0') && (c <= '9')) return 1;
	else return 0;
}

int isxdigit(int c) {
	if ((c >= '0') && (c <= '9')) return 1;
	else if ((c >= 'a') && (c <= 'f')) return 1;
	else if ((c >= 'A') && (c <= 'F')) return 1;
	else return 0;
}

int isalpha(int c) {
	if ((c >= 'a') && (c <= 'z')) return 1;
	else if ((c >= 'A') && (c <= 'Z')) return 1;
	else return 0;
}

int isalnum(int c) {
	if ((c >= '0') && (c <= '9')) return 1;
	else if ((c >= 'a') && (c <= 'z')) return 1;
	else if ((c >= 'A') && (c <= 'Z')) return 1;
	else return 0;
}

static int checkDigit(char ch, int base) {
	int n;

	if (ch >= '0' && ch <= '9') n = ch - '0';
	else if (ch >= 'a' && ch <= 'z') n = ch + 10 - 'a';
	else if (ch >= 'A' && ch <= 'Z') n = ch + 10 - 'A';
	else return -1;

	return n >= base ? -1 : n;
}

long strtol(const char *s, char **endptr, int base) {
	const char *sc;
	char sign;
	int actualbase;
	double result;

	/* if the base is illegal, then return 0; */
	if (! (base >= 2 && base <= 36) ){
		if (endptr) *endptr = (char *)s;
		return 0;
	}

	/* skip leading spaces */
 	for (sc = s; *sc == ' '; ++sc);

	/* parse the sign */
	if (*sc == '-' || *sc == '+') {
		sign = *sc;
		sc++;
	}
	else {
		sign = '+';
	}

	/* the default base = 10 */
	actualbase = base == 0 ? 10 : base;

	/* if base is undefined, and number starts '0x', then we have base 16 */
	if (base == 0 && sc[0] == '0' && (sc[1] == 'x' || sc[1] == 'X')) { actualbase = 16; sc += 2; }

	/* else if base is undefined, and number starts '0', then we have base 8 */
	else if (base == 0 && sc[0] == '0') actualbase = 8;

	/* else if base == 16, then skip any leading '0x' */
	else if (base == 16 && sc[0] == '0' && (sc[1] == 'x' || sc[1] == 'X')) sc += 2;

	/* skip leading zeroes */
	for (; *sc == '0'; ++sc);

	/* the result so far is 0. We are going to work with doubles because these give 52 bits of accuracy */
	result = 0.0;

	/* sc now points to the first unprocessed digit. Keep processing until first non digit or overflow */
	for (;;) {
		int d = checkDigit(*sc, actualbase);

		/* if the digit was illegal, then we have terminated */
		if (d < 0) {
			if (endptr) *endptr = (char *)sc;
			return sign == '-' ? -(long)result : (long)result;
		}

		/* roll in the new digit */
		result = result * actualbase + d;

		/* check for overflow */
		if (sign == '+' && result > (double)LONG_MAX) {
			if (endptr) *endptr = (char *)(sc+1);
			return LONG_MAX;
		}
		if (sign == '-' && result > -(double)LONG_MIN) {
			if (endptr) *endptr = (char *)(sc+1);
			return LONG_MIN;
		}

		/* go on to the next character */
		sc++;
	}
}

/*
 *  Implement rename and unlink in win32 procedures. This is complicated
 *  by the fact that these procedures take wide chars.
 */

#define MAX_FILENAME_LEN           128

int rename(const char *oldname, const char *newname) {
	int succ;
	short soldname[MAX_FILENAME_LEN];
	short snewname[MAX_FILENAME_LEN];

	MultiByteToWideChar(CP_ACP, 0, oldname, strlen(oldname)+1, soldname, MAX_FILENAME_LEN);
	MultiByteToWideChar(CP_ACP, 0, newname, strlen(newname)+1, snewname, MAX_FILENAME_LEN);
	succ = MoveFile(soldname, snewname);
	return !succ;
}

int unlink(char *name) {
	int succ;
	short wname[MAX_FILENAME_LEN];

	MultiByteToWideChar(CP_ACP, 0, name, strlen(name)+1, wname, MAX_FILENAME_LEN);

	succ = DeleteFile(wname);
	return !succ;
}

void *calloc(size_t nelem, size_t size) {
	const size_t n = nelem*size;
	char *p = malloc(n);

	if (p) memset(p, '\0', n);
	return p;
}
