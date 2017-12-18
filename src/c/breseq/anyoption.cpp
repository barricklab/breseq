/*
 * AnyOption 1.3
 *
 * kishan at hackorama dot com  www.hackorama.com JULY 2001
 *
 * + Acts as a common facade class for reading
 *   commandline options as well as options from
 *   an optionfile with delimited type value pairs
 *
 * + Handles the POSIX style single character options ( -w )
 *   as well as the newer GNU long options ( --width )
 *
 * + The option file assumes the traditional format of
 *   first character based comment lines and type value
 *   pairs with a delimiter , and flags which are not pairs
 *
 *  	# this is a coment
 *  	# next line is an option value pair
 *  	width : 100
 *     	# next line is a flag
 *      noimages
 *
 * + Supports printing out Help and Usage
 *
 * + Why not just use getopt() ?
 *
 *   getopt() Its a POSIX standard not part of ANSI-C.
 *   So it may not be available on platforms like Windows.
 *
 * + Why it is so long ?
 *
 *   The actual code which does command line parsing
 *   and option file parsing are done in  few methods.
 *   Most of the extra code are for providing a flexible
 *   common public interface to both a resourcefile and
 *   and command line supporting POSIX style and
 *   GNU long option as well as mixing of both.
 *
 * + Please see "anyoption.h" for public method descriptions
 *
 */

/* Updated Auguest 2004
 * Fix from  Michael D Peters (mpeters at sandia.gov)
 * to remove static local variables, allowing multiple instantiations
 * of the reader (for using multiple configuration files).  There is
 * an error in the destructor when using multiple instances, so you
 * cannot delete your objects (it will crash), but not calling the
 * destructor only introduces a small memory leak, so I
 * have not bothered tracking it down.
 *
 * Also updated to use modern C++ style headers, rather than
 * depricated iostream.h (it was causing my compiler problems)
*/

/*
 * Updated September 2006
 * Fix from Boyan Asenov for a bug in mixing up option indexes
 * leading to exception when mixing different options types
 */

#include <sys/ioctl.h>

#include "libbreseq/common.h"

#include "libbreseq/anyoption.h"

namespace breseq {

	AnyOption::AnyOption()
	{
		init();
	}

	AnyOption::AnyOption(const string& usage_prefix)
	{
		init();
		setVerbose(); /* print warnings about unknown options */
		autoUsagePrint(true); /* print usage for bad options */
		addUsage(usage_prefix);
	}

	AnyOption::~AnyOption()
	{
	}

	void
	AnyOption::init()
	{
		argc 		= 0;
		argv 		= NULL;
		posix_style	= true;
		verbose 	= false;
		filename 	= NULL;
		appname 	= NULL;
		new_argv 	= NULL;
		new_argc 	= 0 ;
		max_legal_args 	= 0 ;
		command_set 	= false;
		file_set 	= false;
		g_value_counter = 0;
		command_set 	= false;
		file_set	= false;
		opt_prefix_char     = '-';
		file_delimiter_char = ':';
		file_comment_char   = '#';
		equalsign 	= '=';
		comment       = '#' ;
		delimiter     = ':' ;
		endofline     = '\n';
		whitespace    = ' ' ;
		nullterminate = '\0';
		once = true;
		hasoptions = false;
		autousage = false;

		strcpy( long_opt_prefix , "--" );

		struct winsize ws;
    ioctl(0,TIOCGWINSZ,&ws);
		terminal_width = (ws.ws_col==0) ? 80 : ws.ws_col;
	}

	void
	AnyOption::setCommandPrefixChar( char _prefix )
	{
		opt_prefix_char = _prefix;
	}

	void
	AnyOption::setCommandLongPrefix( char *_prefix )
	{
		if( strlen( _prefix ) > MAX_LONG_PREFIX_LENGTH ){
			*( _prefix + MAX_LONG_PREFIX_LENGTH ) = '\0';
		}

		strcpy (long_opt_prefix,  _prefix);
	}

	void
	AnyOption::setFileCommentChar( char _comment )
	{
		file_delimiter_char = _comment;
	}


	void
	AnyOption::setFileDelimiterChar( char _delimiter )
	{
		file_comment_char = _delimiter ;
	}

	bool
	AnyOption::CommandSet()
	{
		return( command_set );
	}

	bool
	AnyOption::FileSet()
	{
		return( file_set );
	}

	void
	AnyOption::noPOSIX()
	{
		posix_style = false;
	}

	bool
	AnyOption::POSIX()
	{
		return posix_style;
	}


	void
	AnyOption::setVerbose()
	{
		verbose = true ;
	}

	void
	AnyOption::printVerbose()
	{
		if( verbose )
			cout << endl  ;
	}
	void
	AnyOption::printVerbose( const char *msg )
	{
		if( verbose )
			cout << msg  ;
	}

	void
	AnyOption::printVerbose( char *msg )
	{
		if( verbose )
			cout << msg  ;
	}

	void
	AnyOption::printVerbose( char ch )
	{
		if( verbose )
			cout << ch ;
	}

	bool
	AnyOption::hasOptions()
	{
		return hasoptions;
	}

	void
	AnyOption::autoUsagePrint(bool _autousage)
	{
		autousage = _autousage;
	}

	void
	AnyOption::useCommandArgs( int _argc, char** _argv )
	{
		argc = _argc;
		argv = _argv;
		command_set = true;
		appname = argv[0];
		if(argc > 1) hasoptions = true;
	}

	void
	AnyOption::useFilename( const char *_filename )
	{
		filename = _filename;
		file_set = true;
	}

	/*
	 * set methods for options
	 */

	void
	AnyOption::setCommandOption( const char *opt )
	{
		addOption( opt , '\0', COMMAND_OPT );
	}

	void
	AnyOption::setCommandOption( char opt )
	{
		addOption( "", opt , COMMAND_OPT );
	}

	void
	AnyOption::setCommandOption( const char *opt , char optchar )
	{
		addOption( opt , optchar, COMMAND_OPT );
	}

	void
	AnyOption::setCommandFlag( const char *opt )
	{
		addOption( opt , '\0', COMMAND_FLAG );
	}

	void
	AnyOption::setCommandFlag( char opt )
	{
		addOption( "", opt , COMMAND_FLAG );
	}

	void
	AnyOption::setCommandFlag( const char *opt , char optchar )
	{
		addOption( opt , optchar, COMMAND_FLAG );
	}

	void
	AnyOption::setFileOption( const char *opt )
	{
		addOption( opt, '\0', FILE_OPT );
	}

	void
	AnyOption::setFileOption( char opt )
	{
		addOption( "", opt , FILE_OPT );
	}

	void
	AnyOption::setFileOption( const char *opt , char optchar )
	{
		addOption( opt , optchar, FILE_OPT );
	}

	void
	AnyOption::setFileFlag( const char *opt )
	{
		addOption( opt , '\0', FILE_FLAG );
	}

	void
	AnyOption::setFileFlag( char opt )
	{
		addOption( "", opt , FILE_FLAG );
	}

	void
	AnyOption::setFileFlag( const char *opt , char optchar )
	{
		addOption( opt , optchar, FILE_FLAG );
	}

	void
	AnyOption::setOption( const char *opt )
	{
		addOption( opt , '\0', COMMON_OPT );
	}

	void
	AnyOption::setOption( char opt )
	{
		addOption( "", opt , COMMON_OPT );
	}

	void
	AnyOption::setOption( const char *opt , char optchar )
	{
		addOption( opt , optchar, COMMON_OPT );
	}

	void
	AnyOption::setFlag( const char *opt )
	{
		addOption( opt , 0, COMMON_FLAG );
	}

	void
	AnyOption::setFlag( const char opt )
	{
		addOption( "", opt, COMMON_FLAG );
	}

	void
	AnyOption::setFlag( const char *opt , char optchar )
	{
		addOption( opt , optchar, COMMON_FLAG );
	}
  
	void
	AnyOption::addOption( const string opt, char optchar, int type )
	{
    if (opt.size() > 0) {
      ASSERT(!assigned_options.count(opt), "Duplicate option string declared: --" + opt);
      assigned_options.insert(opt);
    }
      
    if (optchar) {
      string optchar_string(1, optchar);
      ASSERT(!assigned_options.count(optchar_string), "Duplicate option char declared: -" + optchar_string);
      assigned_options.insert(optchar_string);
    }
    
		options.push_back(opt);
    optionchars.push_back(optchar);
		optiontype.push_back(type);
		optionindex.push_back(options.size() - 1);
        

	}

	void
	AnyOption::processCommandArgs(int max_args)
	{
		max_legal_args = max_args;
		processCommandArgs();
	}

	void
	AnyOption::processCommandArgs( int _argc, char **_argv, int max_args )
	{
		max_legal_args = max_args;
		processCommandArgs(  _argc, _argv );
	}

	void
	AnyOption::processCommandArgs( int _argc, char **_argv )
	{
		useCommandArgs( _argc, _argv );
		processCommandArgs();
	}

	void
	AnyOption::processCommandArgs()
	{
		if (!CommandSet())
			return;

		if (max_legal_args == 0)
			max_legal_args = argc;

		new_argv = (int*) malloc((max_legal_args + 1) * sizeof(int));

		/* ignore first argv */
		for (int i = 1; i < argc; i++)
		{
			/* long GNU option */
			if (argv[i][0] == long_opt_prefix[0] && argv[i][1] == long_opt_prefix[1])
			{
				int match_at = parseGNU(argv[i] + 2); /* skip -- */
        if (match_at >= 0 && i < argc - 1) /* found match */ {
          //cout << argv[i] << endl;
          //cout << argv[i+1] << endl;
					setValue(options[match_at], argv[++i]);
        }
			}
			/* POSIX char */
			else if (argv[i][0] == opt_prefix_char)
			{
				if (POSIX())
				{
					char ch = parsePOSIX(argv[i] + 1);/* skip - */
          if (ch != '0' && i < argc - 1) /* matching char */ {
            string combined_arg = argv[++i];
            if (i<argc && strchr("\"'", combined_arg[0])) {
              combined_arg.erase(0,1);
              while (i<argc && !strchr("\"'", combined_arg[combined_arg.length()-1])) {
                combined_arg += " ";
                combined_arg + argv[++i];
              }
              combined_arg.erase(combined_arg.length()-1,1);
            }
						setValue(ch, combined_arg);
          }
				}
				/* treat it as GNU option with a - */
				else
				{
					int match_at = parseGNU(argv[i] + 1); /* skip - */
          if (match_at >= 0 && i < argc - 1) /* found match */ {
            string combined_arg = argv[++i];
            if (i<argc && strchr("\"'", combined_arg[0])) {
              combined_arg.erase(0,1);
              while (i<argc && !strchr("\"'", combined_arg[combined_arg.length()-1])) {
                combined_arg += " ";
                combined_arg + argv[++i];
              }
              combined_arg.erase(combined_arg.length()-1,1);
            }
						setValue(options[match_at], combined_arg);
          }
				}
			}
			/* not option but an argument keep index */
			else
			{
				if (new_argc < max_legal_args)
				{
					new_argv[new_argc] = i;
					new_argc++;
				}
				/* ignore extra arguments */
				else
				{
					printVerbose("Ignoring extra argument: ");
					printVerbose(argv[i]);
					printVerbose();
					printAutoUsage();
				}
			}
		}
	}

	char
	AnyOption::parsePOSIX(char* arg)
	{
		for (unsigned int i = 0; i < strlen(arg); i++)
		{
			char ch = arg[i];
			if (matchChar(ch))
			{ /* keep matching flags till an option */
				/*if last char argv[++i] is the value */
				if (i == strlen(arg) - 1)
				{
					return ch;
				}
				else
				{/* else the rest of arg is the value */
					i++; /* skip any '=' and ' ' */
					while (arg[i] == whitespace || arg[i] == equalsign)
						i++;
					setValue(ch, arg + i);
					return '0';
				}
			}
		}
		return '0';
	}

	int
	AnyOption::parseGNU(char *arg)
	{
		int split_at = 0;
		/* if has a '=' sign get value */
		for (unsigned int i = 0; i < strlen(arg); i++)
		{
			if (arg[i] == equalsign)
			{
				split_at = i; /* store index */
				i = strlen(arg); /* get out of loop */
			}
		}
		if (split_at > 0)
		{ /* it is an option value pair */
			char* tmp = (char*) malloc((split_at + 1) * sizeof(char));
			for (int i = 0; i < split_at; i++)
				tmp[i] = arg[i];
			tmp[split_at] = '\0';

			if (matchOpt(tmp) >= 0)
			{
				setValue(options[matchOpt(tmp)], arg + split_at + 1);
				free(tmp);
			}
			else
			{
				printAutoUsage();
        printVerbose("Unknown command argument option: ");
        printVerbose(arg);
        printVerbose();
        printVerbose();
				free(tmp);
        exit(-1);
				//return -1;
			}
		}
		else
		{ /* regular options with no '=' sign  */
			return matchOpt(arg);
		}
		return -1;
	}


	int
	AnyOption::matchOpt( char *opt )
	{
		for( uint32_t i = 0 ; i < options.size() ; i++ ){
			if( options[i].compare(opt) == 0 ){
				if( optiontype[i] ==  COMMON_OPT ||
					optiontype[i] ==  COMMAND_OPT )
				{ /* found option return index */
					return i;
				}else if( optiontype[i] == COMMON_FLAG ||
					   optiontype[i] == COMMAND_FLAG )
				{ /* found flag, set it */
					setFlagOn( opt );
					return -1;
				}
			}
		}
    printAutoUsage();
		printVerbose( "Unknown command argument option: " );
		printVerbose( opt  ) ;
    printVerbose( );
    printVerbose( );
    exit(-1);
		//return  -1;
	}
	bool
	AnyOption::matchChar( char c )
	{
		for( uint32_t i = 0 ; i < optionchars.size() ; i++ ){
			if( optionchars[i] == c ) { /* found match */
				if(optiontype[i] == COMMON_OPT ||
					 optiontype[i] == COMMAND_OPT )
				{ /* an option store and stop scanning */
					return true;
				}else if( optiontype[i] == COMMON_FLAG ||
					  optiontype[i] == COMMAND_FLAG ) { /* a flag store and keep scanning */
					setFlagOn( c );
					return false;
				}
			}
		}
    printAutoUsage();
		printVerbose( "Unknown command argument option: " );
		printVerbose( c ) ;
    printVerbose( );
    printVerbose( );
    exit(-1);
		//return false;
	}

	/*
	 * public get methods
	 */
	string*
	AnyOption::getValue( const string option )
	{
		for( uint32_t i = 0 ; i < options.size() ; i++ )
			if( options[i].compare(option) == 0 )
				return (values.count(optionindex[i]) > 0) ? &values[ optionindex[i] ] : NULL;

		return NULL;
	}

	bool
	AnyOption::getFlag( const string option )
	{
		for( uint32_t i = 0 ; i < options.size() ; i++ )
			if( options[i].compare(option) == 0 )
				return findFlag( optionindex[i] );

		return false;
	}
  
  bool
	AnyOption::isOptUsed( const string option )
	{
		for( int i = 0 ; i < argc ; i++ )
			if( substitute(argv[i], "-", "") == option )
				return true;
    
		return false;
	}

	string*
	AnyOption::getValue( char option )
	{
		for( uint32_t i = 0 ; i < optionchars.size() ; i++ )
			if( optionchars[i] == option )
				return (values.count(optionindex[i]) > 0) ? &values[ optionindex[i] ] : NULL;

		return NULL;
	}

	bool
	AnyOption::getFlag( char option )
	{
		for( uint32_t i = 0 ; i < optionchars.size() ; i++ )
			if( optionchars[i] == option )
				return findFlag( optionindex[i] );

		return false;
	}
  
  bool
	AnyOption::isOptUsed( char option )
	{
		for( int i = 0 ; i < argc ; i++ )
			if( substitute(argv[i], "-", "") == &option )
				return true;
    
		return false;
	}

	bool
	AnyOption::findFlag( int index )
	{
		return (values.count(index) && values[index].compare(TRUE_FLAG) == 0);
	}

	/*
	 * private set methods
	 */
	bool
	AnyOption::setValue( const string option , string value )
	{
		for( uint32_t i = 0 ; i < options.size() ; i++ ){
			if( options[i].compare(option) == 0 ){
        if (values[ optionindex[i] ].size() > 0) values[ optionindex[i] ]+= "\n";
				values[ optionindex[i] ] += value;
				return true;
			}
		}
		return false;
	}

	bool
	AnyOption::setFlagOn( const string option )
	{
		return setValue(option, TRUE_FLAG);
	}

	bool
	AnyOption::setValue( char option , string value )
	{
		for( uint32_t i = 0 ; i < options.size() ; i++ ){
			if( optionchars[i] == option ){
        if (values[ optionindex[i] ].size() > 0) values[ optionindex[i] ]+= "\n";
				values[ optionindex[i] ] += value;
				return true;
			}
		}
		return false;
	}

	bool
	AnyOption::setFlagOn( char option )
	{
		return setValue(option, TRUE_FLAG);
	}


	int
	AnyOption::getArgc( )
	{
		return new_argc;
	}

	char*
	AnyOption::getArgv( int index )
	{
		if( index < new_argc ){
			return ( argv[ new_argv[ index ] ] );
		}
		return NULL;
	}

	/* dotfile sub routines */

	bool
	AnyOption::processFile()
	{
		if( !FileSet() )
			return false;
		return  ( consumeFile(readFile()) );
	}

	bool
	AnyOption::processFile( const char *filename )
	{
		useFilename(filename );
		return ( processFile() );
	}

	char*
	AnyOption::readFile()
	{
		return ( readFile(filename) );
	}

	/*
	 * read the file contents to a character buffer
	 */

	char*
	AnyOption::readFile( const char* fname )
	{
			int length;
			char *buffer;
			ifstream is;
			is.open ( fname , ifstream::in );
			if( ! is.good() ){
					is.close();
					return NULL;
			}
			is.seekg (0, ios::end);
			length = is.tellg();
			is.seekg (0, ios::beg);
			buffer = (char*) malloc(length*sizeof(char));
			is.read (buffer,length);
			is.close();
			return buffer;
	}

	/*
	 * scans a char* buffer for lines that does not
	 * start with the specified comment character.
	 */
	bool
	AnyOption::consumeFile( char *buffer )
	{

			if( buffer == NULL )
			return false;

			char *cursor = buffer;/* preserve the ptr */
			char *pline = NULL ;
			int linelength = 0;
			bool newline = true;
			for( unsigned int i = 0 ; i < strlen( buffer ) ; i++ ){
			if( *cursor == endofline ) { /* end of line */
				if( pline != NULL ) /* valid line */
						processLine( pline, linelength );
						pline = NULL;
						newline = true;
				}else if( newline ){ /* start of line */
						newline = false;
						if( (*cursor != comment ) ){ /* not a comment */
						pline = cursor ;
								linelength = 0 ;
						}
					}
					cursor++; /* keep moving */
					linelength++;
			}
			free (buffer);
		return true;
	}


	/*
	 *  find a valid type value pair separated by a delimiter
	 *  character and pass it to valuePairs()
	 *  any line which is not valid will be considered a value
	 *  and will get passed on to justValue()
	 *
	 *  assuming delimiter is ':' the behaviour will be,
	 *
	 *  width:10    - valid pair valuePairs( width, 10 );
	 *  width : 10  - valid pair valuepairs( width, 10 );
	 *
	 *  ::::        - not valid
	 *  width       - not valid
	 *  :10         - not valid
	 *  width:      - not valid
	 *  ::          - not valid
	 *  :           - not valid
	 *
	 */

	void
	AnyOption::processLine( char *theline, int length  )
	{
			bool found = false;
			char *pline = (char*) malloc( (length+1)*sizeof(char) );
			for( int i = 0 ; i < length ; i ++ )
					pline[i]= *(theline++);
			pline[length] = nullterminate;
			char *cursor = pline ; /* preserve the ptr */
			if( *cursor == delimiter || *(cursor+length-1) == delimiter ){
					justValue( pline );/* line with start/end delimiter */
			}else{
					for( int i = 1 ; i < length-1 && !found ; i++){/* delimiter */
							if( *cursor == delimiter ){
									*(cursor-1) = nullterminate; /* two strings */
									found = true;
									valuePairs( pline , cursor+1 );
							}
							cursor++;
					}
					cursor++;
					if( !found ) /* not a pair */
							justValue( pline );
			}
			free (pline);
	}

	void
	AnyOption::valuePairs( char *type, string value )
	{
		if ( breseq::chomp(type).size() == 1  ){ /* this is a char option */
			for( uint32_t i = 0 ; i < options.size() ; i++ ){
				if(  optionchars[i] == type[0]  ){ /* match */
					if( optiontype[i] == COMMON_OPT ||
						optiontype[i] == FILE_OPT )
					{
						setValue( type[0] , breseq::chomp(value) );
						return;
					}
				}
			}
		}
		/* if no char options matched */
		for( uint32_t i = 0 ; i < options.size() ; i++ ){
			if( options[i].compare(type) == 0 ){ /* match */
				if( optiontype[i] == COMMON_OPT ||
					optiontype[i] == FILE_OPT )
				{
					setValue( type , breseq::chomp(value) );
					return;
				}
			}
		}
			printVerbose( "Unknown option in resourcefile : " );
		printVerbose( type );
		printVerbose( );
	}

	void
	AnyOption::justValue( char *type )
	{

		if ( breseq::chomp(type).size() == 1  ){ /* this is a char option */
			for( uint32_t i = 0 ; i < options.size() ; i++ ){
				if(  optionchars[i] == type[0]  ){ /* match */
					if( optiontype[i] == COMMON_FLAG ||
						optiontype[i] == FILE_FLAG )
					{
						setFlagOn( type[0] );
						return;
					}
				}
			}
		}
		/* if no char options matched */
		for( uint32_t i = 0 ; i < options.size() ; i++ ){
			if( options[i].compare(type) == 0 ){ /* match */
				if( optiontype[i] == COMMON_FLAG ||
					optiontype[i] == FILE_FLAG )
				{
					setFlagOn( type );
					return;
				}
			}
		}
			printVerbose( "Unknown option in resourcefile : " );
		printVerbose( type  );
		printVerbose( );
	}

	/*
	 * usage and help
	 */


	void
	AnyOption::printAutoUsage()
	{
		if( autousage ) printUsage();
	}

	void
	AnyOption::printUsage()
	{

		if( once ) {
			once = false ;
			cout << endl ;
			for( size_t i = 0 ; i < usage_lines.size() ; i++ )
				cout << usage_lines[i] << endl ;
			cout << endl ;
		}
	}

	void
	AnyOption::addUsage( string line , bool advanced)
	{
    if(!advanced){
      usage_lines.push_back(line);
    }
    advanced_lines.push_back(line);
	}
  
  void
  AnyOption::addUsageSameLine( string line , bool advanced)
  {
    assert(usage_lines.size() > 0);
    assert(advanced_lines.size() > 0);

    if(!advanced){
      // add space if needed...
      if ( usage_lines.back()[usage_lines.back().size()-1] != ' ')
        usage_lines.back() += " ";
      usage_lines.back() += line;
    }
    
    if ( advanced_lines.back()[advanced_lines.back().size()-1] != ' ')
      advanced_lines.back() += " ";
    advanced_lines.back() += line;
  }
  
  void
	AnyOption::printAdvancedUsage()
	{
    if( once ) {
			once = false ;
      cout << endl ;
      for( size_t i = 0 ; i < advanced_lines.size() ; i++ )
        cout << advanced_lines[i] << endl ;
      cout << endl ;
    }
	}

	void
	AnyOption::addUsageError( string line )
	{
		cout << endl ;
		cout << "OPTIONS ERROR : Failed allocating extra memory " << endl ;
		cout << "While adding the usage/help  : \""<< line << "\"" << endl;
		cout << "Exiting." << endl ;
		cout << endl ;
		exit(0);

	}


	string AnyOption::operator[](const string& option_name)
	{
		if (getFlag(option_name))
			return TRUE_FLAG;
		if (getValue(option_name) != NULL)
      return *getValue(option_name);
		else if (default_values.count(option_name))
			return default_values[option_name];
		else
			return "";
	}

	bool AnyOption::count(const string& option_name)
	{
		return (getFlag(option_name) || getValue(option_name) != NULL);
	}

	string AnyOption::word_wrap(string sentence, int width)
	{
		string::iterator it = sentence.begin();
		string::iterator last_space = sentence.begin();
    
    int accumulated_width = 0;
		while (it != sentence.end())
		{
      accumulated_width++;
      if (*it == ' ') last_space = it;

			if (accumulated_width > width)
			{
				// Go back to letter after space
				it = last_space;
				*last_space = '\n';
        accumulated_width = 0;
			}
      it++;
		}

		return sentence;
	}

} // namespace breseq
