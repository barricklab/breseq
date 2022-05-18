/*****************************************************************************
 
 AUTHORS
 
 Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 David B. Knoester
 
 LICENSE AND COPYRIGHT
 
 Copyright (c) 2008-2010 Michigan State University
 Copyright (c) 2011-2022 The University of Texas at Austin
 
 breseq is free software; you can redistribute it and/or modify it under the
 terms the GNU General Public License as published by the Free Software
 Foundation; either version 1, or (at your option) any later version.
 
 *****************************************************************************/

#ifndef _BRESEQ_FILE_PARSE_ERRORS_H_
#define _BRESEQ_FILE_PARSE_ERRORS_H_

namespace breseq {

  /*! Parse errors
   
   Helper class for aggregating all of the errors encountered in a file for output at one time
   
   */
  
  class cFileParseErrors {
  public:
    
    class sFileParseError {
    public:
      uint32_t _line_number;
      string _line;
      string _error_message;
      
      sFileParseError(uint32_t line_number, const string& line, const string& error_message )
      :_line_number(line_number), _line(line), _error_message(error_message)
      {}
      
      bool operator <(const sFileParseError& compare) const
      { return this->_line_number < compare._line_number; }
      
      void print()
      {
        string tab_line = substitute(_line, "\t", "<tab>");
        
        cerr << endl << ">>ERROR: " << _error_message << endl;
        if ((_line_number != 0) || (_line != "")) {
          cerr << ">>ON LINE: " << setw(5) <<  _line_number << endl << tab_line << endl;
        }
      }
    };
    
    // List pairs of line number and error message
    list<sFileParseError> _errors;
    string _filename;
    bool _fatal;
    
    cFileParseErrors(const string& filename = "")
    : _filename(filename), _fatal(false)
    { }
    
    void add_line_error(const uint32_t line_number, const string& line, const string& message, bool fatal)
    {
      _fatal = _fatal || fatal;
      _errors.push_back( sFileParseError(line_number, line, message) );
    }
    
    void print_errors(bool print_file = true)
    {
      if (_errors.size() == 0) return;
      
      if (print_file) {
        cerr << endl;
        cerr << ">>> Error(s) in GenomeDiff format. FILE: " << _filename << " <<<" << endl;
      }
      
      _errors.sort();
      
      for(list<sFileParseError>::iterator it = _errors.begin(); it != _errors.end(); it++) {
        it->print();
      }
      cerr << endl;
    }
    
    bool fatal() {return _fatal;}
  };

  
  
}
#endif
