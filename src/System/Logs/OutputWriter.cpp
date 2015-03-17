/////////////////////////////////////////////////////////////////////////////
/// OutputWriter.cpp
/// Implementation of OutputWriter class
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author vgil
/// \date 18/07/2011
/////////////////////////////////////////////////////////////////////////////

#include "OutputWriter.h"

#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib>

using namespace std;

///////////////////////////////////////////////////////////////
/// \remarks
/// Class constructor.
///
/// \param path [In] Path of the log file
///
/// \author vgil
/// \date 18/07/2011
///////////////////////////////////////////////////////////////
OutputWriter::OutputWriter(std::string path){

	if( path==" "){
		cout<<"Empty path. Exiting..."<<endl;
		exit(EXIT_FAILURE);
	}

	file = fopen(path.c_str(),"w");

	if(!file ){
		cout<<"Impossible to open log "<<path<<endl<<"Exiting..."<<endl;
		exit(EXIT_FAILURE);
	}

}

OutputWriter::~OutputWriter() {
	fclose(file);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It writes a string to the log file
///
/// \param s [In] String to write in the log file
///
/// \author vgil
/// \date 18/07/2011
///////////////////////////////////////////////////////////////
void OutputWriter::write(std::string s){
	fprintf(file,"%s\n",s.c_str());
}

///////////////////////////////////////////////////////////////
/// \remarks
/// It writes a string to the log file adding first the current date
///
/// \param s [In] String to write in the log file
///
/// \author vgil
/// \date 18/07/2011
///////////////////////////////////////////////////////////////
void OutputWriter::writeWithDate(std::string s){

	// From http://www.dreamincode.net/code/snippet1102.htm
	//Find the current time
	time_t curtime = time(0);

	//convert it to tm
	tm now=*localtime(&curtime);

	//BUFSIZ is standard macro that expands to a integer constant expression
	//that is greater then or equal to 256. It is the size of the stream buffer
	//used by setbuf()
	char dest[BUFSIZ]={0};

	//Format string determines the conversion specification's behavior
 	const char format[]="%A, %B %d %Y. %X: ";

	//strftime - converts date and time to a string
	strftime(dest, sizeof(dest)-1, format, &now);

	std::string date(dest,strlen(dest));

	write(date+" "+s);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Forces this log to be written to disk.
///
/// \author vgil
/// \date 19/09/2011
///////////////////////////////////////////////////////////////
void OutputWriter::flush(){
	fflush(file);
}
