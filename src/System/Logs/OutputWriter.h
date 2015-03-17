/////////////////////////////////////////////////////////////////////////////
/// \file OutputWriter.h
/// \brief Class used to write output in log files
///
///
/// \copyright Copyright 2010-2014 BARCELONA SUPERCOMPUTING CENTER. See the COPYRIGHT file at the top-level directory of this distribution.
///
/// \author vgil
/// \date 18/07/2011
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef OUTPUTWRITER_H_
#define OUTPUTWRITER_H_
#include <string>

/////////////////////////////////////////////////////////////////////////////
/// \brief Class used to write output in log files
///
/// \author vgil
/// \date 18/07/2011
/////////////////////////////////////////////////////////////////////////////
class OutputWriter {

	public:
		OutputWriter(std::string path);
		virtual ~OutputWriter();

		// Redefine in a new class if you want to get different behaviors
		virtual	void 	write(std::string s);
		void 			writeWithDate(std::string s);
		void			flush();

	private:
		FILE* file;

		OutputWriter(const OutputWriter&);
		OutputWriter& operator=(const OutputWriter&);

};

#endif /* OUTPUTWRITER_H_ */
