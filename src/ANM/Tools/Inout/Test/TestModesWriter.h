/////////////////////////////////////////////////////////////////////////////
/// \file TestModesWriter.h
/// \author vgil
/// \date 26/02/2015
/////////////////////////////////////////////////////////////////////////////

#ifndef TESTMODESWRITER_H_
#define TESTMODESWRITER_H_

#include <iosfwd>
#include "../../../../Tools/TestTools.h"


class TestModesWriter : public Test {
	public:
		TestModesWriter(std::string name);
		virtual ~TestModesWriter();

		void run();
		void init();

	  private:
		void finish();

		bool testCAFiltering(const char* prot_path,
				const char* expected_indices);

		bool testEigenFiltering(const char* indices_to_keep,
				const char* eigenvectors_file,
				const char* expected_eigenvectors);
};

#endif /* TESTMODESWRITER_H_ */
