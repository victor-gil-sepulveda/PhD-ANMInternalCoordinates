/////////////////////////////////////////////////////////////////////////////
/// \file bla,bla
///
/// \brief bla,bla
///
/// \author myName
/// \date 18/02/2015
/////////////////////////////////////////////////////////////////////////////

#pragma once

#ifndef MODEWRITERSHANDLER_H_
#define MODEWRITERSHANDLER_H_

#include <string>
#include <map>
#include <vector>
class ModesWriter;
class AtomSet;

class ModesWriterHandler {
	public:
				 ModesWriterHandler();
		virtual ~ModesWriterHandler();

		void setDirectory(std::string directory);
		void setCounter(unsigned int c);
		void setAtomSet(AtomSet* aset);
		std::string generatePath(std::string name);

		ModesWriter* getWriter(std::string name);

		void saveCoordinates(std::string coordset_name, std::vector<double>& coords);
		void saveCACoordinates(std::string coordset_name);

		std::vector<double>& getCoordinates(std::string coordset_name);

	private:
		unsigned int counter;
		std::string directory;
		std::map<std::string,ModesWriter*> writers; // Used writers
		std::map<std::string, std::vector<double> > coordinates_bank; // Temporary container for coordinates
		AtomSet* atomset;
};

#endif /* MODEWRITERSHANDLER_H_ */
