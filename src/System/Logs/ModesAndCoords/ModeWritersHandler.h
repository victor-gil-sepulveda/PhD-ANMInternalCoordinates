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
#include "../../../ANM/Tools/Inout/ModesWriter.h"
class AtomSet;
class Unit;

class ModesWriterHandler {
	public:
				 ModesWriterHandler();
		virtual ~ModesWriterHandler();

		void setDirectory(std::string directory);
		void setCounter(unsigned int c);
		void setAtomSet(AtomSet* aset);

		std::string generatePath(std::string name);

		ModesWriter* getWriter(std::string name);
		void logStepAndVector(std::string logger_name, std::vector<double>& v);
		void logStepAndVector(std::string logger_name, std::vector<unsigned int>& v);
		void logStepAndCurrentCACoordinates(std::string logger_name);
		void logStepAndCAProposal(std::string logger_name, std::vector<double>& coord_increments);
		void logStepAndCurrentDihedralAngles(std::string logger_name, std::vector<Unit*>& units);
		void logStepAndDihedralAnglesProposal(std::string logger_name, std::vector<double>& target_angle_increments,
				std::vector<Unit*>& units);
		void logStepAndDihedralToCartesianProposal(std::string logger_name, std::vector<double>& target_angle_increments,
				std::vector<Unit*>& units);

		void saveCoordinates(std::string coordset_name, std::vector<double>& coords);
		void saveCACoordinates(std::string coordset_name);
		std::vector<double>& getCoordinates(std::string coordset_name);
		void removeCoordinates(std::string coordset_name);

	private:
		unsigned int counter;
		std::string directory;
		ModesWriter writer; // only one writer for each handler ! (we set the path again)

		std::map<std::string, std::vector<double> > coordinates_bank; // Temporary container for coordinates
		AtomSet* atomset;
};

#endif /* MODEWRITERSHANDLER_H_ */
