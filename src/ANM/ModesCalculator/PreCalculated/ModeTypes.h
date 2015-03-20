
#ifndef __MODE_TYPES_H__
#define __MODE_TYPES_H__

#include <string>

namespace ModeTypes{

	enum ModeType{
		ANM_CC,
		ANM_IC,
		PCA_CC,
		PCA_IC,
		UNDEF
	};

	ModeType guessType(std::string type_str);

	std::string toString(ModeType type);

	bool isCartesian(ModeType type);

	bool isInternals(ModeType type);

	bool isANM(ModeType type);

	bool isPCA(ModeType type);
};

#endif
