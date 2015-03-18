
#ifndef __MODE_TYPES_H_
#define __MODE_TYPES_H_

#include <string>

namespace ModeTypes{

	enum ModeType{
		ANM_CC,
		ANM_IC,
		PCA_CC,
		PCA_IC,
		UNDEF
	};

	ModeType guessType(std::string type_str){
		if (type_str == "cc:anm")
			return ANM_CC;

		if (type_str == "ic:anm")
			return ANM_IC;

		if (type_str == "cc:pca")
			return PCA_CC;

		if (type_str == "ic:pca")
			return PCA_IC;

		return UNDEF;
	}

	std::string toString(ModeType type){
		if (type ==  ANM_CC) return "cc:anm";
		if (type ==  ANM_IC) return "ic:anm";
		if (type ==  PCA_CC) return "cc:pca";
		if (type ==  PCA_IC) return "ic:pca";
	}

	bool isCartesian(ModeType type){
		if (type == ANM_CC || type == PCA_CC){
			return true;
		}
		else{
			return false;
		}
	}

	bool isInternals(ModeType type){
		if (type == ANM_IC || type == PCA_IC){
			return true;
		}
		else{
			return false;
		}
	}

	bool isANM(ModeType type){
		if (type == ANM_CC || type == ANM_IC){
			return true;
		}
		else{
			return false;
		}
	}

	bool isPCA(ModeType type){
		if (type == PCA_CC || type == PCA_IC){
			return true;
		}
		else{
			return false;
		}
	}
};

#endif
