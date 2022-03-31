/*
 * ParamReader.h
 *
 *  Created on: Dec 14, 2017
 *      Author: vincent
 */

#ifndef PARAMREADER_H_
#define PARAMREADER_H_
#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <functional>
#include <locale>
#include <map>
using namespace std;
class ParamReader {
public:
	ParamReader();
	virtual ~ParamReader();
	ParamReader &open(const char* filename);
	std::map<string,double> numbers();
	std::map<string,string>strings();
private:
	string removeSpaces(string input);
	//int get_value(string s);
    string inputStr;
    std::vector<string> inputContent;
    bool is_number(const std::string& s);

};

template<class T>
class Params
{
public:
	Params(map<string,T> params)
	{
		this->params = params;
	};

	T get(string key,T def)
	{
		if(this->params.count(key)==0)
		{
			cout<<"missing param \""<<key<<"\",use default value "<<def<<endl;
			return def;
		}
		return this->params[key];
	};

private:
	map<string,T> params;
};


#endif /* PARAMREADER_H_ */
