  /*! \file ConfigFileParser.cpp
      \author Henri Hugo Sieber <henri.hugo.sieber@cern.ch>
      \brief Source for an input string parsing class
  */

  #include "ConfigFileParser.hpp"

  // STD Library
  #include <iostream>
  #include <utility>
  #include <stdlib.h>
  #include <iomanip>
  #include <typeinfo>

  // Boost
  #include <boost/property_tree/json_parser.hpp>
  #include <boost/optional.hpp>
  #include <boost/foreach.hpp>
  #include <boost/throw_exception.hpp>

  ConfigFileParser* ConfigFileParser::fParser = nullptr;

  /*! \brief Printing utility for ConfigFileParser member
  
      \param key The tree key
      \param tree The Boost property tree 
      \param shift Tabulation shift
  */
  void PrintValues( const std::string& key
                  , const boost::property_tree::ptree& tree
                  , bool shift = false)
  {
    std::string tab = shift ? "\t" : "";

    // "key : entry"
    if (tree.empty()) 
      std::cout << tab << std::setw(35) << std::left << std::setfill('.') << key << ": " 
        << std::right << tree.data() << std::endl;  
    // "key : array"
    else if (tree.size() > 0 && tree.begin()->first == "") {
      std::cout << tab << std::setw(35) << std::left << std::setfill('.') << key << ": "; 
      for (auto it = tree.begin(); it != tree.end(); it++)
      {
        if (std::distance(tree.begin(), it) < static_cast<int>(tree.size()-1)) 
          std::cout << std::right << it->second.data() << ", ";
        else break;
      }
      std::cout << std::right << tree.rbegin()->second.data() << std::endl;
    }
  }

  ConfigFileParser::ConfigFileParser()
  {
    // ...
  }

  ConfigFileParser::ConfigFileParser(std::string ConfigFile)
  {
    fConfigFile = ConfigFile;
    boost::property_tree::read_json(fConfigFile, fPtree);
  }

  ConfigFileParser::~ConfigFileParser()
  {
    // ...
  }

  ConfigFileParser* ConfigFileParser::GetConfigFileParser()
  {
    if (!fParser)
      fParser = new ConfigFileParser();
    return fParser;
  }

  void 
  ConfigFileParser::SetConfigFile(std::string ConfigFile)
  {
    fConfigFile = ConfigFile;
    boost::property_tree::read_json(fConfigFile, fPtree);
  }

  void 
  ConfigFileParser::PrintConfigFile()
  {
    // NOTE: only supports one sub-tree for now.

    std::cout << "Using configuration file \033[1;33m" << fConfigFile 
      << "\033[0m with inputs:" << std::endl;
    BOOST_FOREACH(boost::property_tree::ptree::value_type& value, fPtree) {
      const std::string& key = value.first;
      const boost::property_tree::ptree& tree = value.second;

      // comment
      if (key.find("//") != std::string::npos) continue;

      PrintValues(key, tree);
      // "key : sub-tree"
      if (tree.size() > 0 && tree.begin()->first != "") {
        std::cout << std::setw(35) << std::left << std::setfill('.') << key << ". \n"; 
        for (auto it = tree.begin(); it != tree.end(); it++)
        {
          const std::string& subkey = it->first;
          const boost::property_tree::ptree& subtree = it->second;

          // comment
          if (subkey.find("//") != std::string::npos) continue;
          PrintValues(subkey, subtree, false);
        }
      }
    }

  }
  #if ENABLE_DEFAULT == 1
  template<typename T> T
  ConfigFileParser::GetConfigValue( std::string name, T _default )
  {
    if (boost::optional<T> value = fPtree.get_optional<T>(name)) {
      return *value;
    }
    else {

      std::cout << "Info in <ConfigFileParser::GetConfigValue>: "
                << "Using default value for <" << name << "> " << _default << "\n";

      return _default;
    }
  }

  template<typename T> std::vector<T>
  ConfigFileParser::GetConfigValues( std::string name
                                  , std::vector<T> _default )
  {
    std::vector<T> values;
    try{
      for (boost::property_tree::ptree::value_type& val : fPtree.get_child(name)) 
      {
        // TODO: add test for safe type conversion
        values.push_back(val.second.get_value<T>());
      }
      return values;
    }
    //catch (const boost::wrapexcept<boost::property_tree::ptree_bad_path>& e) {  
    catch (const std::exception& e) {

      std::cout << "Info in <ConfigFileParser::GetConfigValues>: "
                << "Using default value for <" << name << "> ";
      for( typename std::vector<T>::iterator it = values.begin();
          it != values.end(); it++ ) std::cout << *it << " "; 
      std::cout << "\n";

      return _default;
    }
  }

  template<typename T> T
  ConfigFileParser::GetConfigValueFromList( std::string name
                                          , std::string subname
                                          , T _default )
  {

    boost::property_tree::ptree tree; 
    try {
      T value;
      for (boost::property_tree::ptree::value_type& val : fPtree.get_child(name)) 
      {
        // get the label
        if (val.first == subname) {
          value = val.second.get_value<T>();
          break;
        }
      }
      return value;
    }
    //catch (const boost::wrapexcept<boost::property_tree::ptree_bad_path>& e) {  
    catch (const std::exception& e) {

      std::cout << "Info in <ConfigFileParser::GetConfigValue>: "
                << "Using default value for <" << name << "::" << subname << "> "
                << _default << "\n";

      return _default;
    }

  }

  template<typename T> std::vector<T>
  ConfigFileParser::GetConfigValuesFromList( std::string name
                                          , std::string subname
                                          , std::vector<T> _default )
  {
    std::vector<T> values;

    boost::property_tree::ptree tree; 
    try {
      tree = fPtree.get_child(name);
      const boost::property_tree::ptree& subtree = tree.get_child(subname);

      for (auto it = subtree.begin(); it != subtree.end(); it++) {
        values.push_back(it->second.get_value<T>());
      }

      return values;
    }
    //catch (const boost::wrapexcept<boost::property_tree::ptree_bad_path>& e) {  
    catch (const std::exception& e) {

      std::cout << "Info in <ConfigFileParser::GetConfigValuesFromList>: "
                << "Using default value for <" << name << "::" << subname << "> ";
      for( typename std::vector<T>::iterator it = values.begin();
          it != values.end(); it++ ) std::cout << *it << " "; 
      std::cout << "\n";

      return _default;
    }

  }
  #else
  template<typename T> T
  ConfigFileParser::GetConfigValue( std::string name )
  {
    if (boost::optional<T> value = fPtree.get_optional<T>(name)) {
      return *value;
    }
    else {
      std::cerr << "Error in <ConfigFileParser::GetConfigValue>:\n "
                << name << " does not referred to a node of the tree." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  template<typename T> std::vector<T>
  ConfigFileParser::GetConfigValues( std::string name )
  {
    std::vector<T> values;
    try{
      for (boost::property_tree::ptree::value_type& val : fPtree.get_child(name)) 
      {
        // TODO: add test for safe type conversion
        values.push_back(val.second.get_value<T>());
      }
    }
    //catch (const boost::wrapexcept<boost::property_tree::ptree_bad_path>& e) {  
    catch (const std::exception& e) {
      std::cerr << e.what() << std::endl;
      std::cerr << "Error in <ConfigFileParser::GetConfigValues>: \n"
                << name << " does not correspond to a node of the tree." << std::endl;
      exit(EXIT_FAILURE);
    }

    return values;
  }

  template<typename T> T
  ConfigFileParser::GetConfigValueFromList( std::string name
                                          , std::string subname)
  {

    boost::property_tree::ptree tree; 
    try {
      T value;
      for (boost::property_tree::ptree::value_type& val : fPtree.get_child(name)) 
      {
        // get the label
        if (val.first == subname) {
          value = val.second.get_value<T>();
          break;
        }
      }
      return value;
    }
    //catch (const boost::wrapexcept<boost::property_tree::ptree_bad_path>& e) {  
    catch (const std::exception& e) {
      std::cerr << e.what() << std::endl;
      std::cerr << "Error in <ConfigFileParser::GetConfigValueFromList>: \n"
                << subname << " does not correspond to a node of the tree." << std::endl;
      exit(EXIT_FAILURE);
    }

  }

  template<typename T> std::vector<T>
  ConfigFileParser::GetConfigValuesFromList( std::string name
                                          , std::string subname)
  {
    std::vector<T> values;

    boost::property_tree::ptree tree; 
    try {
      tree = fPtree.get_child(name);
      const boost::property_tree::ptree& subtree = tree.get_child(subname);

      for (auto it = subtree.begin(); it != subtree.end(); it++) {
        values.push_back(it->second.get_value<T>());
      }

      return values;
    }
    //catch (const boost::wrapexcept<boost::property_tree::ptree_bad_path>& e) {  
    catch (const std::exception& e) {
      std::cerr << e.what() << std::endl;
      std::cerr << "Error in <ConfigFileParser::GetConfigValuesFromList>: \n"
                << subname << " does not correspond to a node of the tree." << std::endl;
      exit(EXIT_FAILURE);
    }

  }
  #endif 

  // Explicit template instantation
  #if ENABLE_DEFAULT == 1
  template bool                     ConfigFileParser::GetConfigValue(std::string, bool);
  template int                      ConfigFileParser::GetConfigValue(std::string, int);
  template double                   ConfigFileParser::GetConfigValue(std::string, double);
  template std::string              ConfigFileParser::GetConfigValue(std::string, std::string);

  template std::vector<bool>        ConfigFileParser::GetConfigValues(std::string, std::vector<bool>);
  template std::vector<int>         ConfigFileParser::GetConfigValues(std::string, std::vector<int>);
  template std::vector<double>      ConfigFileParser::GetConfigValues(std::string, std::vector<double>);
  template std::vector<std::string> ConfigFileParser::GetConfigValues(std::string, std::vector<std::string>);

  template bool                     ConfigFileParser::GetConfigValueFromList( std::string
                                                                            , std::string
                                                                            , bool );
  template int                      ConfigFileParser::GetConfigValueFromList( std::string
                                                                            , std::string
                                                                            , int );
  template double                   ConfigFileParser::GetConfigValueFromList( std::string
                                                                            , std::string
                                                                            , double );
  template std::string              ConfigFileParser::GetConfigValueFromList( std::string
                                                                            , std::string
                                                                            , std::string );

  template std::vector<bool>        ConfigFileParser::GetConfigValuesFromList( std::string
                                                                            , std::string
                                                                            , std::vector<bool> );
  template std::vector<int>         ConfigFileParser::GetConfigValuesFromList( std::string
                                                                            , std::string
                                                                            , std::vector<int> );
  template std::vector<double>      ConfigFileParser::GetConfigValuesFromList( std::string
                                                                            , std::string
                                                                            , std::vector<double> );
  template std::vector<std::string> ConfigFileParser::GetConfigValuesFromList( std::string
                                                                            , std::string
                                                                            , std::vector<std::string> );
  #else 
  template bool                     ConfigFileParser::GetConfigValue(std::string);
  template int                      ConfigFileParser::GetConfigValue(std::string);
  template double                   ConfigFileParser::GetConfigValue(std::string);
  template std::string              ConfigFileParser::GetConfigValue(std::string);

  template std::vector<bool>        ConfigFileParser::GetConfigValues(std::string);
  template std::vector<int>         ConfigFileParser::GetConfigValues(std::string);
  template std::vector<double>      ConfigFileParser::GetConfigValues(std::string);
  template std::vector<std::string> ConfigFileParser::GetConfigValues(std::string);

  template bool                     ConfigFileParser::GetConfigValueFromList(std::string,std::string);
  template int                      ConfigFileParser::GetConfigValueFromList(std::string,std::string);
  template double                   ConfigFileParser::GetConfigValueFromList(std::string,std::string);
  template std::string              ConfigFileParser::GetConfigValueFromList(std::string,std::string);

  template std::vector<bool>        ConfigFileParser::GetConfigValuesFromList(std::string,std::string);
  template std::vector<int>         ConfigFileParser::GetConfigValuesFromList(std::string,std::string);
  template std::vector<double>      ConfigFileParser::GetConfigValuesFromList(std::string,std::string);
  template std::vector<std::string> ConfigFileParser::GetConfigValuesFromList(std::string,std::string);
  #endif
