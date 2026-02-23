    /*! \file ConfigFileParser.hpp
        \author Henri Hugo Sieber <henri.hugo.sieber@cern.ch>
        \brief Header for an input string parsing class
    */

    #ifndef ConfigFileParser_hh
    #define ConfigFileParser_hh

    #define ENABLE_DEFAULT 1

    // STD Library
    #include <string>
    #include <map>
    #include <array>
    #include <vector>

    // Boost
    #include <boost/property_tree/ptree.hpp>
    #include <boost/property_tree/ptree_fwd.hpp>

    /*! \class ConfigFileParser
        \brief A class for parsing inputs from the command line.

        This class parses inputs from the command line. It is based on the Boost
        library and supports as input JSON format files. In its current version it
        only supports Boost property trees with one level of node with childrens which
        is referred to as a list within this class.
    */
    class ConfigFileParser{

    public:

        /*! \brief The default constructor. */
        ConfigFileParser();
        /*! \brief A constructor for direct initialization. 
            \param name The name of the input JSON configuration file (with .json extension).
        */
        ConfigFileParser(std::string name);
        /*! \brief The destructor. */
        ~ConfigFileParser();
        
        /*! \brief Getting the static instance of the class. */
        static ConfigFileParser* GetConfigFileParser();

        /*! \brief Setting the JSON configuration file.
            \param name The name of the input JSON configuration file (with .json extension).
        */
        void SetConfigFile(std::string name);
        /*! \brief Printing the content of the JSON configuration file. */
        void PrintConfigFile();

        // non-list types
    #if ENABLE_DEFAULT == 1
        /*! \brief Template function for parsing single value at node level 0
            
            \param name The name of the zero-th level node
            \param _default The default value in case not provided in the parsing file
        */
        template<typename T> T GetConfigValue(std::string name, T _default = T());

        /*! \brief Template function for parsing values at node level 0
            
            \param name The name of the zero-th level node
            \param _default The default value in case not provided in the parsing file
        */
        template<typename T> std::vector<T> GetConfigValues( std::string name
                                                        , std::vector<T> _default = std::vector<T>());

        /*! \brief Template function for parsing single value at note level 1
            
            \param name The name of the zero-th level node
            \param subname The name of the level 1 node
            \param _default The default value in case not provided in the parsing file
        */
        template<typename T> T GetConfigValueFromList( std::string name
                                                    , std::string subname
                                                    , T _default = T() );

        /*! \brief Template function for parsing values at note level 1
            
            \param name The name of the zero-th level node
            \param subname The name of the level 1 node
            \param _default The default value in case not provided in the parsing file
        */
        template<typename T> std::vector<T> GetConfigValuesFromList( std::string name
                                                                , std::string subname
                                                                , std::vector<T> _default = std::vector<T>() );
    #else
        /*! \brief Template function for parsing single value at node level 0
            
            \param name The name of the zero-th level node
        */
        template<typename T> T GetConfigValue(std::string name );

        /*! \brief Template function for parsing values at node level 0
            
            \param name The name of the zero-th level node
        */
        template<typename T> std::vector<T> GetConfigValues(std::string name);

        /*! \brief Template function for parsing single value at note level 1
            
            \param name The name of the zero-th level node
            \param subname The name of the level 1 node
        */
        template<typename T> T GetConfigValueFromList(std::string name, std::string subname);

        /*! \brief Template function for parsing values at note level 1
            
            \param name The name of the zero-th level node
            \param subname The name of the level 1 node
        */
        template<typename T> std::vector<T> GetConfigValuesFromList(std::string name, std::string subname);
    #endif

    private:
        /*! The static instance of the class.  */
        static ConfigFileParser* fParser; 
        /*! The JSON configuration file name (with .json extension). */
        std::string fConfigFile;     
        /*! The Boost property tree. */
        boost::property_tree::ptree fPtree; 
    };

    #endif // ConfigFileParser_hh
