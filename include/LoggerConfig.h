#ifndef LOGGER_CONFIG_H
#define LOGGER_CONFIG_H


#include "Utils.h"
#include "XmlConfig.h"
#include "Logger.h"

#include <string>
using namespace std;

namespace jdb{

	class LoggerConfig
	{
	public:
		LoggerConfig(){  }
		~LoggerConfig(){ }
		
		static Logger * makeLogger( XmlConfig * config, string nodePath ) {


			if ( config && config->nodeExists( nodePath ) ){
				string ll = config->getString( nodePath + ".logLevel" );
				string outputStream = config->getString( nodePath + ".outputStream" );
				return (new Logger( Logger::logLevelFromString( ll ) ) );
			}

			return (new Logger() );
		}

	};

}




#endif