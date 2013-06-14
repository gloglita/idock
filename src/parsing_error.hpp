#pragma once
#ifndef IDOCK_PARSING_ERROR_HPP
#define IDOCK_PARSING_ERROR_HPP

#include <string>
#include <stdexcept>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include "common.hpp"

using std::domain_error;
using boost::filesystem::path;
using boost::lexical_cast;

/// Represents a parsing error.
class parsing_error : public domain_error
{
public:
	/// Constructs a parsing error.
	parsing_error(const path& file, const size_t line, const string& reason) : domain_error("Error parsing \"" + file.filename().string() + "\" on line " + lexical_cast<string>(line) + ": " + reason) {}
};

#endif
