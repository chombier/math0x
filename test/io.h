#ifndef MATH0X_TEST_IO_H
#define MATH0X_TEST_IO_H

#include <string>

namespace math0x {
	namespace test {

		inline static std::string ok() { return "\033[32;1mok\033[0m"; }
		inline static std::string fail() { return "\033[31;1mfail\033[0m"; }

	}
}


#endif
