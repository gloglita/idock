#include "qtn4.hpp"

qtn4::qtn4(const float q0, const float q1, const float q2, const float q3)
{
	(*this)[0] = q0;
	(*this)[1] = q1;
	(*this)[2] = q2;
	(*this)[3] = q3;
}
