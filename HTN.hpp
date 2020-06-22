#pragma once

#include "HTN0.hpp"
#include "HTN1.hpp"
#include "HTN2.hpp"

namespace HTN {
	/* expression */

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t inverse(const htn_t& src) {
		return src.inverse();
	}

	template<typename lhs_t, typename rhs_t,
		typename std::enable_if < std::is_arithmetic < rhs_t > ::value > ::type* = nullptr,
		typename std::enable_if < std::is_base_of < HTN::htn_tag, lhs_t > ::value > ::type* = nullptr >
		lhs_t pow(const lhs_t& base, const rhs_t& index) {
		return base.pow(index);
	}

	template<typename lhs_t, typename rhs_t,
		typename std::enable_if < std::is_base_of < HTN::htn_tag, lhs_t > ::value > ::type* = nullptr,
		typename std::enable_if < std::is_base_of < HTN::htn_tag, rhs_t > ::value > ::type* = nullptr >
		lhs_t pow(const lhs_t& base, const rhs_t& index) {
		return base.pow(index);
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t sqrt(const htn_t& src) {
		return src.sqrt();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t exp(const htn_t& src) {
		return src.exp();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t log(const htn_t& src) {
		return src.log();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t sin(const htn_t& src) {
		return src.sin();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t cos(const htn_t& src) {
		return src.cos();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t tan(const htn_t& src) {
		return src.tan();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t asin(const htn_t& src) {
		return src.asin();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t acos(const htn_t& src) {
		return src.acos();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t atan(const htn_t& src) {
		return src.atan();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t sinh(const htn_t& src) {
		return src.sinh();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t cosh(const htn_t& src) {
		return src.cosh();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t tanh(const htn_t& src) {
		return src.tanh();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t asinh(const htn_t& src) {
		return src.asinh();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t acosh(const htn_t& src) {
		return src.acosh();
	}

	template<typename htn_t,
		typename std::enable_if < std::is_base_of<HTN::htn_tag, htn_t>::value>::type* = nullptr >
		htn_t atanh(const htn_t& src) {
		return src.atanh();
	}
}

/* expression with emmbedded arithmetic type */

template<typename lhs_t, typename rhs_t, 
	typename std::enable_if < std::is_arithmetic<lhs_t>::value>::type* = nullptr,
	typename std::enable_if < std::is_base_of<HTN::htn_tag, rhs_t>::value>::type* = nullptr >
	rhs_t operator+(const lhs_t& lhs, const rhs_t& rhs) {
	return rhs_t(lhs) + rhs;
}

template<typename lhs_t, typename rhs_t,
	typename std::enable_if < std::is_arithmetic < rhs_t > ::value > ::type* = nullptr,
	typename std::enable_if < std::is_base_of < HTN::htn_tag, lhs_t > ::value > ::type* = nullptr >
	lhs_t operator+(const lhs_t& lhs, const rhs_t& rhs) {
	return lhs + lhs_t(rhs);
}



template<typename lhs_t, typename rhs_t,
	typename std::enable_if < std::is_arithmetic<lhs_t>::value>::type* = nullptr,
	typename std::enable_if < std::is_base_of<HTN::htn_tag, rhs_t>::value>::type* = nullptr >
	rhs_t operator-(const lhs_t& lhs, const rhs_t& rhs) {
	return rhs_t(lhs) - rhs;
}

template<typename lhs_t, typename rhs_t,
	typename std::enable_if < std::is_arithmetic < rhs_t > ::value > ::type* = nullptr,
	typename std::enable_if < std::is_base_of < HTN::htn_tag, lhs_t > ::value > ::type* = nullptr >
	lhs_t operator-(const lhs_t& lhs, const rhs_t& rhs) {
	return lhs - lhs_t(rhs);
}



template<typename lhs_t, typename rhs_t,
	typename std::enable_if < std::is_arithmetic<lhs_t>::value>::type* = nullptr,
	typename std::enable_if < std::is_base_of<HTN::htn_tag, rhs_t>::value>::type* = nullptr >
	rhs_t operator*(const lhs_t& lhs, const rhs_t& rhs) {
	return rhs_t(lhs) * rhs;
}

template<typename lhs_t, typename rhs_t,
	typename std::enable_if < std::is_arithmetic < rhs_t > ::value > ::type* = nullptr,
	typename std::enable_if < std::is_base_of < HTN::htn_tag, lhs_t > ::value > ::type* = nullptr >
	lhs_t operator*(const lhs_t& lhs, const rhs_t& rhs) {
	return lhs * lhs_t(rhs);
}



template<typename lhs_t, typename rhs_t,
	typename std::enable_if < std::is_arithmetic<lhs_t>::value>::type* = nullptr,
	typename std::enable_if < std::is_base_of<HTN::htn_tag, rhs_t>::value>::type* = nullptr >
	rhs_t operator/(const lhs_t& lhs, const rhs_t& rhs) {
	return rhs_t(lhs) / rhs;
}

template<typename lhs_t, typename rhs_t,
	typename std::enable_if < std::is_arithmetic < rhs_t > ::value > ::type* = nullptr,
	typename std::enable_if < std::is_base_of < HTN::htn_tag, lhs_t > ::value > ::type* = nullptr >
	lhs_t operator/(const lhs_t& lhs, const rhs_t& rhs) {
	return lhs / lhs_t(rhs);
}


