#pragma once

#include <type_traits>

//#include "Eigen/Core"

#define VIENNACL_WITH_EIGEN
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
//#include "viennacl/hankel_matrix.hpp"

namespace HTN {

	constexpr static double PI() {
		return 3.1415926535897932384;
	}
	constexpr static double NAPIER() {
		return 2.7182818284590452353;
	}

	struct htn_tag {};

	/*
		No Differentiation
	    template
		- s_type : "float" or "double"
		- dim    : The size of the variable you want to use (0 ... dim-1)
	*/
	template<typename s_type, unsigned int dim = 1,
		typename std::enable_if<std::is_floating_point<s_type>::value>::type* = nullptr>
	class htn_r0 : htn_tag {
	public:
		using this_type_s = htn_r0<s_type, dim>;
		using inner_type_s = s_type;
		friend class htn_r0<s_type, dim>;

		template<unsigned int id = 0>
		static this_type_s x(s_type scalar) {
			return this_type_s(scalar);
		}

		static this_type_s x(s_type scalar, unsigned int id = 0) {
			return this_type_s(scalar);
		}

		htn_r0() = default;
		htn_r0(s_type init) : x0(init) {}
		//htn_r0(const inner_type_s& init) : x0(init) {}

		//operator type() {
		//	return (type)x0;
		//}

		this_type_s& operator=(const this_type_s& src) {
			x0 = src.x0;
			return *this;
		}

		this_type_s& operator=(s_type src) {
			x0 = src;
			return *this;
		}

		this_type_s operator+()  const {
			return *this;
		}

		this_type_s operator-() const {
			return create_this_type_s(minus0r());
		}

		this_type_s operator+(const this_type_s& rhs)  const {
			return create_this_type_s(add0r(rhs));
		}

		this_type_s operator-(const this_type_s& rhs) const {
			return create_this_type_s(sub0r(rhs));
		}

		this_type_s operator*(const this_type_s& rhs)  const {
			return create_this_type_s(mul0r(rhs));
		}

		this_type_s operator/(const this_type_s& rhs)  const {
			return create_this_type_s(mul0r(rhs.inverse()));
		}

		this_type_s& operator+=(const this_type_s& rhs) {
			return (*this = *this + rhs);
		}

		this_type_s& operator-=(const this_type_s& rhs) {
			return (*this = *this - rhs);
		}

		this_type_s& operator*=(const this_type_s& rhs) {
			return (*this = *this * rhs);
		}

		this_type_s& operator/=(const this_type_s& rhs) {
			return (*this = *this / rhs);
		}

		this_type_s& operator+=(s_type rhs) {
			return (*this = *this + this_type_s(rhs));
		}

		this_type_s& operator-=(s_type rhs) {
			return (*this = *this - this_type_s(rhs));
		}

		this_type_s& operator*=(s_type rhs) {
			return (*this = *this * this_type_s(rhs));
		}

		this_type_s& operator/=(s_type rhs) {
			return (*this = *this / this_type_s(rhs));
		}

		this_type_s inverse() const {
			return create_this_type_s(inv0r());
		}

		this_type_s pow(s_type index) const {
			return create_this_type_s(pow0r(index));
		}

		this_type_s pow(const this_type_s& index) const {
			return create_this_type_s(pow0r(index.x0));
		}

		this_type_s sqrt() const {
			return this->pow(static_cast<s_type>(0.5));
		}

		this_type_s exp() const {
			return create_this_type_s(exp0r());
		}

		this_type_s log() const {
			return create_this_type_s(log0r());
		}

		this_type_s sin() const {
			return create_this_type_s(sin0r());
		}

		this_type_s cos() const {
			return create_this_type_s(cos0r());
		}

		this_type_s tan() const {
			return create_this_type_s(tan0r());
		}

		this_type_s asin() const {
			return create_this_type_s(asin0r());
		}

		this_type_s acos() const {
			return create_this_type_s(acos0r());
		}

		this_type_s atan() const {
			return create_this_type_s(atan0r());
		}

		this_type_s sinh() const {
			return this_type_s(std::sinh(x0));
		}

		this_type_s cosh() const {
			return this_type_s(std::cosh(x0));
		}

		this_type_s tanh() const {
			return this_type_s(std::tanh(x0));
		}

		this_type_s asinh() const {
			return this_type_s(std::asinh(x0));
		}

		this_type_s acosh() const {
			return this_type_s(std::acosh(x0));
		}

		this_type_s atanh() const {
			return this_type_s(std::atanh(x0));
		}


		inner_type_s& X0() {
			return x0;
		}

		inner_type_s X0() const {
			return x0;
		}

		s_type X0_to_type() const {
			return x0;
		}

		constexpr static unsigned int get_dim() {
			return dim;
		}

	protected:
		inner_type_s x0 = 0;

		this_type_s create_this_type_s(const inner_type_s& src) const {
			return this_type_s(src);
		}

		/*Scalar minus operator of ViennaCL do not work. Bug? -> Substitute with arithmetic sub operator.*/
		inner_type_s minus0r()  const {
			return 0 - x0;
		}

		inner_type_s add0r(const this_type_s& rhs)  const {
			return x0 + rhs.x0;
		}

		inner_type_s sub0r(const this_type_s& rhs)  const {
			return x0 - rhs.x0;
		}

		inner_type_s mul0r(const this_type_s& rhs)  const {
			return x0 * rhs.x0;
		}

		inner_type_s inv0r()  const {
			assert(x0 != .0);

			return 1.0 / x0;
		}

		inner_type_s pow0r(s_type index) const {
			/*infinity error*/
			assert(!(x0 == .0 && index < 0));
			/*imaginary number error*/
			assert(!(x0 < .0 && std::floor(index) != index));

			if (index == 0) return inner_type_s(1);
			return std::pow((s_type)x0, index);
		}

		inner_type_s pow0r2(const this_type_s& index) const {
			/*infinity error*/
			assert(!(x0 == .0 && index.x0 < 0));
			/*imaginary number error*/
			assert(!(x0 < .0 && std::floor(index.x0) != index.x0));

			if ((s_type)index.x0 == 0) return inner_type_s(1);
			return std::pow((s_type)x0, (s_type)index.x0);
		}

		inner_type_s exp0r() const {
			return std::exp(x0);
		}

		inner_type_s log0r() const {
			assert(x0 > .0);
			return std::log((s_type)x0);
		}

		inner_type_s sin0r() const {
			return std::sin((s_type)x0);
		}

		inner_type_s cos0r() const {
			return std::cos((s_type)x0);
		}

		inner_type_s tan0r() const {
			return std::tan((s_type)x0);
		}

		inner_type_s asin0r() const {
			return std::asin((s_type)x0);
		}

		inner_type_s acos0r() const {
			return std::acos((s_type)x0);
		}

		inner_type_s atan0r() const {
			return std::atan((s_type)x0);
		}

	};

}



