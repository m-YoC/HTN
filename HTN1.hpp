#pragma once

#include "HTN0.hpp"

namespace HTN {

	/*
		Automatic Calculation of "Gradient" 
		template
		- s_type : "float" or "double"
		- dim    : The size of the variable you want to use (0 ... dim-1)
	*/
	template<typename s_type, unsigned int dim = 1,
		typename std::enable_if<std::is_floating_point<s_type>::value>::type* = nullptr>
	class htn_r1 : public htn_r0<s_type, dim>, htn_tag {
	public:
		using this_type_v = htn_r1<s_type, dim>;
		using inner_type_s = viennacl::scalar<s_type>;
		using inner_type_v = viennacl::vector<s_type>;
		using inner_type_v_with_bool = std::pair<bool, inner_type_v>;
		friend class  htn_r1<s_type, dim>;

		template<unsigned int id = 0>
		static this_type_v x(s_type scalar) {
			assert(id < dim);
			return this_type_v(inner_type_s(scalar), viennacl::unit_vector<s_type>(dim, id));
		}

		static this_type_v x(s_type scalar, unsigned int id = 0) {
			assert(id < dim);
			return this_type_v(inner_type_s(scalar), viennacl::unit_vector<s_type>(dim, id));
		}

		htn_r1() = default;
		using htn_r0<s_type, dim>::htn_r0;

		htn_r1(const inner_type_s& init0, const inner_type_v& init1) {
			this->x0 = init0;
			this->x1 = init1;
			this->is_zero_vector = false;
		}

		this_type_v& operator=(const this_type_v& src) {
			this->x0 = src.x0;
			this->x1 = src.x1;
			this->is_zero_vector = src.is_zero_vector;
			return *this;
		}

		this_type_v& operator=(s_type src) {
			this->x0 = src;
			if (!is_zero_vector) {
				this->x1.empty();
				is_zero_vector = true;
			}
			return *this;
		}

		this_type_v operator+()  const {
			return *this;
		}

		this_type_v operator-()  const {
			return create_this_type_v(this->minus0r(), this->minus1r());
		}

		this_type_v operator+(const this_type_v& rhs)  const {
			return create_this_type_v(this->add0r(rhs), this->add1r(rhs));
		}

		this_type_v operator-(const this_type_v& rhs)  const {
			return create_this_type_v(this->sub0r(rhs), this->sub1r(rhs));
		}

		this_type_v operator*(const this_type_v& rhs)  const {
			return create_this_type_v(this->mul0r(rhs), this->mul1r(rhs));
		}

		this_type_v operator/(const this_type_v& rhs)  const {
			auto rhsinv = rhs.inverse();
			return create_this_type_v(this->mul0r(rhsinv), this->mul1r(rhsinv));
		}

		this_type_v& operator+=(const this_type_v& rhs) {
			return (*this = *this + rhs);
		}

		this_type_v& operator-=(const this_type_v& rhs) {
			return (*this = *this - rhs);
		}

		this_type_v& operator*=(const this_type_v& rhs) {
			return (*this = *this * rhs);
		}

		this_type_v& operator/=(const this_type_v& rhs) {
			return (*this = *this / rhs);
		}

		this_type_v& operator+=(s_type rhs) {
			return (*this = *this + this_type_v(rhs));
		}

		this_type_v& operator-=(s_type rhs) {
			return (*this = *this - this_type_v(rhs));
		}

		this_type_v& operator*=(s_type rhs) {
			return (*this = *this * this_type_v(rhs));
		}

		this_type_v& operator/=(s_type rhs) {
			return (*this = *this / this_type_v(rhs));
		}

		this_type_v inverse() const {
			return create_this_type_v(this->inv0r(), this->inv1r());
		}

		this_type_v pow(s_type index) const {
			return create_this_type_v(this->pow0r(index), this->pow1r(index));
		}

		this_type_v pow(const this_type_v& index) const {
			if (index.is_zero_vector) {
				return create_this_type_v(this->pow0r((s_type)index.x0), this->pow1r((s_type)index.x0));
			}
			else {
				return create_this_type_v(this->pow0r((s_type)index.x0), this->pow1r2(index));
			}
		}

		this_type_v sqrt() const {
			return (*this).pow(static_cast<s_type>(0.5));
		}

		this_type_v exp() const {
			return create_this_type_v(this->exp0r(), this->exp1r());
		}

		this_type_v log() const {
			return create_this_type_v(this->log0r(), this->log1r());
		}

		this_type_v sin() const {
			return create_this_type_v(this->sin0r(), this->sin1r());
		}

		this_type_v cos() const {
			return create_this_type_v(this->cos0r(), this->cos1r());
		}

		this_type_v tan() const {
			return create_this_type_v(this->tan0r(), this->tan1r());
		}

		this_type_v asin() const {
			return create_this_type_v(this->asin0r(), this->asin1r());
		}

		this_type_v acos() const {
			return create_this_type_v(this->acos0r(), this->acos1r());
		}

		this_type_v atan() const {
			return create_this_type_v(this->atan0r(), this->atan1r());
		}

		this_type_v sinh() const {
			return ((*this).exp() - (-(*this)).exp()) * this_type_v(static_cast<s_type>(0.5));
		}

		this_type_v cosh() const {
			return ((*this).exp() + (-(*this)).exp()) * this_type_v(static_cast<s_type>(0.5));
		}

		this_type_v tanh() const {
			const auto exp_p = (*this).exp();
			const auto exp_m = (-(*this)).exp();
			return (exp_p - exp_m) / (exp_p + exp_m);
		}

		this_type_v asinh() const {
			const auto _1 = this_type_v(static_cast<s_type>(1));
			return ((*this) + ( this->pow(2) + _1 ).sqrt()).log();
		}

		this_type_v acosh(s_type pm = 1) const {
			const auto _1 = this_type_v(static_cast<s_type>(1));
			const auto _sgn = this_type_v(static_cast<s_type>((pm >= 0) - (pm < 0)));
			return ((*this) + _sgn * ( this->pow(2) - _1 ).sqrt() ).log();
		}

		this_type_v atanh() const {
			const auto _05 = this_type_v(static_cast<s_type>(0.5));
			const auto _1 = this_type_v(static_cast<s_type>(1));
			return ( (_1 + *this) / (_1 - *this) ).log() * _05;
		}

		inner_type_v& X1() {
			return x1;
		}

		inner_type_v X1() const {
			return x1;
		}

		Eigen::Matrix<s_type, -1, 1> X1_to_eig() const {
			Eigen::Matrix<s_type, -1, 1> res(dim);
			if (is_zero_vector) {
				res.setZero();
			}
			else {
				viennacl::copy(x1, res);
			}
			return res;
		}

	protected:
		inner_type_v x1;
		bool is_zero_vector = true;

		this_type_v create_this_type_v(const inner_type_s& src0, const inner_type_v_with_bool& src1) const {
			if (!src1.first) {
				return this_type_v(src0, src1.second);
			}
			else {
				return this_type_v(src0);
			}
		}

		inner_type_v_with_bool minus1r()  const {
			if (!is_zero_vector) {
				return inner_type_v_with_bool(false, -x1);
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool add1r(const this_type_v& rhs)  const {
			if (is_zero_vector == false && rhs.is_zero_vector == false) {
				return inner_type_v_with_bool(false, x1 + rhs.x1);
			}
			else if (is_zero_vector == false && rhs.is_zero_vector == true) {
				return inner_type_v_with_bool(false, x1);
			}
			else if (is_zero_vector == true && rhs.is_zero_vector == false) {
				return inner_type_v_with_bool(false, rhs.x1);
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool sub1r(const this_type_v& rhs)  const {
			if (is_zero_vector == false && rhs.is_zero_vector == false) {
				return inner_type_v_with_bool(false, x1 - rhs.x1);
			}
			else if (is_zero_vector == false && rhs.is_zero_vector == true) {
				return inner_type_v_with_bool(false, x1);
			}
			else if (is_zero_vector == true && rhs.is_zero_vector == false) {
				return inner_type_v_with_bool(false, -rhs.x1);
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool mul1r(const this_type_v& rhs)  const {
			if (is_zero_vector == false && rhs.is_zero_vector == false) {
				return inner_type_v_with_bool(false, rhs.x0 * x1 + this->x0 * rhs.x1);
			}
			else if (is_zero_vector == false && rhs.is_zero_vector == true) {
				return inner_type_v_with_bool(false, rhs.x0 * x1);
			}
			else if (is_zero_vector == true && rhs.is_zero_vector == false) {
				return inner_type_v_with_bool(false, this->x0 * rhs.x1);
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool inv1r()  const {
			assert(this->x0 != .0);
			if (!is_zero_vector) {
				return inner_type_v_with_bool(false, -x1 / (this->x0 * this->x0));
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool pow1r(s_type index)  const {
			/*infinity error*/
			assert(!(this->x0 == .0 && index < 1));
			/*imaginary number error*/
			assert(!(this->x0 < .0 && std::floor(index) != index));

			if (index == 0) return inner_type_v_with_bool(true, inner_type_v());

			if (!is_zero_vector) {
				return inner_type_v_with_bool(false, x1 * ( index * std::pow((s_type)this->x0, index - static_cast<s_type>(1)) ) );
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool pow1r2(const this_type_v& index)  const {
			/*infinity error*/
			assert(!(this->x0 == .0 && index.x0 < 1));
			/*imaginary number error*/
			assert(!(this->x0 < .0 && std::floor(index.x0) != index.x0));

			if (is_zero_vector == false && index.is_zero_vector == false) {
				const s_type c2 = (s_type)index.x0 * std::pow((s_type)this->x0, (s_type)index.x0 - static_cast<s_type>(1));
				const s_type c1 = std::pow((s_type)this->x0, (s_type)index.x0) * std::log((s_type)this->x0);
				return inner_type_v_with_bool(false, c2 * x1 + c1 * index.x1);
			}
			else if (is_zero_vector == false && index.is_zero_vector == true) {
				const s_type c2 = (s_type)index.x0 * std::pow((s_type)this->x0, (s_type)index.x0 - static_cast<s_type>(1));
				return inner_type_v_with_bool(false, c2 * x1);
			}
			else if (is_zero_vector == true && index.is_zero_vector == false) {
				const s_type c1 = std::pow((s_type)this->x0, (s_type)index.x0) * std::log((s_type)this->x0);
				return inner_type_v_with_bool(false, c1 * index.x1);
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool exp1r()  const {

			if (!is_zero_vector) {
				return inner_type_v_with_bool(false, std::exp((s_type)this->x0) * x1);
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool log1r()  const {
			assert(this->x0 > .0);

			if (!is_zero_vector) {
				return inner_type_v_with_bool(false, x1 / (s_type)this->x0);
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool sin1r()  const {

			if (!is_zero_vector) {
				return inner_type_v_with_bool(false, std::cos((s_type)this->x0) * x1);
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool cos1r()  const {

			if (!is_zero_vector) {
				return inner_type_v_with_bool(false, -std::sin((s_type)this->x0) * x1);
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool tan1r()  const {

			if (!is_zero_vector) {
				return inner_type_v_with_bool(false, (static_cast<s_type>(1) + std::pow(std::tan((s_type)this->x0), 2)) * x1);
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool asin1r()  const {

			if (!is_zero_vector) {
				return inner_type_v_with_bool(false, x1 / std::sqrt(static_cast<s_type>(1) - std::pow((s_type)this->x0, 2)) );
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool acos1r()  const {

			if (!is_zero_vector) {
				return inner_type_v_with_bool(false, -x1 / std::sqrt(static_cast<s_type>(1) - std::pow((s_type)this->x0, 2)));
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}

		inner_type_v_with_bool atan1r()  const {

			if (!is_zero_vector) {
				return inner_type_v_with_bool(false, x1 / (static_cast<s_type>(1) + std::pow((s_type)this->x0, 2)) );
			}
			else {
				return inner_type_v_with_bool(true, inner_type_v());
			}
		}


	};

	
}
