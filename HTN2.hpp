#pragma once

#include "HTN1.hpp"

namespace HTN {

	/*
		Automatic Calculation of "Gradient" and "Hessian Matrix" 
		template
		- s_type : "float" or "double"
		- dim    : The size of the variable you want to use (0 ... dim-1)
	*/
	template<typename s_type, unsigned int dim = 1,
		typename std::enable_if<std::is_floating_point<s_type>::value>::type* = nullptr>
	class htn_r2 : public htn_r1<s_type, dim>, htn_tag {
	public:
		using this_type_m = htn_r2<s_type, dim>;
		using inner_type_s = viennacl::scalar<s_type>;
		using inner_type_v = viennacl::vector<s_type>;
		using inner_type_v_with_bool = std::pair<bool, inner_type_v>;
		using inner_type_m = viennacl::matrix<s_type>;
		using inner_type_m_with_bool = std::pair<bool, inner_type_m>;
		friend class htn_r2<s_type, dim>;

		template<unsigned int id = 0>
		static this_type_m x(s_type scalar) {
			assert(id < dim);
			return this_type_m(inner_type_s(scalar), viennacl::unit_vector<s_type>(dim, id));
		}

		static this_type_m x(s_type scalar, unsigned int id = 0) {
			assert(id < dim);
			return this_type_m(inner_type_s(scalar), viennacl::unit_vector<s_type>(dim, id));
		}

		htn_r2() = default;
		//using htn0d<type, var_num>::htn0d;
		using htn_r1<s_type, dim>::htn_r1;

		htn_r2(const inner_type_s& init0, const inner_type_v& init1, const inner_type_m& init2) {
			this->x0 = init0;
			this->x1 = init1;
			this->x2 = init2;
			this->is_zero_vector = false;
			this->is_zero_matrix = false;
		}

		this_type_m& operator=(const this_type_m& src) {
			this->x0 = src.x0;
			this->x1 = src.x1;
			this->x2 = src.x2;
			this->is_zero_vector = src.is_zero_vector;
			this->is_zero_matrix = src.is_zero_matrix;
			return *this;
		}

		this_type_m& operator=(s_type src) {
			this->x0 = src;
			if (!this->is_zero_vector) {
				this->x1.empty();
				this->is_zero_vector = true;
			}
			if (!this->is_zero_matrix) {
				this->x2.empty();
				this->is_zero_matrix = true;
			}
			return *this;
		}

		this_type_m operator+()  const {
			return *this;
		}

		this_type_m operator-()  const {
			return create_this_type_m(this->minus0r(), this->minus1r(), this->minus2r());
		}

		this_type_m operator+(const this_type_m& rhs)  const {
			return create_this_type_m(this->add0r(rhs), this->add1r(rhs), this->add2r(rhs));
		}

		this_type_m operator-(const this_type_m& rhs)  const {
			return create_this_type_m(this->sub0r(rhs), this->sub1r(rhs), this->sub2r(rhs));
		}

		this_type_m operator*(const this_type_m& rhs)  const {
			return create_this_type_m(this->mul0r(rhs), this->mul1r(rhs), this->mul2r(rhs));
		}

		this_type_m operator/(const this_type_m& rhs)  const {
			auto rhsinv = rhs.inverse();
			return create_this_type_m(this->mul0r(rhsinv), this->mul1r(rhsinv), this->mul2r(rhsinv));
		}

		this_type_m& operator+=(const this_type_m& rhs) {
			return (*this = *this + rhs);
		}

		this_type_m& operator-=(const this_type_m& rhs) {
			return (*this = *this - rhs);
		}

		this_type_m& operator*=(const this_type_m& rhs) {
			return (*this = *this * rhs);
		}

		this_type_m& operator/=(const this_type_m& rhs) {
			return (*this = *this / rhs);
		}

		this_type_m& operator+=(s_type rhs) {
			return (*this = *this + this_type_m(rhs));
		}

		this_type_m& operator-=(s_type rhs) {
			return (*this = *this - this_type_m(rhs));
		}

		this_type_m& operator*=(s_type rhs) {
			return (*this = *this * this_type_m(rhs));
		}

		this_type_m& operator/=(s_type rhs) {
			return (*this = *this / this_type_m(rhs));
		}

		this_type_m inverse() const {
			return create_this_type_m(this->inv0r(), this->inv1r(), this->inv2r());
		}

		this_type_m pow(s_type index) const {
			return create_this_type_m(this->pow0r(index), this->pow1r(index), this->pow2r(index));
		}

		this_type_m pow(const this_type_m& index) const {
			if (index.is_zero_vector && index.is_zero_matrix) {
				return create_this_type_m(this->pow0r((s_type)index.x0), this->pow1r((s_type)index.x0), this->pow2r((s_type)index.x0));
			}
			else {
				return create_this_type_m(this->pow0r((s_type)index.x0), this->pow1r2(index), this->pow2r2(index));
			}
			
		}

		this_type_m sqrt() const {
			return (*this).pow(static_cast<s_type>(0.5));
		}

		this_type_m exp() const {
			return create_this_type_m(this->exp0r(), this->exp1r(), this->exp2r());
		}

		this_type_m log() const {
			return create_this_type_m(this->log0r(), this->log1r(), this->log2r());
		}

		this_type_m sin() const {
			return create_this_type_m(this->sin0r(), this->sin1r(), this->sin2r());
		}

		this_type_m cos() const {
			return create_this_type_m(this->cos0r(), this->cos1r(), this->cos2r());
		}

		this_type_m tan() const {
			return create_this_type_m(this->tan0r(), this->tan1r(), this->tan2r());
		}

		this_type_m asin() const {
			return create_this_type_m(this->asin0r(), this->asin1r(), this->asin2r());
		}

		this_type_m acos() const {
			return create_this_type_m(this->acos0r(), this->acos1r(), this->acos2r());
		}

		this_type_m atan() const {
			return create_this_type_m(this->atan0r(), this->atan1r(), this->atan2r());
		}

		this_type_m sinh() const {
			return ((*this).exp() - (-(*this)).exp()) * this_type_m(static_cast<s_type>(0.5));
		}

		this_type_m cosh() const {
			return ((*this).exp() + (-(*this)).exp()) * this_type_m(static_cast<s_type>(0.5));
		}

		this_type_m tanh() const {
			const auto exp_p = (*this).exp();
			const auto exp_m = (-(*this)).exp();
			return (exp_p - exp_m) / (exp_p + exp_m);
		}

		this_type_m asinh() const {
			const auto _1 = this_type_m(static_cast<s_type>(1));
			return ((*this) + (this->pow(2) + _1).sqrt()).log();
		}

		this_type_m acosh(s_type pm = 1) const {
			const auto _1 = this_type_m(static_cast<s_type>(1));
			const auto _sgn = this_type_m(static_cast<s_type>((pm >= 0) - (pm < 0)));
			return ((*this) + _sgn * (this->pow(2) - _1).sqrt()).log();
		}

		this_type_m atanh() const {
			const auto _05 = this_type_m(static_cast<s_type>(0.5));
			const auto _1 = this_type_m(static_cast<s_type>(1));
			return ((_1 + *this) / (_1 - *this)).log() * _05;
		}

		inner_type_m& X2() {
			return x2;
		}

		inner_type_m X2() const {
			return x2;
		}

		Eigen::Matrix<s_type, -1, -1> X2_to_eig() const {
			Eigen::Matrix<s_type, -1, -1> res(dim, dim);
			if (is_zero_matrix) {
				res.setZero();
			}
			else {
				viennacl::copy(x2, res);
			}
			return res;
		}

	protected:
		inner_type_m x2;
		bool is_zero_matrix = true;

		this_type_m create_this_type_m(const inner_type_s& src0, const inner_type_v_with_bool& src1, const inner_type_m_with_bool& src2)  const {
			if (!src1.first && !src2.first) {
				return this_type_m(src0, src1.second, src2.second);
			}
			else if (!src1.first && src2.first) {
				return this_type_m(src0, src1.second);
			}
			else if (src1.first && !src2.first) {
				return this_type_m(src0, viennacl::zero_vector<s_type>(dim), src2.second);
			}
			else {
				return this_type_m(src0);
			}
		}

		inner_type_m_with_bool minus2r()  const {
			if (!is_zero_matrix) {
				return inner_type_m_with_bool(false, -x2);
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
		}

		inner_type_m_with_bool add2r(const this_type_m& rhs)  const {
			if (is_zero_matrix == false && rhs.is_zero_matrix == false) {
				return inner_type_m_with_bool(false, x2 + rhs.x2);
			}
			else if (is_zero_matrix == false && rhs.is_zero_matrix == true) {
				return inner_type_m_with_bool(false, x2);
			}
			else if (is_zero_matrix == true && rhs.is_zero_matrix == false) {
				return inner_type_m_with_bool(false, rhs.x2);
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
		}

		inner_type_m_with_bool sub2r(const this_type_m& rhs)  const {
			if (is_zero_matrix == false && rhs.is_zero_matrix == false) {
				return inner_type_m_with_bool(false, x2 - rhs.x2);
			}
			else if (is_zero_matrix == false && rhs.is_zero_matrix == true) {
				return inner_type_m_with_bool(false, x2);
			}
			else if (is_zero_matrix == true && rhs.is_zero_matrix == false) {
				return inner_type_m_with_bool(false, -rhs.x2);
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
		}

		inner_type_m_with_bool mul2r(const this_type_m& rhs)  const {
			if (this->is_zero_vector && is_zero_matrix && rhs.is_zero_vector && rhs.is_zero_matrix) {
				return inner_type_m_with_bool(true, inner_type_m());
			}

			inner_type_m res(dim, dim);
			res.clear();

			if (!is_zero_matrix) {
				res += rhs.x0 * x2;
			}
			if (!rhs.is_zero_matrix) {
				res += this->x0 * rhs.x2;
			}
			if (!this->is_zero_vector && !rhs.is_zero_vector) {
				inner_type_m buf = viennacl::linalg::outer_prod(this->x1, rhs.x1);
				res += buf + viennacl::trans(buf);
			}

			return inner_type_m_with_bool(false, res);
		}

		inner_type_m_with_bool inv2r()  const {
			assert(this->x0 != .0);
			if (!this->is_zero_vector && !is_zero_matrix) {
				return inner_type_m_with_bool(false, (2.0 * viennacl::linalg::outer_prod(this->x1, this->x1) - this->x0 * x2) / (this->x0 * this->x0 * this->x0));
			}
			else if (!this->is_zero_vector && is_zero_matrix) {
				return inner_type_m_with_bool(false, (2.0 * viennacl::linalg::outer_prod(this->x1, this->x1)) / (this->x0 * this->x0 * this->x0));
			}
			else if (this->is_zero_vector && !is_zero_matrix) {
				return inner_type_m_with_bool(false, x2 / (-(s_type)this->x0 * (s_type)this->x0));
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
		}

		inner_type_m_with_bool pow2r(s_type index)  const {
			/*infinity error*/
			assert(!(this->x0 == .0 && index < 2));
			/*imaginary number error*/
			assert(!(this->x0 < .0 && std::floor(index) != index));

			if (index == 0) return inner_type_m_with_bool(true, inner_type_m());
			if (this->is_zero_vector && is_zero_matrix) return inner_type_m_with_bool(true, inner_type_m());
			if(is_zero_matrix && index == static_cast<s_type>(1)) return inner_type_m_with_bool(true, inner_type_m());

			inner_type_m res(dim, dim);

			if (!this->is_zero_vector) {
				res += viennacl::linalg::outer_prod(this->x1, this->x1) * (index * (index - static_cast<s_type>(1)) * std::pow((s_type)this->x0, index - static_cast<s_type>(2)));
			}

			if (!is_zero_matrix) {
				res += x2 * (index * std::pow((s_type)this->x0, index - static_cast<s_type>(1)));
			}

			return inner_type_m_with_bool(false, res);
		}

		inner_type_m_with_bool pow2r2(const this_type_m& index)  const {
			/*infinity error*/
			assert(!(this->x0 == .0 && index.x0 < 2));
			/*imaginary number error*/
			assert(!(this->x0 < .0 && std::floor(index.x0) != index.x0));

			bool no_add = true;
			inner_type_m res(dim, dim);

			if ((s_type)index.x0 != 0) {

				const s_type c1 = (s_type)index.x0 * std::pow((s_type)this->x0, (s_type)index.x0 - static_cast<s_type>(2));

				if (!this->is_zero_matrix) {
					res += c1 * (s_type)this->x0 * x2;
					no_add &= false;
				}
				if (!this->is_zero_vector && (s_type)index.x0 != static_cast<s_type>(1)) {
					res += c1 * ((s_type)index.x0 - static_cast<s_type>(1)) * viennacl::linalg::outer_prod(this->x1, this->x1);
					no_add &= false;
				}

			}

			const s_type logx = std::log((s_type)this->x0);
			const s_type c2 = std::pow((s_type)this->x0, (s_type)index.x0 * logx);
			if (!index.is_zero_matrix) {
				res += c2 * index.x2;
				no_add &= false;
			}
			if (!index.is_zero_vector) {
				res += logx * c2 * viennacl::linalg::outer_prod(index.x1, index.x1);
				no_add &= false;
			}

			if (!this->is_zero_vector && !index.is_zero_vector) {
				const s_type c3 = (static_cast<s_type>(1) + (s_type)index.x0 * logx) * std::pow((s_type)this->x0, (s_type)index.x0 - static_cast<s_type>(1));
				inner_type_m buf = viennacl::linalg::outer_prod(this->x1, index.x1);
				res += c3 * (buf + viennacl::trans(buf));
				no_add &= false;
			}

			if (!no_add) {
				return inner_type_m_with_bool(false, res);
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
			
		}

		inner_type_m_with_bool exp2r()  const {
			
			if (!this->is_zero_vector && !is_zero_matrix) {
				return inner_type_m_with_bool(false, std::exp((s_type)this->x0) * (x2 + viennacl::linalg::outer_prod(this->x1, this->x1)) );
			}
			else if (!this->is_zero_vector && is_zero_matrix) {
				return inner_type_m_with_bool(false, std::exp((s_type)this->x0) * viennacl::linalg::outer_prod(this->x1, this->x1));
			}
			else if (this->is_zero_vector && !is_zero_matrix) {
				return inner_type_m_with_bool(false, std::exp((s_type)this->x0) * x2);
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
		}

		inner_type_m_with_bool log2r()  const {
			assert(this->x0 > .0);

			if (!this->is_zero_vector && !is_zero_matrix) {
				inner_type_m outprpd = viennacl::linalg::outer_prod(this->x1, this->x1);
				return inner_type_m_with_bool(false, x2 / (s_type)this->x0 - outprpd / ((s_type)this->x0 * (s_type)this->x0));
			}
			else if (!this->is_zero_vector && is_zero_matrix) {
				inner_type_m outprpd = viennacl::linalg::outer_prod(this->x1, this->x1);
				return inner_type_m_with_bool(false, outprpd / (-(s_type)this->x0 * (s_type)this->x0) );
			}
			else if (this->is_zero_vector && !is_zero_matrix) {
				return inner_type_m_with_bool(false, x2 / (s_type)this->x0);
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
		}

		inner_type_m_with_bool sin2r()  const {

			if (!this->is_zero_vector && !is_zero_matrix) {
				return inner_type_m_with_bool(false, std::cos((s_type)this->x0) * x2 - std::sin((s_type)this->x0) * viennacl::linalg::outer_prod(this->x1, this->x1));
			}
			else if (!this->is_zero_vector && is_zero_matrix) {
				return inner_type_m_with_bool(false, -std::sin((s_type)this->x0) * viennacl::linalg::outer_prod(this->x1, this->x1));
			}
			else if (this->is_zero_vector && !is_zero_matrix) {
				return inner_type_m_with_bool(false, std::cos((s_type)this->x0) * x2);
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
		}

		inner_type_m_with_bool cos2r()  const {

			if (!this->is_zero_vector && !is_zero_matrix) {
				return inner_type_m_with_bool(false, -std::sin((s_type)this->x0) * x2 - std::cos((s_type)this->x0) * viennacl::linalg::outer_prod(this->x1, this->x1));
			}
			else if (!this->is_zero_vector && is_zero_matrix) {
				return inner_type_m_with_bool(false, -std::cos((s_type)this->x0) * viennacl::linalg::outer_prod(this->x1, this->x1));
			}
			else if (this->is_zero_vector && !is_zero_matrix) {
				return inner_type_m_with_bool(false, -std::sin((s_type)this->x0) * x2);
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
		}

		inner_type_m_with_bool tan2r()  const {

			if (!this->is_zero_vector && !is_zero_matrix) {
				const s_type tanx = std::tan((s_type)this->x0);
				const s_type c = static_cast<s_type>(1) + tanx * tanx;
				return inner_type_m_with_bool(false, c * (x2 + static_cast<s_type>(2) * tanx * viennacl::linalg::outer_prod(this->x1, this->x1)) );
			}
			else if (!this->is_zero_vector && is_zero_matrix) {
				const s_type tanx = std::tan((s_type)this->x0);
				const s_type c = static_cast<s_type>(1) + tanx * tanx;
				return inner_type_m_with_bool(false, c * static_cast<s_type>(2) * tanx * viennacl::linalg::outer_prod(this->x1, this->x1) );
			}
			else if (this->is_zero_vector && !is_zero_matrix) {
				const s_type tanx = std::tan((s_type)this->x0);
				const s_type c = static_cast<s_type>(1) + tanx * tanx;
				return inner_type_m_with_bool(false, c * x2);
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
		}

		inner_type_m_with_bool asin2r()  const {

			if (!this->is_zero_vector && !is_zero_matrix) {
				const s_type denom_buf = std::sqrt(static_cast<s_type>(1) - std::pow((s_type)this->x0, 2));
				return inner_type_m_with_bool(false, x2 / denom_buf + this->x0 * viennacl::linalg::outer_prod(this->x1, this->x1) / std::pow(denom_buf, 3) );
			}
			else if (!this->is_zero_vector && is_zero_matrix) {
				const s_type denom_buf = std::sqrt(static_cast<s_type>(1) - std::pow((s_type)this->x0, 2));
				return inner_type_m_with_bool(false, this->x0 * viennacl::linalg::outer_prod(this->x1, this->x1) / std::pow(denom_buf, 3) );
			}
			else if (this->is_zero_vector && !is_zero_matrix) {
				const s_type denom_buf = std::sqrt(static_cast<s_type>(1) - std::pow((s_type)this->x0, 2));
				return inner_type_m_with_bool(false, x2 / denom_buf);
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
		}

		inner_type_m_with_bool acos2r()  const {

			if (!this->is_zero_vector && !is_zero_matrix) {
				const s_type denom_buf = std::sqrt(static_cast<s_type>(1) - std::pow((s_type)this->x0, 2));
				return inner_type_m_with_bool(false, x2 / (-denom_buf) - (s_type)this->x0 * viennacl::linalg::outer_prod(this->x1, this->x1) / std::pow(denom_buf, 3));
			}
			else if (!this->is_zero_vector && is_zero_matrix) {
				const s_type denom_buf = std::sqrt(static_cast<s_type>(1) - std::pow((s_type)this->x0, 2));
				return inner_type_m_with_bool(false, -(s_type)this->x0 * viennacl::linalg::outer_prod(this->x1, this->x1) / std::pow(denom_buf, 3));
			}
			else if (this->is_zero_vector && !is_zero_matrix) {
				const s_type denom_buf = std::sqrt(static_cast<s_type>(1) - std::pow((s_type)this->x0, 2));
				return inner_type_m_with_bool(false, x2 / (-denom_buf));
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
		}

		inner_type_m_with_bool atan2r()  const {

			if (!this->is_zero_vector && !is_zero_matrix) {
				const s_type denom_buf = static_cast<s_type>(1) + std::pow((s_type)this->x0, 2);
				return inner_type_m_with_bool(false, x2 / denom_buf - static_cast<s_type>(2) * (s_type)this->x0 * viennacl::linalg::outer_prod(this->x1, this->x1) / std::pow(denom_buf, 2));
			}
			else if (!this->is_zero_vector && is_zero_matrix) {
				const s_type denom_buf = static_cast<s_type>(1) + std::pow((s_type)this->x0, 2);
				return inner_type_m_with_bool(false, -static_cast<s_type>(2) * (s_type)this->x0 * viennacl::linalg::outer_prod(this->x1, this->x1) / std::pow(denom_buf, 2));
			}
			else if (this->is_zero_vector && !is_zero_matrix) {
				const s_type denom_buf = static_cast<s_type>(1) + std::pow((s_type)this->x0, 2);
				return inner_type_m_with_bool(false, x2 / denom_buf);
			}
			else {
				return inner_type_m_with_bool(true, inner_type_m());
			}
		}
	};

}


