# HTN: Hyper Triple Numbers Library
 Automatic forward-mode differential library.  
 We can calculate gradient vector and hessian matrix automatically.  
 This library uses [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [ViennaCL](http://viennacl.sourceforge.net/doc/index.html) library.

 The library is MIT Lisence but Lisences of Eigen and ViennaCL library take over each one.  

 HTN's calculation is based on maltivariate Taylor expansion version of dual number.

 ## Usage Code

```cpp:main
//#define VIENNACL_WITH_OPENMP
//#define VIENNACL_WITH_OPENCL
//#define VIENNACL_WITH_CUDA
#include "HTN.hpp"

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/StdVector"
#include "Eigen/SVD"
#include "Eigen/Dense"
#include "Eigen/LU"

#include <iostream> 

void main() {

	using f_type = double;
	const static unsigned int dim = 3;

	/* scalar */
	using htn0 = HTN::htn_r0<f_type, dim>;

	/* scalar & gradient vector */
	using htn1 = HTN::htn_r1<f_type, dim>;

	/* scalar & gradient vector & hessian matrix */
	using htn2 = HTN::htn_r2<f_type, dim>;

	/* x_0 <= 3 */
	auto x0 = htn2::x(3); /* or htn2::x<0>(3); or htn2::x(3, 0); */
	/* x_1 <= -pi */
	auto x1 = htn2::x<1>(-HTN::PI()); /* or htn2::x(-HTN::PI(), 1) */

	/* htn is initialized to 0(constant). */
	/* x2 <= 0 (not variable) */
	htn2 x2; 

	/* !warning: One variable can represent multiple states.
	   example; auto x0 = htn2::x(3); and x0 = htn2::x(-2);  */

	/* function: 
	inverse, pow, sqrt, exp, log, 
	sin, cos, tan, asin, acos, atan, 
	sinh, cosh, tanh, asinh, acosh, atanh  */
	auto y = HTN::exp(htn2::x<1>(0.5));
	auto z = htn2::x<2>(0.5).acos();
	float a = 3.2;
	auto f = HTN::pow(2.5 * htn2::x(3), a) / y;
	f -= z;

	/* convert to scalar, Eigen vector(Dynamic size) and Eigen matrix(Dynamic size) */
	std::cout << f.X0_to_type() << std::endl << std::endl;
	std::cout << f.X1_to_eig() << std::endl << std::endl;
	std::cout << f.X2_to_eig() << std::endl << std::endl;

	/*-----------------------------------------------------------------------------*/


	auto Rosenbrock = [](const auto& data) {
		htn2 f;
		for (int i = 0; i <= htn2::get_dim() - 2; ++i) {
			f += 100 * HTN::pow(htn2::x(data[i + 1], i + 1) - htn2::x(data[i], i).pow(2), 2) + HTN::pow(1 - htn2::x(data[i], i), 2);
		}
		return f;
	};

	/*Newton Method*/
	std::cout << "---Minimization of Rosenbrock function by Newton method---" << std::endl;

	Eigen::Matrix<f_type, dim, 1> x(4, 3, -4);
	for (int i = 0; i < 15; ++i) {
		auto rb = Rosenbrock(x);
		x -= rb.X2_to_eig().inverse() * rb.X1_to_eig();

		std::cout << x.transpose() << std::endl;
	}
	std::cout << std::endl;


	std::cout << std::endl;
}

```
