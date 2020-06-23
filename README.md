# HTN: Hyper Triple Numbers Library
 Automatic forward-mode differential library.  
 We can calculate gradient vector and hessian matrix automatically.  
 This library uses [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [ViennaCL](http://viennacl.sourceforge.net/doc/index.html) library.

 The library is MIT Lisence but Lisences of Eigen and ViennaCL library take over each one.  

 HTN's calculation is based on maltivariate Taylor expansion version of dual number.

 ## Example Code

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

## Base Math

#### values
All HTN values are as follows;

<img src="https://latex.codecogs.com/gif.latex?\chi&space;=&space;x&space;&plus;&space;\mathbf{\epsilon}^{\rm&space;T}\mathbf{x}&space;&plus;&space;\frac{1}{2}\mathbf{\epsilon}^{\rm&space;T}&space;X&space;\mathbf{\epsilon}" />  

, where 

<img src="https://latex.codecogs.com/gif.latex?\mathbf{\epsilon}&space;=&space;(\epsilon_1~\epsilon_2~\cdots~\epsilon_D)^{\rm&space;T},~{\bf&space;x}&space;\in&space;\mathbb{R}^{D&space;\times&space;1},~X&space;\in&space;\mathbb{R}^{D&space;\times&space;D},~X&space;=&space;X^{\rm&space;T}" title="\mathbf{\epsilon} = (\epsilon_1~\epsilon_2~\cdots~\epsilon_D)^{\rm T},~{\bf x} \in \mathbb{R}^{D \times 1},~X \in \mathbb{R}^{D \times D},~X = X^{\rm T}" />

and dual units;

<img src="https://latex.codecogs.com/gif.latex?\inline&space;\epsilon_i&space;\neq&space;0,&space;~i&space;\neq&space;j&space;\rightarrow&space;\epsilon_i&space;\neq&space;\epsilon_j,&space;~\epsilon_i&space;\epsilon_j&space;=&space;\epsilon_j&space;\epsilon_i,&space;~\epsilon_i&space;\epsilon_j&space;\epsilon_k&space;=&space;0" title="\epsilon_i \neq 0, ~i \neq j \rightarrow \epsilon_i \neq \epsilon_j, ~\epsilon_i \epsilon_j = \epsilon_j \epsilon_i, ~\epsilon_i \epsilon_j \epsilon_k = 0" />
.

If HTN value is a i-th base variable, 
<img src="https://latex.codecogs.com/gif.latex?\inline&space;\chi&space;=&space;x_i&space;&plus;&space;\epsilon_i" title="\chi_i = x_i + \epsilon_i" />.

Else, if constant, <img src="https://latex.codecogs.com/gif.latex?\inline&space;\chi&space;=&space;c" title="\chi = c" />.

#### functions

All functions are based on  

<img src="https://latex.codecogs.com/gif.latex?f(\mathbf{\chi}=(\chi_1,~\chi_2,~...))&space;=&space;f&space;&plus;&space;\mathbf{\epsilon}^{\rm&space;T}&space;\sum_{i}&space;{\bf&space;x}_i&space;f_{x_i}&space;&plus;&space;\frac{1}{2}&space;\mathbf{\epsilon}^{\rm&space;T}&space;\left(&space;\sum_i&space;X_i&space;f_{x_i}&space;&plus;&space;\sum_i&space;\sum_j&space;{\bf&space;x}_i&space;{\bf&space;x}_j^{\rm&space;T}&space;f_{x_i&space;x_j}&space;\right)&space;\mathbf{\epsilon}" title="f(\mathbf{\chi}=(\chi_1,~\chi_2,~...)) = f + \mathbf{\epsilon}^{\rm T} \sum_{i} {\bf x}_i f_{x_i} + \frac{1}{2} \mathbf{\epsilon}^{\rm T} \left( \sum_i X_i f_{x_i} + \sum_i \sum_j {\bf x}_i {\bf x}_j^{\rm T} f_{x_i x_j} \right) \mathbf{\epsilon}" />  

, where <img src="https://latex.codecogs.com/gif.latex?\inline&space;f&space;=&space;f(x_1,&space;\cdots,&space;x_D)" title="f = f(x_1, \cdots, x_D)" />.

In particular,  

<img src="https://latex.codecogs.com/gif.latex?f(&space;\chi&space;)&space;=&space;f&space;&plus;&space;\mathbf{\epsilon}^{\rm&space;T}&space;{\bf&space;x}&space;f'&space;&plus;&space;\frac{1}{2}&space;\mathbf{\epsilon}^{\rm&space;T}&space;\left(&space;X&space;f'&space;&plus;&space;{\bf&space;x}{\bf&space;x}^{\rm&space;T}&space;f''&space;\right)&space;\mathbf{\epsilon}" title="f( \chi ) = f + \mathbf{\epsilon}^{\rm T} {\bf x} f' + \frac{1}{2} \mathbf{\epsilon}^{\rm T} \left( X f' + {\bf x}{\bf x}^{\rm T} f'' \right) \mathbf{\epsilon}" />  
<img src="https://latex.codecogs.com/gif.latex?f(\chi_1,&space;\chi_2)&space;=&space;f&space;&plus;&space;\mathbf{\epsilon}^{\rm&space;T}&space;({\bf&space;x}_1&space;f_{x_1}&space;&plus;&space;{\bf&space;x}_2&space;f_{x_2})&space;&plus;&space;\frac{1}{2}&space;\mathbf{\epsilon}^{\rm&space;T}&space;\left(&space;X_1&space;f_{x_1}&space;&plus;&space;X_2&space;f_{x_2}&space;&plus;&space;{\bf&space;x}_1&space;{\bf&space;x}_1^{\rm&space;T}&space;f_{x_1&space;x_1}&space;&plus;&space;({\bf&space;x}_1&space;{\bf&space;x}_2^{\rm&space;T}&space;&plus;&space;{\bf&space;x}_2&space;{\bf&space;x}_1^{\rm&space;T})&space;f_{x_1&space;x_2}&space;&plus;&space;{\bf&space;x}_2&space;{\bf&space;x}_2^{\rm&space;T}&space;f_{x_2&space;x_2}&space;\right)&space;\mathbf{\epsilon}" title="f(\chi_1, \chi_2) = f + \mathbf{\epsilon}^{\rm T} ({\bf x}_1 f_{x_1} + {\bf x}_2 f_{x_2}) + \frac{1}{2} \mathbf{\epsilon}^{\rm T} \left( X_1 f_{x_1} + X_2 f_{x_2} + {\bf x}_1 {\bf x}_1^{\rm T} f_{x_1 x_1} + ({\bf x}_1 {\bf x}_2^{\rm T} + {\bf x}_2 {\bf x}_1^{\rm T}) f_{x_1 x_2} + {\bf x}_2 {\bf x}_2^{\rm T} f_{x_2 x_2} \right) \mathbf{\epsilon}" />
.
