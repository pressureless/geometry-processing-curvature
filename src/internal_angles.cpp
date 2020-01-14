#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
    // Add with your code
	A.resize(l_sqr.rows(), l_sqr.cols());
	for (int i = 0; i < l_sqr.rows(); ++i)
	{
		double cos0 = (l_sqr(i,1)+l_sqr(i,2)-l_sqr(i,0))/(2*std::sqrt(l_sqr(i,1)*l_sqr(i,2)));
		A(i,0) = std::acos(cos0);

		double cos1 = (l_sqr(i,0)+l_sqr(i,2)-l_sqr(i,1))/(2*std::sqrt(l_sqr(i,0)*l_sqr(i,2)));
		A(i,1) = std::acos(cos1);

		double cos2 = (l_sqr(i,1)+l_sqr(i,0)-l_sqr(i,2))/(2*std::sqrt(l_sqr(i,1)*l_sqr(i,0)));
		A(i,2) = std::acos(cos2);
	}
}
