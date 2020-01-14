#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include <igl/squared_edge_lengths.h>  


void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
    D = Eigen::VectorXd::Zero(V.rows());
	  Eigen::MatrixXd L, A;
	  igl::squared_edge_lengths(V, F, L);

	  internal_angles(L, A);

	  Eigen::VectorXd S = Eigen::VectorXd::Zero(V.rows());
    for (int i = 0; i < F.rows(); ++i)
    {
      S(F(i,0)) += A(i,0);
      S(F(i,1)) += A(i,1);
      S(F(i,2)) += A(i,2);
    }
    
    for (int i = 0; i < V.rows(); ++i)
    {
    	D(i) = 2*M_PI - S(i);
    }
}
