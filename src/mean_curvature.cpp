#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h> 
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>  
#include <igl/per_vertex_normals.h>  

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
    // Replace with your code
 	H = Eigen::VectorXd::Zero(V.rows());
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);

	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

    Eigen::SparseMatrix<double> M_INV;
	igl::invert_diag(M, M_INV);

	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);

	Eigen::MatrixXd Hn = M_INV * L * V;
 	for (int i = 0; i < V.rows(); ++i)
 	{
 		H(i) = Hn.row(i).norm() * (Hn(i,0)/N(i,0) > 0 ? 1 : -1);
 		// std::cout<<"n:"<<N.row(i)<<", Hn:"<<Hn.row(i)<<", Hi:"<<H(i)<<", sign:"<<Hn(i,0)/N(i,0)<<std::endl;
 	}
}
