#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>  
#include <igl/sort.h>  
#include <igl/pinv.h>  
#include <igl/per_vertex_normals.h>  
#include <igl/sort.h>
 
void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{ 
    K1 = Eigen::VectorXd::Zero(V.rows());
    K2 = Eigen::VectorXd::Zero(V.rows());
    D1 = Eigen::MatrixXd::Zero(V.rows(),3);
    D2 = Eigen::MatrixXd::Zero(V.rows(),3); 
    // Replace with your code
    Eigen::SparseMatrix<double> A;
    igl::adjacency_matrix(F, A);
    A.makeCompressed(); 

    for (int i = 0; i < V.rows(); ++i)
    {
        //1. calculate second neighbors and P
        Eigen::VectorXd judge = Eigen::VectorXd::Zero(V.rows());
        Eigen::VectorXd secNgbs = Eigen::VectorXd::Zero(V.rows());
        Eigen::MatrixXd P = Eigen::MatrixXd::Zero(V.rows(), 3);
        int nums = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it)
        {
            int firNeighbor = it.row();   // row index   
            for (Eigen::SparseMatrix<double>::InnerIterator its(A, firNeighbor); its; ++its)
            {
                int secNeighbor = its.row();   // row index   
                if (!judge(secNeighbor) && secNeighbor != i)
                {
                    judge(secNeighbor) = 1;
                    secNgbs(nums) = secNeighbor;
                    P.row(nums) = V.row(secNeighbor) - V.row(i);  //v_i - v
                    nums++;
                }
            }
        } 
        secNgbs.conservativeResize(nums, 1);
        P.conservativeResize(nums, 3);  
        //2. eigen decomposition
        Eigen::JacobiSVD<Eigen::MatrixXd> svd;
        svd.compute(P, Eigen::ComputeThinU | Eigen::ComputeThinV); 
        Eigen::MatrixXd matrixV = svd.matrixV();
        Eigen::MatrixXd PP = matrixV.leftCols(2);
        Eigen::VectorXd minVec = matrixV.col(2); 

        Eigen::MatrixXd S = P * PP;
        Eigen::VectorXd B = P * minVec; 

        //3. least-squares fitting 
        Eigen::MatrixXd AA(P.rows(), 5);
        for (int j = 0; j < P.rows(); ++j)
        {
            AA(j, 0) = S(j, 0);
            AA(j, 1) = S(j, 1);
            AA(j, 2) = S(j, 0) * S(j, 0);
            AA(j, 3) = S(j, 0) * S(j, 1);
            AA(j, 4) = S(j, 1) * S(j, 1);
        }

        Eigen::MatrixXd AA_INV;
        igl::pinv(AA, AA_INV);
        Eigen::VectorXd X = AA_INV * B; 

        //4. The shape operator
        Eigen::MatrixXd firFunForm(2, 2), secFunForm(2, 2);
        firFunForm(0, 0) = 1 + X(0)*X(0);
        firFunForm(0, 1) = X(0) * X(1);
        firFunForm(1, 0) = firFunForm(0, 1);
        firFunForm(1, 1) = 1 + X(1)*X(1);
        //
        double denominator = std::sqrt(X(0)*X(0) + 1 + X(1)*X(1));
        secFunForm(0, 0) = 2 * X(2) / denominator;
        secFunForm(0, 1) = X(3) / denominator;
        secFunForm(1, 0) = secFunForm(0, 1);
        secFunForm(1, 1) = 2 * X(4) / denominator;

        Eigen::MatrixXd SS = - secFunForm * firFunForm.inverse(); 

        //5. the principal curvatures
        Eigen::EigenSolver<Eigen::MatrixXd> kSolve(SS); 
        Eigen::VectorXd kValues = kSolve.eigenvalues().real(); 
        Eigen::MatrixXd kVectors = kSolve.eigenvectors().real(); 

        //6. Lift the principal tangent directions back
        if (kValues(0) > kValues(1))
        {
            K1(i) = kValues(0);
            D1.row(i) = kVectors.col(0).transpose() * PP.transpose();

            K2(i) = kValues(1);
            D2.row(i) = kVectors.col(1).transpose() * PP.transpose();
        }
        else{
            K1(i) = kValues(1);
            D1.row(i) = kVectors.col(1).transpose() * PP.transpose();

            K2(i) = kValues(0);
            D2.row(i) = kVectors.col(0).transpose() * PP.transpose();
        } 
    }
}
