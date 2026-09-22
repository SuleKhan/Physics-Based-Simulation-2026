#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <memory>

template <typename T, int dim>
using Vector = Eigen::Matrix<T, dim, 1, 0, dim, 1>;

template <typename T, int n, int m>
using Matrix = Eigen::Matrix<T, n, m, 0, n, m>;

using Vector2F = Vector<double, 2>;
using Vector3F = Vector<double, 3>;
using Vector4F = Vector<double, 4>;
using Vector2I = Vector<int, 2>;
using Vector3I = Vector<int, 3>;
using Vector4I = Vector<int, 4>;

using Matrix2F = Matrix<double, 2, 2>;
using Matrix3F = Matrix<double, 3, 3>;

using VectorXF = Matrix<double, Eigen::Dynamic, 1>;
using MatrixXF = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXI = Vector<int, Eigen::Dynamic>;
using MatrixXI = Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;

using TripletF = Eigen::Triplet<double>;
using TripletListF = std::vector<TripletF>;
using SparseVectorF = Eigen::SparseVector<double>;
using SparseMatrixF = Eigen::SparseMatrix<double>;
using TripletI = Eigen::Triplet<int>;
using TripletListI = std::vector<TripletI>;
using SparseVectorI = Eigen::SparseVector<int>;
using SparseMatrixI = Eigen::SparseMatrix<int>;
