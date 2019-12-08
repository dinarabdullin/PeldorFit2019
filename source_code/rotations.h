#ifndef ROTATIONS_H
#define ROTATIONS_H

#include <cmath>
#include <vector>

class RotationMatrix
{
public:
	std::vector<std::vector<double>> R;
	std::vector<std::vector<double>> RT;

	// Initialize as a unit matrix
	RotationMatrix(const double& alpha, const double& beta, const double& gamma) {

		double c1 = cos(alpha);
		double s1 = sin(alpha);
		double c2 = cos(beta);
		double s2 = sin(beta);
		double c3 = cos(gamma);
		double s3 = sin(gamma);
		// Rotation matrix(active rotation)
		R.reserve(3);
		std::vector<double> row; row.reserve(3);
		row.push_back( c1*c3 - c2*s1*s3 );
		row.push_back(-c1*s3 - c2*c3*s1 );
		row.push_back( s1*s2            );
		R.push_back(row);
		row.clear();
		row.push_back( c3*s1 + c1*c2*s3 );
		row.push_back( c1*c2*c3 - s1*s3 );
		row.push_back(-c1*s2            );
		R.push_back(row);
		row.clear();
		row.push_back( s2*s3            );
		row.push_back( c3*s2            );
		row.push_back( c2               );
		R.push_back(row);
		row.clear();
		// Transposed rotation matrix (passive rotation)
		RT.reserve(3);
		row.push_back(R[0][0]);
		row.push_back(R[1][0]);
		row.push_back(R[2][0]);
		RT.push_back(row);
		row.clear();
		row.push_back(R[0][1]);
		row.push_back(R[1][1]);
		row.push_back(R[2][1]);
		RT.push_back(row);
		row.clear();
		row.push_back(R[0][2]);
		row.push_back(R[1][2]);
		row.push_back(R[2][2]);
		RT.push_back(row);
		row.clear();
	}

	// Returns the rotation matrix for the Euler 'z-x-z' active rotation
	void reset_angles(const double& alpha, const double& beta, const double& gamma) {

		double c1 = cos(alpha);
		double s1 = sin(alpha);
		double c2 = cos(beta);
		double s2 = sin(beta);
		double c3 = cos(gamma);
		double s3 = sin(gamma);
		// Rotation matrix (active rotation)
		R[0][0] = c1*c3 - c2*s1*s3;
		R[0][1] = -c1*s3 - c2*c3*s1;
		R[0][2] = s1*s2;
		R[1][0] = c3*s1 + c1*c2*s3;
		R[1][1] = c1*c2*c3 - s1*s3;
		R[1][2] = -c1*s2;
		R[2][0] = s2*s3;
		R[2][1] = c3*s2;
		R[2][2] = c2;
		// Transposed rotation matrix (passive rotation)
		RT[0][0] = R[0][0];
		RT[0][1] = R[1][0];
		RT[0][2] = R[2][0];
		RT[1][0] = R[0][1];
		RT[1][1] = R[1][1];
		RT[1][2] = R[2][1];
		RT[2][0] = R[0][2];
		RT[2][1] = R[1][2];
		RT[2][2] = R[2][2];
	}

	// A dot product of the rotation matrix and a vector 
	std::vector<double> dot_product(std::vector<double> const& v, bool const& transposed) const 
	{
		std::vector<double> p; p.reserve(3);
		for (size_t i = 0; i < 3; i++) {
			if (transposed) {
				p.push_back(RT[i][0] * v[0] + RT[i][1] * v[1] + RT[i][2] * v[2]);
			}
			else {
				p.push_back(R[i][0] * v[0] + R[i][1] * v[1] + R[i][2] * v[2]);
			}
		}
		return p;
	}

	// A product of the rotation matrix and a 3x3 matrix 
	std::vector<std::vector<double>> matrix_product(std::vector<std::vector<double>> const& M, bool const& transposed) const
	{
		std::vector<std::vector<double>> P; P.reserve(3);
		std::vector<double> row; row.reserve(3);
		for (size_t i = 0; i < 3; ++i) {
			for (size_t j = 0; j < 3; ++j) {
				if (transposed) {
					row.push_back(RT[i][0] * M[0][j] + RT[i][1] * M[1][j] + RT[i][2] * M[2][j]);
				}
				else {
					row.push_back(R[i][0] * M[0][j] + R[i][1] * M[1][j] + R[i][2] * M[2][j]);
				}
			}
			P.push_back(row);
			row.clear();
		}
		return P;
	}
};

#endif