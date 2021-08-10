#include <vector>
#include <cmath>

template <class T>
void eigenJacobiSolve(const std::vector<std::vector<T>>& input,
	std::vector<T>& eigenvalues,
	std::vector<std::vector<T>>& eigenvectors)
{
	int n = static_cast<int>(input.size());
	if (n <= 1) return;
	for (auto& vec : input)
		if (static_cast<int>(vec.size()) != n) return;
	
	//initi
	std::vector<std::vector<T>> S(input);
	std::vector<std::vector<T>> V(n, std::vector<T>(n, 0));
	const double eps = 1e-7;
	for (int i = 0; i < n; ++i)
		V[i][i] = (T)1;

	int iter = 0;
	const int maxIter = n * n * 30;
	while (1) {
		//search max Sij
		int i = 0;
		int j = 1;
		T Sij = S[i][j];
		T Sii, Sjj;

		for (int r = 0; r < n; ++r)
			for (int c = 0; c < n; ++c)
				if (r != c && fabs(S[r][c]) > fabs(Sij)) {
					Sij = S[r][c];
					i = r;
					j = c;
				}
		Sii = S[i][i];
		Sjj = S[j][j];
		iter++;

		if (fabs(Sij) < eps || iter > maxIter)
			break;

		double theta = atan2(-2.0 * Sij, Sjj - Sii) * 0.5;
		double cos_theta = cos(theta);
		double sin_theta = sin(theta);

		//update
		S[i][i] = static_cast<T>(cos_theta * cos_theta * Sii + 2.0 * cos_theta * sin_theta * Sij + sin_theta * sin_theta * Sjj);
		S[j][j] = static_cast<T>(sin_theta * sin_theta * Sii - 2.0 * cos_theta * sin_theta * Sij + cos_theta * cos_theta * Sjj);
		S[i][j] = static_cast<T>((cos_theta * cos_theta - sin_theta * sin_theta) * Sij - (Sii - Sjj) * sin_theta * cos_theta);
		S[j][i] = S[i][j];

		for (int k = 0; k < n; ++k) {
			if (k != i && k != j) {
				T Sik = S[i][k];
				T Sjk = S[j][k];
				S[i][k] = static_cast<T>(cos_theta * Sik + sin_theta * Sjk);
				S[k][i] = S[i][k];
				S[j][k] = static_cast<T>(-sin_theta * Sik + cos_theta * Sjk);
				S[k][j] = S[j][k];
			}
		}

		for (int k = 0; k < n; ++k) {
			T Ski = V[k][i];
			V[k][i] = static_cast<T>(V[k][j] * sin_theta + Ski * cos_theta);
			V[k][j] = static_cast<T>(V[k][j] * cos_theta - Ski * sin_theta);
		}
	}

	eigenvalues.resize(n);
	eigenvectors.resize(n, std::vector<T>(n));
	for (int k = 0; k < n; ++k)
		eigenvalues[k] = S[k][k];

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			eigenvectors[i][j] = V[j][i];
}