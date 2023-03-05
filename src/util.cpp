#include <vector>
#include <numeric>
#include <cstdio>
#include <cstdlib>

size_t factorial(size_t k) {
	if (k <= 1) {
		return 1;
	} else {
		return k * factorial(k - 1);
	}
}

std::vector<std::vector<int>> nchoosek(int n, int k) {
	std::vector<std::vector<int>> rc;
	const int N = factorial(n) / factorial(n - k) / factorial(k);
	std::vector<int> combo(k);
	std::iota(combo.begin(), combo.end(), 0);
	if (n == k) {
		rc.push_back(combo);
	} else {
		for (int q = 0; q < N - 1; q++) {
			rc.push_back(combo);
			for (int i = 0; i < combo.size(); i++) {
				printf("%i ", combo[i]);
			}
			printf("\n");
			int dim = k - 1;
			while (combo[dim] == n - 1) {
				dim--;
			}
			REDO: combo[dim]++;
			for (int i = dim + 1; i < combo.size(); i++) {
				combo[i] = combo[i - 1] + 1;
			}
			if (combo.back() >= n) {
				dim--;
				goto REDO;
			}

		}
	}
	rc.push_back(combo);
	return rc;
}
